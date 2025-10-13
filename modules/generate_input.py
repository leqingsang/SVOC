import re, pandas as pd
import gzip

def get_canonical_transcript(gene, exon_number, gtf_file_path):
    """
    Retrieve RefSeq values for the specified gene name and exon number(s) from the GTF file.
    
    Parameters:
        gene (str): Gene name, such as "BRCA1".
        exon_number (int or list): Exon number(s), such as 9 or [6, 7].
        gtf_file_path (str): The path to the GTF file.
    
    Returns:
        list: A list of RefSeq values. If no RefSeq value is found, returns an empty list.
    """
    # Compile regular expressions to match gene names and RefSeq
    gene_pattern = re.compile(r'gene\s+"' + re.escape(gene) + r'"')
    refseq_pattern = re.compile(r'transcript_id\s+"(NM_\d+)')
    tag_pattern = re.compile(r'tag\s+"([^"]+)"')
    
    # Initialize a list to store all matched RefSeq values
    refseq_tags = {}

    # Handle exon_number as either an integer or a list
    if isinstance(exon_number, int):
        # If exon_number is an integer, compile a single pattern
        exon_patterns = [re.compile(r'exon_number\s+"' + str(exon_number) + r'"')]
    elif isinstance(exon_number, list):
        # If exon_number is a list, compile a pattern for each exon number
        exon_patterns = [re.compile(r'exon_number\s+"' + str(num) + r'"') for num in exon_number]
    else:
        raise ValueError("exon_number must be an integer or a list of integers.")

    try:
        with gzip.open(gtf_file_path, 'rt') as file:
            for line in file:
                # Check if the gene name is included
                if gene_pattern.search(line):
                    # Check each exon pattern
                    for pattern in exon_patterns:
                        if pattern.search(line):
                            match_refseq = refseq_pattern.search(line)
                            match_tag = tag_pattern.search(line)
                            if match_refseq and match_tag:
                                refseq = match_refseq.group(1)
                                tag = match_tag.group(1)
                                refseq_tags[refseq] = tag
    except FileNotFoundError:
        print(f"Error: The file {gtf_file_path} was not found.")
    except Exception as e:
        print(f"Error: {e}")

    return refseq_tags


def transvar_transcript_filter(input_file, output_file, gtf_file_path, mark_file):
    """
    Filter transcripts, but keep all transcripts for a genomic coordinate if none of them are canonical.
    Also, mark each line and output a marked file.
    """
    # First pass: Identify genomic coordinates with canonical transcripts
    coordinate_canonical_status = {}
    with open(input_file, "r") as infile:
        headers = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info']
        #infile.readline()  # Skip header
        print("Start scanning transcripts...")
        for line in infile:
            fields = line.strip().split('\t')
            row = dict(zip(headers, fields))
            
            # Skip intergenic variations or rows without transcript information
            if "Intergenic" in row['info'] or row['transcript'] == '.':
                continue
            
            # Extract genomic coordinate identifier
            coordinate_id = f"{row['vChr']}_{row['vPos']}_{row['vrsID']}_{row['vRef']}_{row['vAlt']}"
            match1 = re.search(r'\w+_\d+', row['transcript'])
            if match1:
                transcript = match1.group(0)
            else:
                continue  # 或者根据需要设置一个默认值
            # Extract transcript and region information
            #transcript = re.search(r'\w+_\d+', row['transcript']).group(0)
            match = re.search(r'_exon_(\d+)|_exons_\[(\d+(?:,\d+)*)\]', row['region'])
            if match:
                exon_number = match.group(1) or match.group(2).split(',')
                exon_number = [int(num) for num in exon_number]
            else:
                exon_number = None
            
            # Get canonical transcripts for the gene and exon
            canonical_transcripts = get_canonical_transcript(row['gene'], exon_number, gtf_file_path)
            
            # Update the coordinate's canonical status
            if coordinate_id not in coordinate_canonical_status:
                coordinate_canonical_status[coordinate_id] = False
            if transcript in list(canonical_transcripts.keys()):
                coordinate_canonical_status[coordinate_id] = True
    
    # Second pass: Generate the mark file
    with open(input_file, "r") as infile, open(mark_file, "w") as markfile:
        headers = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info']
        headers_with_mark = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info', 'mark', 'transcript_tag']
        markfile.write('\t'.join(headers_with_mark) + '\n')
        infile.seek(0)  # Ensure we start from the beginning of the file
        #infile.readline()  # Skip header
        
        for line in infile:
            fields = line.strip().split('\t')
            row = dict(zip(headers, fields))
            
            # Skip intergenic variations or rows without transcript information
            if "Intergenic" in row['info'] or row['transcript'] == '.':
                mark = '1'  # Mark as 1 for intergenic or no transcript
                transcript_tag = '.'
                markfile.write(line.strip() + f'\t{mark}\t{transcript_tag}\n')
                continue
            
            # Extract genomic coordinate identifier
            coordinate_id = f"{row['vChr']}_{row['vPos']}_{row['vrsID']}_{row['vRef']}_{row['vAlt']}"
            
            # Extract transcript and region information
            match1 = re.search(r'\w+_\d+', row['transcript'])
            if match1:
                transcript = match1.group(0)
            else:
                continue  # 或者根据需要设置一个默认值
            match = re.search(r'_exon_(\d+)|_exons_\[(\d+(?:,\d+)*)\]', row['region'])
            if match:
                exon_number = match.group(1) or match.group(2).split(',')
                exon_number = [int(num) for num in exon_number]
            else:
                exon_number = None
            
            # Determine the mark for the current line
            if coordinate_canonical_status.get(coordinate_id, False):
                transcript_tags = get_canonical_transcript(row['gene'], exon_number, gtf_file_path)
                if transcript in list(transcript_tags.keys()):
                    mark = '2'  # Mark as 2 for canonical transcript
                    transcript_tag = transcript_tags[transcript]
                else:
                    mark = '3'  # Mark as 3 for non-canonical transcript in a coordinate with canonical transcripts
                    transcript_tag = '.'
            else:
                mark = '1'  # Mark as 1 for non-canonical transcript in a coordinate with no canonical transcripts
                transcript_tag = '.'
            
            markfile.write(line.strip() + f'\t{mark}\t{transcript_tag}\n')
        print("Transcript scan completed.")

    # Third pass: Read the mark file and write lines to the output file based on the marks
    with open(mark_file, "r") as markfile, open(output_file, "w") as outfile:
        headers_with_mark = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info', 'mark', 'transcript_tag']
        outfile.write('\t'.join(headers_with_mark) + '\n')
        markfile.readline()  # Skip header
        
        for line in markfile:
            fields = line.strip().split('\t')
            row = dict(zip(headers_with_mark, fields))
            mark = row['mark']
            
            if mark in ['1', '2']:
                outfile.write('\t'.join(fields) + '\n')
    
def transvar_transcript_filter_web(input_file, output_file, gtf_file_path, mark_file):
    """
    Filter transcripts, but keep all transcripts for a genomic coordinate if none of them are canonical.
    Also, mark each line and output a marked file.
    """
    # First pass: Identify genomic coordinates with canonical transcripts
    coordinate_canonical_status = {}
    with open(input_file, "r") as infile:
        headers = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info']
        #infile.readline()  # Skip header
        print("Start scanning transcripts...")
        for line in infile:
            fields = line.strip().split('\t')
            row = dict(zip(headers, fields))
            
            # Skip intergenic variations or rows without transcript information
            if "Intergenic" in row['info'] or row['transcript'] == '.':
                continue
            
            # Extract genomic coordinate identifier
            coordinate_id = f"{row['vChr']}_{row['vPos']}_{row['vrsID']}_{row['vRef']}_{row['vAlt']}"
            
            # Extract transcript and region information
            match1 = re.search(r'\w+_\d+', row['transcript'])
            if match1:
                transcript = match1.group(0)
            else:
                continue  # 或者根据需要设置一个默认值
            #transcript = re.search(r'\w+_\d+', row['transcript']).group(0)
            match = re.search(r'_exon_(\d+)|_exons_\[(\d+(?:,\d+)*)\]', row['region'])
            if match:
                exon_number = match.group(1) or match.group(2).split(',')
                exon_number = [int(num) for num in exon_number]
            else:
                exon_number = None
            
            # Get canonical transcripts for the gene and exon
            canonical_transcripts = get_canonical_transcript(row['gene'], exon_number, gtf_file_path)
            
            # Update the coordinate's canonical status
            if coordinate_id not in coordinate_canonical_status:
                coordinate_canonical_status[coordinate_id] = False
            if transcript in list(canonical_transcripts.keys()):
                coordinate_canonical_status[coordinate_id] = True
    
    # Second pass: Generate the mark file
    with open(input_file, "r") as infile, open(mark_file, "w") as markfile:
        headers = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info']
        headers_with_mark = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info', 'mark', 'transcript_tag']
        markfile.write('\t'.join(headers_with_mark) + '\n')
        infile.seek(0)  # Ensure we start from the beginning of the file
        #infile.readline()  # Skip header
        
        for line in infile:
            fields = line.strip().split('\t')
            row = dict(zip(headers, fields))
            
            # Skip intergenic variations or rows without transcript information
            if "Intergenic" in row['info'] or row['transcript'] == '.':
                mark = '1'  # Mark as 1 for intergenic or no transcript
                transcript_tag = '.'
                markfile.write(line.strip() + f'\t{mark}\t{transcript_tag}\n')
                continue
            
            # Extract genomic coordinate identifier
            coordinate_id = f"{row['vChr']}_{row['vPos']}_{row['vrsID']}_{row['vRef']}_{row['vAlt']}"
            
            # Extract transcript and region information
            match1 = re.search(r'\w+_\d+', row['transcript'])
            if match1:
                transcript = match1.group(0)
            else:
                continue  # 或者根据需要设置一个默认值
            match = re.search(r'_exon_(\d+)|_exons_\[(\d+(?:,\d+)*)\]', row['region'])
            if match:
                exon_number = match.group(1) or match.group(2).split(',')
                exon_number = [int(num) for num in exon_number]
            else:
                exon_number = None
            
            # Determine the mark for the current line
            if coordinate_canonical_status.get(coordinate_id, False):
                transcript_tags = get_canonical_transcript(row['gene'], exon_number, gtf_file_path)
                if transcript in list(transcript_tags.keys()):
                    mark = '2'  # Mark as 2 for canonical transcript
                    transcript_tag = transcript_tags[transcript]
                else:
                    mark = '3'  # Mark as 3 for non-canonical transcript in a coordinate with canonical transcripts
                    transcript_tag = '.'
            else:
                mark = '1'  # Mark as 1 for non-canonical transcript in a coordinate with no canonical transcripts
                transcript_tag = '.'
            
            markfile.write(line.strip() + f'\t{mark}\t{transcript_tag}\n')
        print("Transcript scan completed.")

    # Third pass: Read the mark file and write lines to the output file based on the marks
    with open(mark_file, "r") as markfile, open(output_file, "w") as outfile:
        headers_with_mark = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info', 'mark', 'transcript_tag']
        outfile.write('\t'.join(headers_with_mark) + '\n')
        markfile.readline()  # Skip header
        
        for line in markfile:
            fields = line.strip().split('\t')
            row = dict(zip(headers_with_mark, fields))
            mark = row['mark']
            
            if mark in ['2']:
                outfile.write('\t'.join(fields) + '\n')


def merge_transvar_annovar(transvar_path, annovar_path, output_path):
    try:
        transvar = pd.read_csv(transvar_path, sep='\t', dtype=str)
    except Exception as e:
        print(f"Error reading {transvar_path}: {e}")
        exit()
    try:
        annovar = pd.read_csv(annovar_path, sep='\t', dtype=str)
    except Exception as e:
        print(f"Error reading {annovar_path}: {e}")
        exit()
    
    transvar_keys = ['vChr', 'vPos', 'vrsID', 'vRef', 'vAlt']
    annovar_keys = ['Otherinfo4', 'Otherinfo5', 'Otherinfo6', 'Otherinfo7', 'Otherinfo8']

    # Ensure that both files have these fields
    if not all(key in transvar.columns for key in transvar_keys):
        print(f"Error: Missing columns in {transvar_path}. Required columns: {transvar_keys}")
        exit()
    if not all(key in annovar.columns for key in annovar_keys):
        print(f"Error: Missing columns in {annovar_path}. Required columns: {annovar_keys}")
        exit()

    # Rename the matching fields of the second file to merge with the first file
    annovar_renamed = annovar.rename(columns=dict(zip(annovar_keys, transvar_keys)))

    # Merge two files
    merged = pd.merge(transvar, annovar_renamed, on=transvar_keys, how='inner')
    merged.to_csv(output_path, sep='\t', index=False)
