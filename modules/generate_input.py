import re, pandas as pd

def get_canonical_transcript(gene, exon_number, gtf_file_path):
    """
    Retrieve RefSeq values for the specified gene name and exon number from the GTF file.
    
    Parameters:
		gene (str):  Gene name, such as "BRCA1".
		exon_number (int):  Exon number, such as 9.
		Gtf_file_cath (str): The path to the GTF file.
    
    return:
		str:  If the RefSeq value is not found, return None.
    """
    # Compile regular expressions to match gene names, exon numbers, and RefSeq
    gene_pattern = re.compile(r'gene_name\s+"' + re.escape(gene) + r'"')
    exon_pattern = re.compile(r'exon_number\s+' + str(exon_number) + r';')
    refseq_pattern = re.compile(r'db_xref\s+"RefSeq:(NM_\d+)')
    refseq_matches = []  # Store all matched RefSeq
    try:
        with open(gtf_file_path, 'r') as file:
            for line in file:
                # Check if the gene name and exon number are included
                if gene_pattern.search(line) and exon_pattern.search(line):
                    # Extract RefSeq
                    match = refseq_pattern.findall(line)
                    if match:
                        refseq_matches.extend(match)
    except FileNotFoundError:
        print(f"Error: The file {gtf_file_math} was not found.")
    except Exception as e:
        print(f"Error：{e}")
    return refseq_matches

def transvar_transcript_filter(input_file, output_file, gtf_file_path):
    """
    Filter transcripts.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        headers = ['vChr','vPos','vrsID','vRef','vAlt','QUAL','FILTER','INFO','transcript','gene','strand','coordinate','region','info']
        outfile.write('\t'.join(headers) + '\n')  # Separate headers with tabs and write them to a file
        for line in infile:
            fields = line.strip().split('\t')  # Split fields
            row = dict(zip(headers, fields)) # Use the zip function to combine key value pairs and convert them into dictionaries using the dict function
            gene = row['gene']  
            var_type = row['info'].split('CSQN=')[1].split(';')[0] if ';' in row['info'] else row['info'].split('CSQN=')[1]
            if var_type == "IntergenicSNV": # There is no transcript of intergenic variation, no filtering!!!
                outfile.write(line)
            else:
                transcript = re.search(r'\w+_\d+', row['transcript']).group(0) # group(0) Return the entire matching string.
                exon_number = re.search(r'_exon_(\d+)', row['region']).group(1) # group(1) Return the content that matches the first capture group (i.e., (\ d+)).
                canonical_transcripts = get_canonical_transcript(gene, exon_number, gtf_file_path)
                if transcript in canonical_transcripts:
                    outfile.write(line)

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