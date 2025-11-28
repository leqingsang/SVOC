import os, asyncio, re, getopt, sys, pandas as pd, textwrap, subprocess, copy, logging, io, time, platform, optparse, gzip, glob, json, configparser
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1' 
import modules.db_handler as db_handler
import modules.check_evidence as check_evidence
import modules.spliceai_prediction as spliceai
from modules.query_gnomad import query_gnomAD_local
from modules.generate_input import transvar_transcript_filter, merge_transvar_annovar
from modules.output_filter import svoc_output_filter

prog="SVOC"
version = """Version: 20250606
Written by Leqing SANG, lqsang25@m.fudan.edu.cn. 
SVOC is free for non-commercial use without warranty.
Please contact the author for commercial use.
Copyright (C) 2025 Ying YU Lab.
"""
usage = """Usage: ./%prog -c config.ini   # Run the examples in config.ini
       ./%prog  -i your_vcf_input -o prefix_of_your_output
"""
begin_description = """=============================================================================
Classification of the oncogenicity of somatic variants using python scripts. ^_^

 ########  ##    ##  ########  ########
 ##        ##    ##  ##    ##  ##      
 ##        ##    ##  ##    ##  ##      
 ########  ##    ##  ##    ##  ##      
       ##   ##  ##   ##    ##  ##      
       ##   ##  ##   ##    ##  ##      
 ########    ####    ########  ########

=============================================================================
"""
end_description = """
=============================================================================

 ########  ##    ##  ########  ########
 ##        ##    ##  ##    ##  ##      
 ##        ##    ##  ##    ##  ##      
 ########  ##    ##  ##    ##  ##      
       ##   ##  ##   ##    ##  ##      
       ##   ##  ##   ##    ##  ##      
 ########    ####    ########  ########

Thanks for using SVOC!
Report bugs to yqsang23@m.fudan.edu.cn.
=============================================================================
"""
paras = {}

def ConfigSectionMap(config,section):
    """
    Read the configuration items of a specified section from the configuration file and store these configuration items in a global variable paras.
    """
    global paras
    options = config.options(section)
    for option in options: # Traverse each option.
        try:
            paras[option] = config.get(section, option) # Attempt to read the value of the configuration item.
            if paras[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            paras[option] = None
    return

def getAAchangePart(AAchange_single):
    refAA = ''
    posAA = ''
    altAA = ''
    startAA = ''
    endAA = ''
    AAchange_position = ''
    if ('del' in AAchange_single or 'delins' in AAchange_single or 'ins' in AAchange_single or 'dup' in AAchange_single) and '_' in AAchange_single:
        startAA = re.match(r"p\.[A-Z](\d+)_[A-Z](\d+)", AAchange_single).group(1)
        endAA = re.search(r"p\.[A-Z](\d+)_[A-Z](\d+)", AAchange_single).group(2)
        posAA = startAA + '_' + endAA
        refAA = re.match(r"p\.([A-Z])\d+_([A-Z])\d+", AAchange_single).group(1) + re.match(r"p\.([A-Z])\d+_([A-Z])\d+", AAchange_single).group(2)
        altAA = ''
        AAchange_position = re.match(r"p\.[A-Z]\d+_[A-Z]\d+", AAchange_single).group(0)
    elif AAchange_single:
        refAA = re.search(r"p\.([A-Za-z]|\*)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(1)
        posAA = re.search(r"p\.([A-Za-z]|\*|)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(2)
        altAA = re.search(r"p\.([A-Za-z]|\*|)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(3)
        startAA = posAA
        endAA = posAA
        AAchange_position = 'p.'+ refAA + posAA
    return refAA, posAA, altAA, startAA, endAA, AAchange_position

def main():
    config=configparser.ConfigParser()
    parser = optparse.OptionParser(version=version, usage=usage)
    parser.add_option("-c", "--config", dest="config", action="store",
                  help="The config file of all options. It is for your own configure file.You can edit all the options in the configure and if you use this options,you can ignore all the other options bellow.", metavar="config.ini")
    parser.add_option("-b", "--buildver", dest="buildver", action="store",
                  help="The genomic build version, it can be hg19 , will support other version later.", metavar="hg19")
    parser.add_option("-i", "--input", dest="input", action="store",
                  help="The input file contains your variants.", metavar="example/example1.vcf")
    parser.add_option("-o", "--output", dest="output", action="store",
                  help="The prefix of output file which contains the results. ", metavar="example/example1")
    parser.add_option("--output_type", dest="output_type", action="store",
                  help="The output file type, it can be json, csv and txt.", metavar="txt")
    parser.add_option("-d", "--database_svoc", dest="database_svoc", action="store",
                  help="The  database location/dir for the SVOC dataset files.", metavar="svocdb")
    
    group = optparse.OptionGroup(parser, "Annovar Options",
                                "Caution: Check these options from manual of Annovar.")
    group.add_option("--table_annovar", action="store", help="The Annovar perl script of table_annovar.pl.",dest="table_annovar", metavar="./table_annovar.pl")
    group.add_option("--annovar_database_locat", dest="annovar_database_locat", action="store",
            help="The  database location/dir for the annotation datasets.", metavar="humandb")
    parser.add_option_group(group)
    (options, args) = parser.parse_args()
    
    if len(sys.argv)==1: # If the user does not provide any parameters, print help information and exit the script.
        parser.print_help()
        sys.exit()

    print("%s" %begin_description)
    print("%s" %version)
    print("Notice: Your command of SVOC is %s" % sys.argv[:])
    
    config_file = os.path.join(os.path.dirname(__file__),"config.ini") 
    if os.path.isfile(config_file):
        config.read(config_file)
        sections = config.sections()
        for section in sections:
            ConfigSectionMap(config,section)    
    else:
        print("Error: The default configure file of [ config.ini ] is not here, exit! Please redownload the SVOC.")
        sys.exit()

    # begin to process user's options:
    if options.config != None: # User defined configuration file.
        if os.path.isfile(options.config):
            config.read(options.config)
            sections = config.sections()
            for section in sections:
                ConfigSectionMap(config,section)
        else:
            print("Error: The config file [ %s ] is not here,please check the path of your config file." % options.config)
            sys.exit()
    if options.buildver != None:
        paras['buildver']=options.buildver
    if options.table_annovar != None:
        if os.path.isfile(options.table_annovar):
            paras['table_annovar']=options.table_annovar
        else:
            print("Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar" 
                    % options.table_annovar)
    if options.annovar_database_locat != None:
        paras['annovar_database_locat']=options.annovar_database_locat
    if options.input != None:
        paras['inputfile']=options.input
    if options.output != None:
        paras['outfile']=options.output
    if options.output_type != None:
        paras['output_type']=options.output_type
    if options.database_svoc != None:
        paras['database_svoc']=options.database_svoc
        paras['mane'] = paras['database_svoc']+'/MANE.GRCh38.v1.4.refseq_genomic.gtf.gz'
        paras['ref_fasta'] = paras['database_svoc']+'/ref_fasta'
        paras['gnomadv2'] = paras['database_svoc']+'/gnomadv2'
        
    if not os.path.isfile(paras['inputfile']):
        print("Error: Your input file [ %s ] is not here,please check the path of your input file." % paras['inputfile'])
        sys.exit()
    print ("INFO: The options are %s " % paras)

    buildver = paras['buildver']
    vcf_input_file = paras['inputfile']
    output_prefix = paras['outfile']
    output_type = paras['output_type']
    table_annovar = paras['table_annovar']
    annovar_db_path = paras['annovar_database_locat'] 
    transvar_file = f"{output_prefix}.{paras['buildver']}_transvar.txt"  
    annovar_file = f"{output_prefix}.{paras['buildver']}_multianno.txt"
    transvar_filter_file = f"{output_prefix}.{paras['buildver']}_transvar.filter.txt"
    mane_gtf_file = paras['mane']
    mark_file = f"{output_prefix}.{paras['buildver']}_transvar_mark.txt"
    if buildver == 'hg19':
        ref_fasta_file = f"{paras['ref_fasta']}/hg19.fa"
        gnomadv2_genomes = f"{paras['gnomadv2']}/gnomad.genomes.r2.1.1.sites.vcf.bgz"
        gnomadv2_exomes = f"{paras['gnomadv2']}/gnomad.exomes.r2.1.1.sites.vcf.bgz"
    elif buildver == 'hg38':
        ref_fasta_file = f"{paras['ref_fasta']}/hg38.fa"
        gnomadv2_genomes = f"{paras['gnomadv2']}/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
        gnomadv2_exomes = f"{paras['gnomadv2']}/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
    svoc_input_file = f"{output_prefix}.{paras['buildver']}_svocinput.txt"
    svoc_output_json_file = f"{output_prefix}.{paras['buildver']}_svocoutput.json"
    svoc_output_csv_file = f"{output_prefix}.{paras['buildver']}_svocoutput.csv"
    svoc_output_txt_file = f"{output_prefix}.{paras['buildver']}_svocoutput.txt"
    
    if not os.path.exists(svoc_input_file):
        print("Start generating SVOC input file...")
        # Annotated by ANNOVAR and TransVar
        if not os.path.exists(annovar_file):
#            annovar_command = (
#                f"perl {table_annovar} {vcf_input_file} {annovar_db_path} "
#                f"-vcfinput -build {buildver} -out {output_prefix} "
#                "-protocol refGene,avsnp151,dbnsfp47a_interpro,clinvar_20240917,intervar_20180118,dbnsfp47a,dbscsnv11 "
#                "-operation g,f,f,f,f,f,f -nastring . -thread 2 -remove "
#                "-arg '-hgvs,\",\",\",\",\",\"'"
#            )
            annovar_command = (
                f"perl {table_annovar} {vcf_input_file} {annovar_db_path} "
                f"-vcfinput -build {buildver} -out {output_prefix} "
                "-protocol refGene,dbnsfp47a,dbscsnv11 "
                "-operation g,f,f -nastring . -thread 1 -remove "
                "-arg '-hgvs,\",\"'"
            )
            # Execute ANNOVAR command
            try:
                print("Running ANNOVAR annotation...")
                subprocess.run(annovar_command, shell=True, check=True)
                print("ANNOVAR annotation completed.")
            except subprocess.CalledProcessError as e:
                print(f"Error running ANNOVAR command: {e}")
        else:
            print("ANNOVAR annotation file already exists.")

        if not os.path.exists(transvar_file):
            transvar_command = (
                f"transvar ganno --vcf {vcf_input_file} --refseq --refversion {buildver}| "
                f"awk '{{gsub(/  +/, \"\\t\", $0); print}}' > {transvar_file}"
            )
            # Execute TransVar command
            try:
                print("Running TransVar annotation...")
                subprocess.run(transvar_command, shell=True, check=True)
                print("TransVar annotation completed.")
            except subprocess.CalledProcessError as e:
                print(f"Error running TransVar command: {e}")
        else:
            print("TransVar annotation file already exists.")
        
        if not os.path.exists(transvar_filter_file):
            # Filter TransVar annotated transcript
            transvar_transcript_filter(transvar_file, transvar_filter_file, mane_gtf_file, mark_file)
            print("Transcript filtering completed.")
        else:
            print("Transcript filtering has already been completed.")
        
        merge_transvar_annovar(transvar_filter_file, annovar_file, svoc_input_file)
        print("SVOC input file generation completed.")
    else:
        print("SVOC input file already exists.")
# ----------------------------Part 1: Processing tsv input information.-----------------------------

    df_output = pd.DataFrame(
        columns=['Variant_ID', 
                'Gene',
                'Transcript',
                'Transcript_tag',
                'Variant_Information',
                'Point', 
                'Classification', 
                'Evidence_Codes',
                "SVOC_Result_Source",
                "Functional_Result_Source",
                'OncoKB_Classification',
                "Expert_Experience",
                'OVS1','OVS1_Basis',
                'OS1','OS1_Basis',
                'OS2','OS2_Basis',
                'OS3','OS3_Basis',
                'OM1','OM1_Basis',
                'OM2','OM2_Basis',
                'OM3','OM3_Basis',
                'OM4','OM4_Basis',
                'OP1','OP1_Basis',
                'OP2','OP2_Basis',
                'OP3','OP3_Basis',
                'OP4','OP4_Basis',
                'SBVS1','SBVS1_Basis',
                'SBS1','SBS1_Basis',
                'SBS2','SBS2_Basis',
                'SBP1','SBP1_Basis',
                'SBP2','SBP2_Basis',
                'vChr','vPos','vrsID','vRef','vAlt'
                ])
    df = pd.read_csv(svoc_input_file, sep='\t')
    for index, row in df.iterrows():
        vChr = row['vChr']
        vPos = row['vPos']
        vrsID = row['vrsID']
        vRef = row['vRef']
        vAlt = row['vAlt']
        Variant_ID = index+1
        Variant_Information = row['coordinate']
        gene = row['Gene.refGene'] # annovar gene
        Transcript = re.search(r'(\w+)_\d+', row['transcript']).group(0) if isinstance(row['transcript'], str) and re.search(r'(\w+)_\d+', row['transcript'])  else ''
        Transcript_tag = row['transcript_tag']
        var_type = row['info'].split('CSQN=')[1].split(';')[0]
        AAchange_single = re.search(r'p\.([^/]+)', row['coordinate']).group() if isinstance(row['coordinate'], str) and re.search(r'p\.([^/]+)', row['coordinate']) else ''
        if AAchange_single == "p.(=)": # Synonymous
            AAchange_single = re.search(r'p\.(.+)', row['AAChange.refGene'].split(',')[0]).group() if isinstance(row['AAChange.refGene'], str) and re.search(r'p\.(.+)', row['AAChange.refGene'].split(',')[0]) else ''
        # row['coordinate'] = chr7:g.140477838_140477852del15/c.1457_1471del15/p.N486_P490delNVTAP
        if 'del' in AAchange_single and 'delins' not in AAchange_single:
            AAchange_single = re.match(r'[^del]*del', AAchange_single).group()
        cDNA = re.search(r'c\.([^/]+)', row['coordinate']).group(0) if isinstance(row['coordinate'], str) and re.search(r'c\.([^/]+)', row['coordinate']) else ''
        gDNA = re.search(r'g\.([^/]+)', row['coordinate']).group(0) if isinstance(row['coordinate'], str) and re.search(r'g\.([^/]+)', row['coordinate']) else ''
        refAA = ''  
        posAA = ''  
        altAA = ''  
        if AAchange_single:
            aa = re.search(r'p\.(.*)', AAchange_single).group(1) # without 'p.', eg：V600E
            refAA, posAA, altAA, startAA, endAA, AAchange_position = getAAchangePart(AAchange_single)
        else:
            aa = ''
            refAA = ''
            posAA = ''
            altAA = ''
            startAA = ''
            endAA = ''
            AAchange_position = ''
        chr = row['Chr'] # int
        pos = '' # str
        ref = row['Ref']
        alt = row['Alt']

        # Obtain the position of base changes, including two formats:：g.140453136_140453137delinsTT or g.140453134T>C
        match_gDNA = re.match(r"g\.(\d+(?:_\d+)?)(?:([A-Z])>([A-Z]))?", gDNA)
        if match_gDNA:
            pos = match_gDNA.group(1)

        # Initialize classification information.
        point = 0 # Point
        classification = '' # Classification
        standard = '' # Evidence_Codes
        oncokb_class = '/' # OncoKB_Classification
        oncokb_Description = '/' # OncoKB_Description
        expert_experience = '/' # Expert_Experience
        Functional_Result_Source = '/' # Expert/ClinGen/OncoKB/none
        expertfunc_class = '/'
        expertfunc_Description = '/'

        OVS1 = 0
        OVS1_Basis = '/'
        OS1 = 0
        OS1_Basis = '/'
        OS2 = 0
        OS2_Basis = '/'
        OS3 = 0
        OS3_Basis = '/'
        OM1 = 0
        OM1_Basis = '/'
        OM2 = 0
        OM2_Basis = '/'
        OM3 = 0
        OM3_Basis = '/'
        OM4 = 0
        OM4_Basis = '/'
        OP1 = 0
        OP1_Basis = '/'
        OP2 = 0
        OP2_Basis = '/'
        OP3 = 0
        OP3_Basis = '/'
        OP4 = 0
        OP4_Basis = '/'
        SBVS1 = 0
        SBVS1_Basis = '/'
        SBS1 = 0
        SBS1_Basis = '/'
        SBS2 = 0
        SBS2_Basis = '/'
        SBP1 = 0
        SBP1_Basis = '/'
        SBP2 = 0
        SBP2_Basis = '/'


# ----------------------Part 2: Processing Information from the Database-----------------------------------
        isCancerHotspots = False
        isCOSMICHotspots = False
        isExpertCuratedHotspots = False
        CancerHotspots_sample = 0
        CancerHotspots_count = 0
        COSMIC_sample = 0
        COSMIC_count = 0

        if aa != '':
            isCancerHotspots, CancerHotspots_sample, CancerHotspots_count= db_handler.getHotspots(gene, refAA, altAA, posAA, aa)
        if AAchange_single != '':
            isCOSMICHotspots, COSMIC_sample, COSMIC_count = db_handler.getCOSMICHotspots(gene, AAchange_single, AAchange_position)
        isExpertCuratedHotspots, ExpertEvidenceCode, ExpertCuratedHotspots_Basis = db_handler.getExpertCuratedHotspots(gene, gDNA, buildver)
        # gnomad_result = query_gnomAD_local(chr, pos, ref, alt, gnomadv2_exomes) if query_gnomAD_local(chr, pos, ref, alt, gnomadv2_exomes) else query_gnomAD_local(chr, pos, ref, alt, gnomadv2_genomes)
        
        inGnomAD = False
        MAF = -1
        Max_AC = 0
        Max_AN = 0
        
        gnomad_result = query_gnomAD_local(chr, pos, ref, alt, gnomadv2_exomes)
        if gnomad_result:
            Continent, Max_AC, Max_AN, MAF = gnomad_result
            if MAF == '.' and Max_AC == '.' and Max_AN == '.' and Continent == '.':
                gnomad_result = query_gnomAD_local(chr, pos, ref, alt, gnomadv2_genomes)
        else:
            gnomad_result = query_gnomAD_local(chr, pos, ref, alt, gnomadv2_genomes)

        if gnomad_result:# not none
            Continent, Max_AC, Max_AN, MAF = gnomad_result
            if MAF != '.' and Max_AC != '.' and Max_AN != '.' and Continent != '.':
                MAF = float(MAF)
                Max_AC = int(Max_AC)
                Max_AN = int(Max_AN)
                inGnomAD = True
            else: # not none and not all '.'
                inGnomAD = False
                MAF = -1
                Max_AC = 0
                Max_AN = 0
        else: # none
            inGnomAD = False
            MAF = -1
            Max_AC = 0
            Max_AN = 0
           

        # functional data information
        vcep_PS3 = False
        vcep_BS3 = False
        oncokb_O = False
        oncokb_N = False
        expertfunc_OS2 = False
        expertfunc_SBS2 = False
        vcep_metcodes = db_handler.getVCEPMetCodes(Transcript, cDNA)
        oncokb_class, oncokb_Description = db_handler.getOncoKBRes(aa,gene)
        if oncokb_class is None:
            oncokb_class = '/'
        expertfunc_class, expertfunc_Description = db_handler.getExpertFuncRes(aa,gene,gDNA)

        if vcep_metcodes:
            vcep_PS3 = "PS3" in vcep_metcodes
            vcep_BS3 = "BS3" in vcep_metcodes
        oncokb_O = (oncokb_class == "Oncogenic" or oncokb_class == "Likely Oncogenic") 
        oncokb_N = oncokb_class == "Likely Neutral"
        expertfunc_OS2 = (expertfunc_class == "OS2")
        expertfunc_SBS2 = (expertfunc_class == "SBS2")

        # Tumor suppressor gene/oncogene
        isONG = False
        isTSG = False
        isONG = db_handler.isOncoGene(gene)
        isTSG = db_handler.isTSG(gene)

        # functional domain information
        domain_name = db_handler.getOncokbDomain(gene,startAA,endAA)

        # Predictive evidence
        isOVS1, OVS1_strength, consequence = check_evidence.isOVS1(chr, pos, ref, alt, gDNA, isTSG, buildver, vPos, vRef, vAlt)
        
        SIFT_pred = row['SIFT_pred']
        MutationAssessor_pred = row['MutationAssessor_pred']
        FATHMM_pred = row['FATHMM_pred']
        Polyphen2_HDIV_pred = row['Polyphen2_HDIV_pred']
        Polyphen2_HVAR_pred = row['Polyphen2_HVAR_pred']
        MutationTaster_pred = row['MutationTaster_pred']
        CADD_phred = row['CADD_phred']
        REVEL_pred = row['REVEL_score']
        # Splicing effect prediction
        DS_AG, DS_AL, DS_DG, DS_DL = spliceai.getSpliceAI('chr'+str(chr), int(pos), str(ref), str(alt), buildver, ref_fasta_file)
        if max(DS_AG, DS_AL, DS_DG, DS_DL) >= 0.5:
            SpliceAI_pred = True
        else:
            SpliceAI_pred = False
        dbscSNV_ada_score = row['dbscSNV_ADA_SCORE']
        dbscSNV_rf_score = row['dbscSNV_RF_SCORE']

        # Single Genetic Etiology
        isSingleGeneticEtiology, SingleGeneticEtiologyDescription = db_handler.getSingelGeneticEtiology(gene)

        # Same amino acid change as a previously established oncogenic variant
        HasSameAAchange = False
        SameAAchangeDescription = ''
        if AAchange_single != '':
            HasSameAAchange,SameAAchangeDescription = db_handler.getSameAAchange(ref, alt, AAchange_single, gene, buildver)

        # Missense variant at an amino acid residue where a different missense variant determined to be oncogenic
        isOM4 = False
        if ('del' not in AAchange_single) and ('delins' not in AAchange_single) and ('ins' not in AAchange_single) and ('dup' not in AAchange_single) and ('_' not in AAchange_single):
            SameAAresidue = None
            GranthmsDistance = 0
            GranthmsDistance_SameAAresidue = 0
            if AAchange_single != '' and refAA not in ['*', '?', '='] and altAA not in ['*', '?', '=']:
                SameAAresidue = db_handler.getSameAAresidue(AAchange_single, gene, AAchange_position)
            # Amino acid difference from reference amino acid should be greater or at least approximately the same as for missense change determined to be oncogenic.
            if SameAAresidue:
                #print(SameAAresidue)
                GranthmsDistance = db_handler.getGranthmsDistance(refAA, altAA)
                if isinstance(SameAAresidue, str):
                    samerefAA, sameposAA, samealtAA, samestartAA, sameendAA, sameAAchange_position= getAAchangePart(SameAAresidue)
                    if samerefAA not in ['*', '?', '='] and samealtAA not in ['*', '?', '=']:
                        GranthmsDistance_SameAAresidue = db_handler.getGranthmsDistance(samerefAA, samealtAA)
                        if GranthmsDistance >= GranthmsDistance_SameAAresidue:
                            isOM4 = True
                        else:
                            isOM4 = False
                elif isinstance(SameAAresidue, list) and all(isinstance(item, str) for item in SameAAresidue):
                    for item in SameAAresidue:
                        samerefAA, sameposAA, samealtAA, samestartAA, sameendAA, sameAAchange_position = getAAchangePart(item)
                        if samerefAA not in ['*', '?', '='] and samealtAA not in ['*', '?', '=']:
                            GranthmsDistance_SameAAresidue = db_handler.getGranthmsDistance(samerefAA, samealtAA)
                            if GranthmsDistance >= GranthmsDistance_SameAAresidue:
                                isOM4 = True
                                SameAAresidue = item
                                break
                            else:
                                isOM4 = False
            else:
                isOM4 = False


# ----------------------Part Three: Scoring-----------------------------------
        # Evidence of computational tools
        if check_evidence.isSBP1(SIFT_pred,
                           MutationAssessor_pred,
                           FATHMM_pred,
                           Polyphen2_HDIV_pred,
                           Polyphen2_HVAR_pred,
                           MutationTaster_pred,
                           SpliceAI_pred,
                           CADD_phred,
                           REVEL_pred,
                           dbscSNV_ada_score, 
                           dbscSNV_rf_score):
            standard = standard+'SBP1;'
            point -= 1
            SBP1 = 1
            SBP1_Basis = 'SIFT_pred:'+str(SIFT_pred)+',MutationAssessor_pred:'+str(MutationAssessor_pred)+',FATHMM_pred:'+str(FATHMM_pred)+',Polyphen2_HDIV_pred:'+str(Polyphen2_HDIV_pred)+',Polyphen2_HVAR_pred:'+str(Polyphen2_HVAR_pred)+',MutationTaster_pred:'+str(MutationTaster_pred)+',SpliceAI_pred:'+str(SpliceAI_pred)+',CADD_phred:'+str(CADD_phred)+',REVEL_pred:'+str(REVEL_pred)+'dbscSNV_ada_score:'+str(dbscSNV_ada_score)+',dbscSNV_rf_score:'+str(dbscSNV_rf_score)
        elif check_evidence.isOP1(SIFT_pred,
                            MutationAssessor_pred,
                            FATHMM_pred,
                            Polyphen2_HDIV_pred,
                            Polyphen2_HVAR_pred,
                            MutationTaster_pred,
                            SpliceAI_pred,
                            CADD_phred,
                            REVEL_pred,
                            dbscSNV_ada_score, 
                            dbscSNV_rf_score):
            standard = standard+'OP1;'
            point += 1
            OP1 = 1
            OP1_Basis = 'SIFT_pred:'+str(SIFT_pred)+',MutationAssessor_pred:'+str(MutationAssessor_pred)+',FATHMM_pred:'+str(FATHMM_pred)+',Polyphen2_HDIV_pred:'+str(Polyphen2_HDIV_pred)+',Polyphen2_HVAR_pred:'+str(Polyphen2_HVAR_pred)+',MutationTaster_pred:'+str(MutationTaster_pred)+',SpliceAI_pred:'+str(SpliceAI_pred)+',CADD_phred:'+str(CADD_phred)+',REVEL_pred:'+str(REVEL_pred)+'dbscSNV_ada_score:'+str(dbscSNV_ada_score)+',dbscSNV_rf_score:'+str(dbscSNV_rf_score)
        
        # Population evidence
        if check_evidence.isSBVS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'SBVS1;'
            point -= 8
            SBVS1 = 1
            SBVS1_Basis = "Allele frequency of " + str(MAF) + "("+ str(Max_AC) + "/" + str(Max_AN) +")" + " in the "+ Continent + " subpopulation of the gnomAD(v2.1.1 controls)."
        elif check_evidence.isSBS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'SBS1;'
            point -= 4
            SBS1 = 1
            SBS1_Basis = "Allele frequency of " + str(MAF) + "("+ str(Max_AC) + "/" + str(Max_AN) +")" + " in the "+ Continent + " subpopulation of the gnomAD(v2.1.1 controls)."
        elif check_evidence.isOP4(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'OP4;'
            if MAF < 0:
                OP4_Basis = "Absent from controls (gnomAD v2.1.1)."
            else:
                OP4_Basis =str(MAF) + ", extremely low frequency in gnomAD."
            point += 1
            OP4 = 1
            
        
        # Others
        if check_evidence.isOM1(domain_name):
            if not check_evidence.isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count) and not check_evidence.isOS1(HasSameAAchange):
                standard = standard+'OM1;'
                point += 2
                OM1 = 1
                # May be located in multiple critical functional domains simultaneously, for example, MSH2's ['Muts_III ','Muts_IV'].
                if isinstance(domain_name, str):
                    OM1_Basis = "Located in well-established functional domain: " + domain_name
                elif isinstance(domain_name, list):
                    OM1_Basis = "Located in well-established functional domain: " + "; ".join(domain_name)


        if check_evidence.isOM2(isONG, isTSG, var_type):
            if not isOVS1:
                standard = standard+'OM2;'
                point += 2
                OM2 = 1
                if (isONG or isTSG) and ("inframe_deletion" in var_type or "inframe_insertion" in var_type):
                    OM2_Basis = "protein length due to in-frame deletions/insertions in known oncogene or tumor suppressor gene"
                elif isTSG and ("stop_lost" in var_type):
                    OM2_Basis = "protein length due to stop-loss variants in a known tumor suppressor gene"
        
        # Cancer hotspots
        if check_evidence.isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            if not check_evidence.isOS1(HasSameAAchange):
                standard = standard+'OS3;'
                OS3 = 1
                if isCancerHotspots:
                    OS3_Basis = "hotspot in cancerhotspots.org, AA position count is >= 50 (" + str(CancerHotspots_sample) + ") and AA change count is ≥ 10 (" + str(CancerHotspots_count) +")"
                elif not isCancerHotspots and isCOSMICHotspots:
                    OS3_Basis = "hotspot in COSMIC, AA position count is >= 50 (" + str(COSMIC_sample) + ") and AA change count is ≥ 10 (" + str(COSMIC_count) +")"
                elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                    OS3_Basis = ExpertCuratedHotspots_Basis
                point += 4


        elif check_evidence.isOM3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            if not check_evidence.isOM1(domain_name) and not isOM4:
                standard = standard+'OM3;'
                point += 2
                OM3 = 1
                if isCancerHotspots:
                    OM3_Basis = "hotspot in cancerhotspots.org, AA position count is < 50 (" + str(CancerHotspots_sample) + ") and AA change count is ≥ 10 (" + str(CancerHotspots_count) +")"
                elif not isCancerHotspots and isCOSMICHotspots:
                    OM3_Basis = "hotspot in COSMIC, AA position count is < 50 (" + str(COSMIC_sample) + ") and AA change count is ≥ 10 (" + str(COSMIC_count) +")"
                elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                    OM3_Basis = ExpertCuratedHotspots_Basis

        elif check_evidence.isOP3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            standard = standard+'OP3;'  
            point += 1
            OP3 = 1
            if isCancerHotspots:
                OP3_Basis = "hotspot in cancerhotspots.org, AA change count is < 10 (" + str(CancerHotspots_count) +")"
            elif not isCancerHotspots and isCOSMICHotspots:
                OP3_Basis = "hotspot in COSMIC, AA change count is < 10 (" + str(COSMIC_count) +")"
            elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                OP3_Basis = ExpertCuratedHotspots_Basis


        # functional evidence
        if check_evidence.isOS2(vcep_PS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
            if not check_evidence.isOS1(HasSameAAchange):
                standard = standard+'OS2;'
                point += 4
                OS2 = 1
                if expertfunc_OS2:
                    OS2_Basis = expertfunc_Description.replace("\n", " ").replace("\r", " ")
                    OS2_Basis = " ".join(OS2_Basis.split())
                    Functional_Result_Source = 'Expert'
                elif oncokb_O:
                    OS2_Basis = oncokb_Description.replace("\n", " ").replace("\r", " ")
                    OS2_Basis = " ".join(OS2_Basis.split())
                    Functional_Result_Source = 'OncoKB'
                else:
                    OS2_Basis = "PS3 in ClinGen VCEP."
                    Functional_Result_Source = 'ClinGen'
        elif check_evidence.isSBS2(vcep_BS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
            standard = standard+'SBS2;'
            point -= 4
            SBS2 = 1
            if expertfunc_SBS2:
                SBS2_Basis = expertfunc_Description.replace("\n", " ").replace("\r", " ")
                SBS2_Basis = " ".join(SBS2_Basis.split())
                Functional_Result_Source = 'Expert'
            elif oncokb_N:
                SBS2_Basis = oncokb_Description.replace("\n", " ").replace("\r", " ")
                SBS2_Basis = " ".join(SBS2_Basis.split())
                Functional_Result_Source = 'OncoKB'
            else:
                SBS2_Basis = "BS3 in ClinGen VCEP."
                Functional_Result_Source = 'ClinGen'

    
        # Predictive evidence
        if check_evidence.isSBP2(var_type, SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
            standard = standard+'SBP2;'
            point -= 1
            SBP2 = 1
            SBP2_Basis = "A synonymous variant for which splicing prediction algorithms predict no impact on splicing."

        if isOVS1:
            standard = standard+'OVS1;'
            point += 8
            OVS1 = 1
            OVS1_Basis = consequence + " in tumor suppressor gene"           
        
        if check_evidence.isOP2(isSingleGeneticEtiology):
            standard = standard+'OP2;'
            point += 1
            OP2 = 1
            OP2_Basis = SingleGeneticEtiologyDescription.replace("\n", " ").replace("\r", " ")
            OP2_Basis = " ".join(OP2_Basis.split())

        if check_evidence.isOS1(HasSameAAchange):
            standard = standard+'OS1;'
            point += 4
            OS1 = 1
            OS1_Basis = SameAAchangeDescription.replace("\n", " ").replace("\r", " ")
            OS1_Basis = " ".join(OS1_Basis.split())

        if isOM4:
            if not check_evidence.isOS1(HasSameAAchange) and not check_evidence.isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count) and not check_evidence.isOM1(domain_name):
                standard = standard+'OM4;'
                point += 2
                OM4 = 1
                OM4_Basis = f'Missense variant at an amino acid residue where a different missense variant determined to be oncogenic (using this standard) has been documented. Amino acid difference({AAchange_single}:{GranthmsDistance}) from reference amino acid is >= missense change determined to be oncogenic({SameAAresidue}:{GranthmsDistance_SameAAresidue}).'

        
# ----------------------Part Four: Classification-----------------------------------

        if point <= -7:
            classification = "Benign"
        elif point >= -6 and point <= -1:
            classification = "Likely Benign"
        elif point >= 0 and point <= 5:
            classification = "VUS"
        elif point >=6 and point <= 9:
            classification = "Likely Oncogenic"
        elif point >= 10:
            classification = "Oncogenic"
        
        
        result = {'Variant_ID': Variant_ID, 
                  'Gene':gene,
                  'Transcript': Transcript,
                  'Transcript_tag': Transcript_tag,
                  'Variant_Information': Variant_Information,
                  'Point': point,
                  'Classification': classification, 
                  'Evidence_Codes': standard,
                  'SVOC_Result_Source': 'Automated', 
                  'Functional_Result_Source': Functional_Result_Source,
                  'OncoKB_Classification': oncokb_class,
                  'Expert_Experience': expert_experience, 
                  'OVS1': OVS1,
                  'OVS1_Basis': OVS1_Basis,
                  'OS1': OS1,
                  'OS1_Basis': OS1_Basis,
                  'OS2': OS2,
                  'OS2_Basis': OS2_Basis,
                  'OS3': OS3,
                  'OS3_Basis': OS3_Basis,
                  'OM1': OM1,
                  'OM1_Basis': OM1_Basis,
                  'OM2': OM2,
                  'OM2_Basis': OM2_Basis,
                  'OM3': OM3,
                  'OM3_Basis': OM3_Basis,
                  'OM4': OM4,
                  'OM4_Basis': OM4_Basis,
                  'OP1': OP1,
                  'OP1_Basis': OP1_Basis,
                  'OP2': OP2,
                  'OP2_Basis': OP2_Basis,
                  'OP3': OP3,
                  'OP3_Basis': OP3_Basis,
                  'OP4': OP4,
                  'OP4_Basis': OP4_Basis,
                  'SBVS1': SBVS1,
                  'SBVS1_Basis': SBVS1_Basis,
                  'SBS1': SBS1,
                  'SBS1_Basis': SBS1_Basis,
                  'SBS2': SBS2,
                  'SBS2_Basis': SBS2_Basis,
                  'SBP1': SBP1,
                  'SBP1_Basis': SBP1_Basis,
                  'SBP2': SBP2,
                  'SBP2_Basis': SBP2_Basis,
                  'vChr': vChr,
                  'vPos': vPos,
                  'vrsID': vrsID,
                  'vRef': vRef,
                  'vAlt': vAlt}

        df_output = pd.DataFrame([result])
        if index == 0:
            # If it is the first variant, write it directly to the file.
            if output_type == 'json':
                df_output.to_json(svoc_output_json_file, orient='records', lines=True, mode='w')
                df_output.to_csv(svoc_output_txt_file, sep='\t', index=False, mode='w')
            elif output_type == 'csv':
                df_output.to_csv(svoc_output_csv_file, index=False, mode='w')
                df_output.to_csv(svoc_output_txt_file, sep='\t', index=False, mode='w')
            else:
                df_output.to_csv(svoc_output_txt_file, sep='\t', index=False, mode='w')
        else:
            # If it's not the first variant, append the written file.
            if output_type == 'json':
                df_output.to_json(svoc_output_json_file, orient='records', lines=True, mode='a')
                df_output.to_csv(svoc_output_txt_file, sep='\t', mode='a', index=False, header=False)
            elif output_type == 'csv':
                df_output.to_csv(svoc_output_csv_file, index=False, mode='a', header=False)
                df_output.to_csv(svoc_output_txt_file, sep='\t', mode='a', index=False, header=False)
            else:
                df_output.to_csv(svoc_output_txt_file, sep='\t', mode='a', index=False, header=False)

    # svocoutput.filter
    svoc_output_txt_filter_file = f"{output_prefix}.{paras['buildver']}_svocoutput.filter.txt"
    svoc_output_filter(svoc_output_txt_file, svoc_output_txt_filter_file)   

    print("%s" %end_description)

if __name__ == "__main__":
    main()
