#!/root/anaconda3/envs/svocenv/bin/python
import os, asyncio, re, getopt, sys, pandas as pd, textwrap, subprocess, copy, logging, io, time, platform, optparse, gzip, glob
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1' 
import modules.db_handler as db_handler
import modules.check_evidence as function
import modules.spliceai_prediction as spliceai
from modules.query_gnomad import query_gnomAD_api, query_gnomAD_local
from modules.generate_input import transvar_transcript_filter, merge_transvar_annovar
if platform.python_version()< '3.0.0' :
    import ConfigParser
else:
    import configparser

prog="SVOC"
version = """Version: 20250218
Written by Leqing SANG, yqsang23@m.fudan.edu.cn. 
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

def main():
    if platform.python_version()< '3.0.0'  :
        config=ConfigParser.ConfigParser()
    else:
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
        paras['mane'] = paras['database_intervar']+'/MANE.GRCh38.v1.4.ensembl_genomic.gtf'
        paras['ref_fasta'] = paras['database_intervar']+'/ref_fasta'
        
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
    svoc_input_file = f"{output_prefix}.{paras['buildver']}_svocinput.txt"
    svoc_output_json_file = f"{output_prefix}.{paras['buildver']}_svocoutput.json"
    svoc_output_csv_file = f"{output_prefix}.{paras['buildver']}_svocoutput.csv"
    svoc_output_txt_file = f"{output_prefix}.{paras['buildver']}_svocoutput.txt"

    # Annotated by ANNOVAR and TransVar
    # Execute ANNOVAR command
    annovar_command = (
        f"perl {table_annovar} {vcf_input_file} {annovar_db_path} "
        f"-vcfinput -build {buildver} -out {output_prefix} "
        "-protocol refGene,avsnp151,dbnsfp47a_interpro,clinvar_20240917,intervar_20180118,dbnsfp47a,dbscsnv11 "
        "-operation g,f,f,f,f,f,f -nastring . -thread 2 -remove "
        "-arg '-hgvs,\",\",\",\",\",\"'"
    )
    #print(annovar_command)
    try:
        print("Running ANNOVAR annotation...")
        subprocess.run(annovar_command, shell=True, check=True)
        print("ANNOVAR annotation completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error running ANNOVAR command: {e}")
    
    # Execute TransVar command
    transvar_command = (
        f"transvar ganno --vcf {vcf_input_file} --refseq | "
        f"awk '{{gsub(/  +/, \"\\t\", $0); print}}' > {transvar_file}"
    )
    #print(transvar_command)
    try:
        print("Running TransVar annotation...")
        subprocess.run(transvar_command, shell=True, check=True)
        print("TransVar annotation completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error running TransVar command: {e}")
    # # Filter TransVar annotated transcript
    transvar_transcript_filter(transvar_file, transvar_filter_file, mane_gtf_file)
    print("TransVar transcripts filtering completed.")

    # Merge ANNOVAR and TransVar into SVOC input file
    merge_transvar_annovar(transvar_filter_file, annovar_file, svoc_input_file)
    print("SVOC input file has been generated.")

# ----------------------------Part 1: Processing tsv input information.-----------------------------

    df_output = pd.DataFrame(
        columns=['Variant_ID', 
                'Gene',
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
                'SBP2','SBP2_Basis'])
    row_index = 0  # Initialize row index.
    df = pd.read_csv(svoc_input_file, sep='\t')
    for index, row in df.iterrows():
        Variant_ID = index+1
        Variant_Information = row['coordinate']
        gene = row['gene']
        Transcript = re.search(r'(\w+)_\d+', row['transcript']).group(0) if isinstance(row['transcript'], str) and re.search(r'(\w+)_\d+', row['transcript'])  else ''
        var_type = row['info'].split('CSQN=')[1].split(';')[0]
        AAchange_single = re.search(r'p\.([^/]+)', row['coordinate']).group() if isinstance(row['coordinate'], str) and re.search(r'p\.([^/]+)', row['coordinate']) else ''
        cDNA = re.search(r'c\.([^/]+)', row['coordinate']).group(0) if isinstance(row['coordinate'], str) and re.search(r'c\.([^/]+)', row['coordinate']) else ''
        gDNA = re.search(r'g\.([^/]+)', row['coordinate']).group(0) if isinstance(row['coordinate'], str) and re.search(r'g\.([^/]+)', row['coordinate']) else ''
        refAA = ''  
        posAA = ''  
        altAA = ''  
        if AAchange_single:
            aa = re.search(r'p\.(.*)', AAchange_single).group(1) # without 'p.', eg：V600E
        else:
            aa = ''
        chr = row['Chr'] # int
        pos = '' # str
        ref = row['Ref']
        alt = row['Alt']
        
        if ('del' in AAchange_single or 'delins' in AAchange_single) and '_' in AAchange_single:
            startAA = re.match(r"p\.[A-Z](\d+)_[A-Z](\d+)", AAchange_single).group(1)
            endAA = re.search(r"p\.[A-Z](\d+)_[A-Z](\d+)", AAchange_single).group(2)
            posAA = startAA + '_' + endAA
            refAA = re.match(r"p\.([A-Z])\d+_([A-Z])\d+", AAchange_single).group(1) + re.match(r"p\.([A-Z])\d+_([A-Z])\d+", AAchange_single).group(2)
            altAA = ''
        elif AAchange_single:
            refAA = re.search(r"p\.([A-Za-z]|\*)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(1)
            posAA = re.search(r"p\.([A-Za-z]|\*|)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(2)
            altAA = re.search(r"p\.([A-Za-z]|\*|)(\d+)([A-Za-z]|\*|\?|\=)", AAchange_single).group(3)
            startAA = posAA
            endAA = posAA
        else:
            refAA = ''
            posAA = ''
            altAA = ''
            startAA = ''
            endAA = ''
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
        Functional_Result_Source = '/' # Functional_Result_Source：Expert/ClinGen/OncoKB/none
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
            isCOSMICHotspots, COSMIC_sample, COSMIC_count = db_handler.getCOSMICHotspots(gene, AAchange_single)
        isExpertCuratedHotspots, ExpertEvidenceCode, ExpertCuratedHotspots_Basis = db_handler.getExpertCuratedHotspots(gene, gDNA)

        result = query_gnomAD_local(chr, pos, ref, alt, "/mnt/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz") if query_gnomAD_local(chr, pos, ref, alt, "/mnt/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz") else query_gnomAD_local(chr, pos, ref, alt, "/mnt/gnomad/gnomad.genomes.r2.1.1.sites.vcf.bgz")
        if result:
            Continent, Max_AC, Max_AN, MAF = result
            inGnomAD = True
        else:
            inGnomAD = False
            MAF = -1
            Max_AC = 0
            Max_AN = 0

        # Functional data information
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

        # Functional domain information
        domain_name = db_handler.getOncokbDomain(gene,startAA,endAA)

        # Predictive evidence
        isOVS1, OVS1_strength, consequence = function.isOVS1(chr, pos, ref, alt, gDNA, isTSG)
        
        SIFT_pred = row['SIFT_pred']
        MutationAssessor_pred = row['MutationAssessor_pred']
        FATHMM_pred = row['FATHMM_pred']
        Polyphen2_HDIV_pred = row['Polyphen2_HDIV_pred']
        Polyphen2_HVAR_pred = row['Polyphen2_HVAR_pred']
        MutationTaster_pred = row['MutationTaster_pred']
        CADD_phred = row['CADD_phred']
        REVEL_pred = row['REVEL_score']
        # Splicing effect prediction
        DS_AG, DS_AL, DS_DG, DS_DL = spliceai.getSpliceAI('chr'+str(chr), int(pos), str(ref), str(alt))
        if max(DS_AG, DS_AL, DS_DG, DS_DL) >= 0.5:
            SpliceAI_pred = True
        else:
            SpliceAI_pred = False
        dbscSNV_ada_score = row['dbscSNV_ADA_SCORE']
        dbscSNV_rf_score = row['dbscSNV_RF_SCORE']


# ----------------------Part Three: Scoring-----------------------------------
        # Evidence of computational tools
        if function.isSBP1(SIFT_pred,
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
        elif function.isOP1(SIFT_pred,
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
        if function.isSBVS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'SBVS1;'
            point -= 8
            SBVS1 = 1
            SBVS1_Basis = str(MAF) + ", "+ str(Max_AC) + "/" + str(Max_AN) +" alleles in the "+ Continent + " subpopulation of the gnomAD"
        elif function.isSBS1(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'SBS1;'
            point -= 4
            SBS1 = 1
            SBS1_Basis = str(MAF) + ", "+ str(Max_AC) + "/" + str(Max_AN) +" alleles in the "+ Continent + " subpopulation of the gnomAD"
        elif function.isOP4(inGnomAD, MAF, gene, Max_AC, Max_AN):
            standard = standard+'OP4;'
            if MAF < 0:
                OP4_Basis = "Absent from controls (gnomAD v2.1.1)"
            else:
                OP4_Basis =str(MAF) + ", extremely low frequency in gnomAD"
            point += 1
            OP4 = 1
            
        
        # Others
        if function.isOM1(domain_name):
            if not function.isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
                standard = standard+'OM1;'
                point += 2
                # May be located in multiple critical functional domains simultaneously, for example, MSH2's ['Muts_III ','Muts_IV'].
                if isinstance(domain_name, str):
                    OM1_Basis = "Located in well-established functional domain: " + domain_name
                elif isinstance(domain_name, list):
                    OM1_Basis = "Located in well-established functional domain: " + "; ".join(domain_name)


        if function.isOM2(isONG, isTSG, var_type):
            if not isOVS1:
                standard = standard+'OM2;'
                point += 2
                OM2 = 1
                if (isONG or isTSG) and ("inframe_deletion" in var_type or "inframe_insertion" in var_type):
                    OM2_Basis = "protein length due to in-frame deletions/insertions in known oncogene or tumor suppressor gene"
                elif isTSG and ("stop_lost" in var_type):
                    OM2_Basis = "protein length due to stop-loss variants in a known tumor suppressor gene"
        
        # Cancer hotspots
        if function.isOS3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            standard = standard+'OS3;'
            OS3 = 1
            if isCancerHotspots:
                OS3_Basis = "hotspot in cancerhotspots.org, AA position count is >= 50 (" + str(CancerHotspots_sample) + ") and AA change count is ≥ 10 (" + str(CancerHotspots_count) +")"
            elif not isCancerHotspots and isCOSMICHotspots:
                OS3_Basis = "hotspot in COSMIC, AA position count is >= 50 (" + str(COSMIC_sample) + ") and AA change count is ≥ 10 (" + str(COSMIC_count) +")"
            elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                OS3_Basis = ExpertCuratedHotspots_Basis
            point += 4


        elif function.isOM3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            if not function.isOM1(domain_name):
                standard = standard+'OM3;'
                point += 2
                OM3 = 1
                if isCancerHotspots:
                    OM3_Basis = "hotspot in cancerhotspots.org, AA position count is < 50 (" + str(CancerHotspots_sample) + ") and AA change count is ≥ 10 (" + str(CancerHotspots_count) +")"
                elif not isCancerHotspots and isCOSMICHotspots:
                    OM3_Basis = "hotspot in COSMIC, AA position count is < 50 (" + str(COSMIC_sample) + ") and AA change count is ≥ 10 (" + str(COSMIC_count) +")"
                elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                    OM3_Basis = ExpertCuratedHotspots_Basis

        elif function.isOP3(isCancerHotspots, isCOSMICHotspots, isExpertCuratedHotspots, ExpertEvidenceCode, CancerHotspots_sample, CancerHotspots_count, COSMIC_sample, COSMIC_count):
            standard = standard+'OP3;'  
            point += 1
            OP3 = 1
            if isCancerHotspots:
                OP3_Basis = "hotspot in cancerhotspots.org, AA change count is < 10 (" + str(CancerHotspots_count) +")"
            elif not isCancerHotspots and isCOSMICHotspots:
                OP3_Basis = "hotspot in COSMIC, AA change count is < 10 (" + str(COSMIC_count) +")"
            elif not isCancerHotspots and not isCOSMICHotspots and isExpertCuratedHotspots:
                OP3_Basis = ExpertCuratedHotspots_Basis


        # Functional evidence
        if function.isOS2(vcep_PS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
            standard = standard+'OS2;'
            point += 4
            OS2 = 1
            if expertfunc_OS2:
                OS2_Basis = expertfunc_Description
                Functional_Result_Source = 'Expert'
            elif oncokb_O:
                OS2_Basis = oncokb_Description
                Functional_Result_Source = 'OncoKB'
            else:
                OS2_Basis = "PS3 in ClinGen VCEP."
                Functional_Result_Source = 'ClinGen'
        elif function.isSBS2(vcep_BS3, oncokb_O, oncokb_N, expertfunc_OS2, expertfunc_SBS2):
            standard = standard+'SBS2;'
            point -= 4
            SBS2 = 1
            if expertfunc_SBS2:
                SBS2_Basis = expertfunc_Description
                Functional_Result_Source = 'Expert'
            elif oncokb_N:
                SBS2_Basis = oncokb_Description
                Functional_Result_Source = 'OncoKB'
            else:
                SBS2_Basis = "BS3 in ClinGen VCEP."
                Functional_Result_Source = 'ClinGen'

    
        # Predictive evidence
        if function.isSBP2(var_type, SpliceAI_pred, dbscSNV_ada_score, dbscSNV_rf_score):
            standard = standard+'SBP2;'
            point -= 1
            SBP2 = 1
            SBP2_Basis = "A synonymous variant for which splicing prediction algorithms predict no impact on splicing."

        if isOVS1:
            standard = standard+'OVS1;'
            point += 8
            OVS1 = 1
            OVS1_Basis = consequence + " in tumor suppressor gene"           
       


        
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
                  'OS2_evidence': OS2_Basis,
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
                  'SBP2_Basis': SBP2_Basis
                  }

        df_output = pd.DataFrame([result])
        if index == 0:
            # If it is the first variant, write it directly to the file.
            if output_type == 'json':
                df_output.to_json(svoc_output_json_file, orient='records', lines=True, mode='w')
                print(f"SVOC annotation result has been successfully stored in {svoc_output_json_file}.")
            elif output_type == 'csv':
                df_output.to_csv(svoc_output_csv_file, index=False, mode='w')
                print(f"SVOC annotation result has been successfully stored in {svoc_output_csv_file}.")
            else:
                df_output.to_csv(svoc_output_txt_file, sep='\t', index=False, mode='w')
                print(f"SVOC annotation result has been successfully stored in {svoc_output_txt_file}.")
        else:
            # If it's not the first variant, append the written file.
            if output_type == 'json':
                df_output.to_json(svoc_output_json_file, orient='records', lines=True, mode='a')
                print(f"SVOC annotation result has been successfully stored in {svoc_output_json_file}.")
            elif output_type == 'csv':
                df_output.to_csv(svoc_output_csv_file, index=False, mode='a', header=False)
                print(f"SVOC annotation result has been successfully stored in {svoc_output_csv_file}.")
            else:
                df_output.to_csv(svoc_output_txt_file, sep='\t', mode='a', index=False)
                print(f"SVOC annotation result has been successfully stored in {svoc_output_txt_file}.")
        
    
    print(f"The score for {Variant_Information} is {point}, and the oncogenicity classification is {classification}.")
    print("%s" %end_description)

if __name__ == "__main__":
    main()