<div align="center">    
    
![Teaser image](./figures/SVOC_LOGO.png)


</div>

## Table of Contents
- [PREREQUISITE](#PREREQUISITE)
- [OPTIONS](#OPTIONS)
- [EXAMPLE](#EXAMPLE)
- [HOW DOES IT WORK](#HOW-DOES-IT-WORK)
- [WEB SERVER](#WEB-SERVER)
- [LICENSE](#LICENSE)
- [REFERENCE](#REFERENCE)
- [CONTACT](#CONTACT)

## PREREQUISITE

1. Clone the repository:
 ```bash
 git clone https://github.com/leqingsang/SVOC.git
 cd SVOC
 ```

2. Create a conda environment and install the required dependencies:
```bash
conda create -n {ENV_NAME} python=3.11
conda activate {ENV_NAME}
pip install -r requirements.txt
```

3. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) (version >= 2016-02-01), [TransVar](https://github.com/zwdzwd/transvar) (version >= 2.5.10.20211024) and [BCFTools](https://github.com/samtools/bcftools).
   
4. SVOC references the results of loss-of-function variant judgments from [AutoPVS1](https://github.com/JiguangPeng/autopvs1). Although the AutoPVS1 algorithm is integrated into the SVOC project (`SVOC/modules/autopvs1`), the required resources for its execution (such as VEP, its cache, and FASTA files) need to be pre-installed.

    - VEP Installation
    ```bash
    cd /modules
    git clone https://github.com/Ensembl/ensembl-vep.git
    cd ensembl-vep
    git pull
    git checkout release/104
    perl INSTALL.pl
    ```

    - VEP cache and faste files
    VEP cache and faste files can be automatically downloaded and configured using [INSTALL.pl](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer). You can also download and set up them manually:

    ```bash
    r=104
    FTP='ftp://ftp.ensembl.org/pub/'

    # indexed vep cache
    cd $HOME/.vep
    wget $FTP/release-${r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${r}_GRCh38.tar.gz
    wget $FTP/release-${r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${r}_GRCh37.tar.gz
    tar xzf homo_sapiens_vep_${r}_GRCh37.tar.gz
    tar xzf homo_sapiens_vep_${r}_GRCh38.tar.gz

    # fasta
    cd $HOME/.vep/homo_sapiens_refseq/${r}_GRCh37/
    wget $FTP/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    tar xzf Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

    cd $HOME/.vep/homo_sapiens_refseq/${r}_GRCh38/
    wget $FTP/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    tar xzf Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    ```
4. Configuration

`SVOC/config.ini`

```ini
[SVOC]
buildver = hg19 
inputfile = example/example1.vcf
outfile = example/example1
output_type = txt
database_svoc = svocdb 
mane = %(database_svoc)s/MANE.GRCh38.v1.4.ensembl_genomic.gtf
ref_fasta = %(database_svoc)s/ref_fasta
gnomadv2 = %(database_svoc)s/gnomadv2

[Annovar]
table_annovar = /annovar/table_annovar.pl
annovar_database_locat = /annovar/humandb
# check if database file exists
database_names = refGene dbnsfp47a dbscsnv11
```
`ref_fasta` is the default directory for reference genome file, or you can customize it. You need to download:    

- **hg19.fa** can be downloaded from UCSC [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and indexed with `samtools faidx`

- **hg38.fa** can be downloaded from NCBI [GRCh38_no_alt_analysis_set.fna.gz](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) and indexed with `samtools faidx`

- **Note:** The reference genome file also needs to be placed in `SVOC/modules/autopvs1/data` to meet the needs of autopvs1 running.

`gnomadv2` is the default directory for gnomad v2.1.1 file, or you can customize it. You need to download: 

- **gnomad.genomes.r2.1.1.sites.vcf.bgz** can be downloaded from [genomes.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz) and TBI from [genomes.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi)

- **gnomad.exomes.r2.1.1.sites.vcf.bgz** can be downloaded from [exomes.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz) and TBI from [exomes.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi)

- **gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz** can be downloaded from [genomes.liftover_grch38.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz) and TBI from [genomes.liftover_grch38.vcf.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi)

- **gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz** can be downloaded from [exomes.liftover_grch38.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz) and TBI from [exomes.liftover_grch38.vcf.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi)


## OPTIONS

- --version             
show program''s version number and exit

- -h, --help              
show this help message and exit  

- -c config.ini, --config=config.ini           
your own configure file           
You can edit all the options in the configure and if you use this options, you can ignore all the other options bellow.

- -b hg19, --buildver=hg19    
genomic build version           
It can be hg19 and hg38, will support other version later.

- -i example/example1.vcf, --input=example/example1.vcf           
input file of variants for classification           
It can be VCF format, will support other famate later.

- -o example1, --output=example1     
prefix of the output file

- --output_type=txt     
type of the output file     
It can be JSON, CSV and TXT(default) format.

- -d svocdb, --database_svoc=svocdb     
database location/dir for the SVOC dataset files

- --table_annovar=./table_annovar.pl     
Annovar perl script of table_annovar.pl

- --annovar_database_locat=humandb     
database location/dir for the Annovar annotation datasets


## EXAMPLE
```
python svoc.py -c config.ini   # Run the examples in config.ini
python svoc.py -i your_vcf_input -o prefix_of_your_output
```

## HOW DOES IT WORK

SVOC accepts unannotated input files in VCF format, where each line corresponds to one somatic variant. SVOC will call ANNOVAR and TransVar to generate necessary annotations. Then, it automatically scores all 17 criteria and calculates the total score of the variant based on the applicable rules defined in the guidelines. The oncogenicity category is finally assigned according to the total score: ≥10 points for Oncogenic, 6–9 points for Likely Oncogenic, 0–5 points for VUS, -1 to -6 points for Likely Benign, and ≤-7 points for Benign. The output file can be formatted as JSON, CSV, or TXT, with each line corresponding to one somatic variant.

## WEP SERVER
https://svoc.premedkb.cn

## LICENSE

SVOC is free for non-commercial use without warranty. Users need to obtain licenses such as SpliceAI, TransVar, AutoPVS1, VEP and ANNOVAR by themselves. Please contact the authors for commercial use.

## REFERENCE
Horak P, Griffith M, Danos AM, et al. [Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity): Joint recommendations of Clinical Genome Resource (ClinGen), Cancer Genomics Consortium (CGC), and Variant Interpretation for Cancer Consortium (VICC)](https://pubmed.ncbi.nlm.nih.gov/36063163/). <i>Genet Med</i>. 2022;24(9):1991. [doi:10.1016/j.gim.2022.07.001](https://pubmed.ncbi.nlm.nih.gov/36063163/)

## CONTACT
For technical questions please open issue, or contact:
- Leqing Sang <lqsang25@m.fudan.edu.cn>
