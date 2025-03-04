# SVOC
SVOC (Somatic Variant Oncogenicity Classifier) is a bioinformatics software tool for the classification of oncogenicity of somatic variants.

## PREREQUISITE

1. You need install Python version >= 3.11.5.
2. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) (version >= 2016-02-01) and [TransVar](https://github.com/zwdzwd/transvar) (version >= 2.5.10.20211024).
3. SVOC uses [AutoPVS1](https://github.com/JiguangPeng/autopvs1) to determine the effect of variants (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single-exon or multiexon deletion) on gene function. Although the AutoPVS1 algorithm is integrated into the SVOC project, the required resources for its execution (such as VEP, its cache, and FASTA files) need to be pre-installed.

    - VEP Installation
    ```bash
    cd SVOC/modules/
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
table_annovar = /mnt/annovar/table_annovar.pl
annovar_database_locat = /mnt/annovar/humandb
# the database location/dir from annnovar   check if database file exists
database_names = refGene avsnp151 dbnsfp47a_interpro clinvar_20240917 intervar_20180118 dbnsfp47a dbscsnv11
```
`ref_fasta` is the default reference genome file directory, or you can customize it. You need to download:    

- **hg19.fa** can be downloaded from UCSC [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and indexed with `samtools faidx`

- **hg38.fa** can be downloaded from NCBI [GRCh38_no_alt_analysis_set.fna.gz](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) and indexed with `samtools faidx`

`gnomadv2` is the default gnomad v2.1.1 file directory, or you can customize it. You need to download: 

- **gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz** can be downloaded from [genomes.liftover_grch38.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz) and TBI from [genomes.liftover_grch38.vcf.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi)

- **gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz** can be downloaded from [exomes.liftover_grch38.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz) and TBI from [exomes.liftover_grch38.vcf.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi)

- **gnomad.genomes.r2.1.1.sites.vcf.bgz** can be downloaded from [genomes.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz) and TBI from [genomes.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi)

- **gnomad.exomes.r2.1.1.sites.vcf.bgz** can be downloaded from [exomes.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz) and TBI from [exomes.bgz.tbi](https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi)

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
./svoc.py -c config.ini   # Run the examples in config.ini
./svoc.py -i your_vcf_input -o prefix_of_your_output
```

## HOW DOES IT WORK

SVOC accepts unannotated input files in VCF format, where each line corresponds to one somatic variant. SVOC will call ANNOVAR and TransVar to generate necessary annotations. Following annotation, the system automatically scores the variants based on the 17 rules stipulated by the ClinGen/CGC/VICC 2022 standards. Variants are then classified into one of five categories according to their scores: Oncogenic, Likely Oncogenic, Variant of Uncertain Significance (VUS), Likely Benign, or Benign. The output file can be formatted as JSON, CSV, or TXT, with each line corresponding to one somatic variant.

## Web server
Under development, please stay tuned. ^_^

## LICENSE

SVOC is free for non-commercial use without warranty. Users need to obtain licenses such as SpliceAI, TransVar, AutoPVS1, VEP and ANNOVAR by themselves. Please contact the authors for commercial use.

## REFERENCE
[The ClinGen/CGC/VICC 2022 standards](https://pubmed.ncbi.nlm.nih.gov/36063163/) Horak P, Griffith M, Danos AM, et al. Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity): Joint recommendations of Clinical Genome Resource (ClinGen), Cancer Genomics Consortium (CGC), and Variant Interpretation for Cancer Consortium (VICC). Genet Med. 2022;24(9):1991. doi:10.1016/j.gim.2022.07.001
