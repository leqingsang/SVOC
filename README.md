# SVOC
SVOC (Somatic Variant Oncogenicity Classifier) is a bioinformatics software tool for the classification of oncogenicity of somatic variants.

## PREREQUISITE

1. You need install Python version >= 3.11.5.
2. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) (version >= 2016-02-01) and [TransVar](https://github.com/zwdzwd/transvar) (version >= 2.5.10.20211024).
3. SVOC uses [AutoPVS1](https://github.com/JiguangPeng/autopvs1) to determine the effect of variants (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single-exon or multiexon deletion) on gene function. Although the AutoPVS1 algorithm is integrated into the SVOC project, the required resources for its execution (such as VEP, its cache, and FASTA files) need to be pre-installed.

VEP Installation
```bash
cd SVOC/modules/
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/104
perl INSTALL.pl
```

VEP cache and faste files
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
