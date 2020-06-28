# Quality check of the reads got from the company

mkdir data
mv *fq.gz data/
cd data
mkdir fastqc
~/softwares/FastQC/fastqc -t 10 -o fastqc *fq.gz

# Trimming adapters and low quality base using the software Trimmomatic-0.38

mkdir trimmed
cd trimmed
cp ~/softwares/Trimmomatic-0.38/adapters/TruSeq3-PE* .
sh trim.sh >> trimlog 2>&1
    # Quality check of the reads trimmed
    mkdir trimmed_qc
    ~/softwares/FastQC/fastqc -t 10 -o trimmed_qc *fq.gz

# Download genome file from PLAZA

mkdir genome
cd genome
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/Genomes/peq.con.gz
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/GO/go.peq.csv.gz
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/GFF/peq/annotation.all_transcripts.exon_features.peq.gff3.gz
wget ftp://ftp.psb.ugent.be/pub/plaza/plaza_public_monocots_04/GFF/peq/annotation.all_transcripts.all_features.peq.gff3.gz
gunzip peq.con.gz
mv peq.con pha.fa

# Mapping the trimmed reads to the orchid genome

hisat2-build pha.fa pha
mkdir mapping
cd mapping
conda activate python27;
sh hisat.sh  >> mappinglog 2>&1

# Prepare the gene count table for DESeq2
mkdir count;cd count;conda activate base;sh stringtie.sh
for i in sl1 sl2 sl3 sf1 sf2 sf3;do echo $i;grep "FPKM" ${i}/${i}.gtf |sed 's#gene_id "##' |sed 's#"; transcript_id "#\t#'|sed 's#"; ref_gene_name "AT[0-9]G[0-9]*"; cov "#\t#'|sed 's#"; FPKM "#\t#'|sed 's#"; TPM "#\t#
'|sed 's#";$##'|sed 's#"; cov "#\t#'|sed 's#transcript_id "#.\t#' > ${i}_cov_fpkm_tpm;done
wget http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
ls ./*/*gtf > sample  # vi sample
conda activate python27
python prepDE.py -i sample
for i in sl1 sl2 sl3 sf1 sf2 sf3;do echo $i;sort ${i}.gene.abundance -o ${i}.gene.abundance
join -1 1 -2 1 --header o sl3.gene.abundance > oo
join -1 1 -2 1 --header oo sf1.gene.abundance > ooo
join -1 1 -2 1 --header ooo sf2.gene.abundance > oooo
join -1 1 -2 1 --header oooo sf3.gene.abundance > all.sample.gene.abundance

# DEG analysis using DESeq2 in R
