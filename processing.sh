#############
# DEMULTIPLEX
#############
cd uf_novaseq_rnaseq_pool1
bcl2fastq -o Data/Intensities --no-lane-splitting
cd uf_novaseq_rnaseq_pool2
bcl2fastq -o Data/Intensities --no-lane-splitting

###############
# DOWNLOAD HOMD
###############
# download homd genomes for reference
mkdir refdb && cd refdb
wget http://www.homd.org/ftp/HOMD_prokka_genomes/gff/ALL_genomes.gff &
wget http://www.homd.org/ftp/HOMD_prokka_genomes/fna/ALL_genomes.fna 
rm wget-log*

####################
# SET UP ENVIRONMENT
####################
cd ..
conda env create -f environment.yml
conda activate uf_rnaseq
# to deactivate:
# conda deactivate 

#########################
# QUALTIY FILTER AND TRIM
#########################
# first generate quality score metrics
cd raw
ls *fastq.gz | parallel 'fastqc {}'
# trim low quality sequences
mkdir cutadapt
ls *_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | parallel 'cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o cutadapt/{}.1.trim.fastq.gz -p cutadapt/{}.2.trim.fastq.gz --trim-n --minimum-length 100 --max-n 0 -q 30,30 {}_R1_001.fastq.gz {}_R2_001.fastq.gz 1>cutadapt/{}.trim.out'










## DON'T MERGE??




########################
# FILTER OUT HUMAN READS
########################
# get human genome and GTF file
mkdir human_map && cd human_map
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz &
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz 
gzip -d *gz
rm wget-log
# index human genome -- NOTE: this takes a long time but only needs to be run once
STAR --runMode genomeGenerate \
	--genomeFastaFiles GRCh38.p13.genome.fa  \
	--runThreadN 2 \
	--sjdbGTFfile gencode.v36.annotation.gtf 
# map to genome
# map to reference genomes
STAR --runThreadN 4 \
	--genomeDir GenomeDir \
	--readFilesIn {}.fastq \
	--outFileNamePrefix {} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped {}.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold 
# convert bam to sam

# filter out reads that map to human from dataset

########################
# REMOVE RIBOSOMAL GENES
########################





#################
# DENOVO ASSEMBLY
#################
# generate denovo assembly 
Trinity --seqType fq --single /////.fq --CPU 4 -o trinity_out

##################################
# REFERENCE BASED ALIGNMENT (HOMD)
##################################

# how many reads retained when you map merged reads back? go with that approach






# convert GFF to GTF format
gffread ALL_genomes.gff -T -o ALL_genomes.gtf
# generate reference index 
# NOTE:
# be careful to check no one is running something else when running this command or the computer will crash
# NOTE:
# this step takes a long time but only has to be run once
STAR --runMode genomeGenerate \
	--genomeFastaFiles ALL_genomes.fna  \
	--runThreadN 8 \
	--limitGenomeGenerateRAM 66959267424 \
	--sjdbGTFfile ALL_genomes.gtf \
	--genomeChrBinNbits 15
# map to reference genomes
STAR --runThreadN 4 \
	--genomeDir GenomeDir \
	--readFilesIn {}.fastq \
	--outFileNamePrefix {} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped {}.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold 




#################
# GET GENE COUNTS
#################
featureCounts -a ALL_genomes.gtf -F GTF -o countMatrix.txt star-resultsAligned.sortedByCoord.out.bam -g transcript_id








