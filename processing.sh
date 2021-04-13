#############
# DEMULTIPLEX
#############
cd uf_novaseq_rnaseq_pool1
bcl2fastq -o Data/Intensities --no-lane-splitting
cd uf_novaseq_rnaseq_pool2
bcl2fastq -o Data/Intensities --no-lane-splitting

####################
# SET UP ENVIRONMENT
####################
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
ls *_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | parallel 'cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --nextseq-trim=20 -o cutadapt/{}.1.trim.fastq.gz -p cutadapt/{}.2.trim.fastq.gz --trim-n --minimum-length 100 --max-n 0 -q 30,30 {}_R1_001.fastq.gz {}_R2_001.fastq.gz 1>cutadapt/{}.trim.out'

#######################
# REMOVE RIBOSOMAL RNAS
#######################
mkdir ribosomalRNA && cd ribosomalRNA
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > silva_rRNA.fa
rm SILVA*
# index database
bowtie2-build silva_rRNA.fa silva_rRNA.db
# align with bowtie2
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel --gnu 'bowtie2 -x ../ribosomalRNA/silva_rRNA.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter -S ../ribosomalRNA/{}.sam 2>../ribosomalRNA/{}.out'
ls *sam | sed 's/.sam//' | while read line; do samtools view -f 4 $line.sam > $line.unmap.bam; done
# clean up files for space
rm ../cutadapt/*fastq.gz 
ls *sam | sed 's/.sam//' | parallel --gnu 'samtools view -S -b {}.sam > {}.bam'
# convert bam to fastq
ls *unmap.bam | sed 's/.unmap.bam//' | parallel --gnu 'bedtools bamtofastq -i {}.unmap.bam -fq {}.1.filt.fq -fq2 {}.2.filt.fq' 






















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
# STAR --runMode genomeGenerate \
# 	--genomeFastaFiles GRCh38.p13.genome.fa  \
# 	--runThreadN 2 \
# 	--sjdbGTFfile gencode.v36.annotation.gtf 




# map to genome
STAR --runThreadN 8 \
	--runRNGseed 549 \
	--readFilesCommand zcat \
	--readFilesIn read1.fq.gz read2.fq.gz \
	--genomeDir GenomeDir \
	--outSAMtype BAM \
	--outReadsUnmapped {}.unmapped.bam 

# convert bam to sam

# filter out reads that map to human from dataset







#################
# DENOVO ASSEMBLY
#################
# merge all left and right files together into two separate files to generate reference assembly
//////
# generate denovo assembly 
Trinity --seqType fq --SS_lib_type RF --left --right  --CPU 4 -o trinity_out --min_kmer_cov 2 1>trinity.out 2>trinity.err

##################################
# REFERENCE BASED ALIGNMENT (HOMD)
##################################
# download homd genomes for reference
mkdir homd_map && cd homd_map
wget http://www.homd.org/ftp/HOMD_prokka_genomes/gff/ALL_genomes.gff &
wget http://www.homd.org/ftp/HOMD_prokka_genomes/fna/ALL_genomes.fna 
rm wget-log*
# convert GFF to GTF format (only run once)
gffread ALL_genomes.gff -T -o ALL_genomes.gtf
# index genomes (only run once)
# STAR --runMode genomeGenerate \
# 	--genomeFastaFiles ALL_genomes.fna \
# 	--runThreadN 8 \
# 	--limitGenomeGenerateRAM 66959267424 \
# 	--sjdbGTFfile ALL_genomes.gtf \
# 	--genomeChrBinNbits 15
# map to HOMD genome database






STAR --runThreadN 4 \
	--genomeDir GenomeDir \
	--readFilesIn {}.fastq \
	--outFileNamePrefix {} \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped {}.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--runRNGseed 549 \
	--chimOutType SeparateSAMold 


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








