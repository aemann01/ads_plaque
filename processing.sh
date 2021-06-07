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

####################
# MAP RIBOSOMAL RNAS
####################
mkdir ribosomalRNA && cd ribosomalRNA
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > silva_rRNA.fa
rm SILVA*
# index database
bowtie2-build silva_rRNA.fa silva_rRNA.db
# align with bowtie2
cd ../cutadapt
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 50 --gnu 'bowtie2 -x ../ribosomalRNA/silva_rRNA.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal --no-head --no-sq -t -S ../ribosomalRNA/{}.sam 2>../ribosomalRNA/{}.out'
# get read ids to remove
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done

#################
# MAP HUMAN READS
#################
cd ..
mkdir human_map && cd human_map
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz
gzip -d GRCh38.p13.genome.fa.gz
# index database
bowtie2-build GRCh38.p13.genome.fa GRCh38.p13.genome.db
# align with bowtie2
cd ../cutadapt
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ../human_map/GRCh38.p13.genome.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal --no-head --no-sq -t -S ../human_map/{}.sam 2>../human_map/{}.out 1>../human_map/{}.err'
# get read ids to remove
cd ../human_map
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done

##############
# FILTER READS
##############
cd ../cutadapt
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 32 'python ../remove_seqs.py -f {}.1.trim.fastq.gz -i {}.filt.ids -o ../filtered/{}.1.fastq'
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 32 'python ../remove_seqs.py -f {}.2.trim.fastq.gz -i {}.filt.ids -o ../filtered/{}.2.fastq'

#############
# MERGE READS
#############
# merge paired end reads
ls *.1.* | sed 's/.1.fastq.gz//' | while read line; do pear -f $line.1.fastq.gz -r $line.2.fastq.gz -o /home/allie/uf_rnaseq/merged/$line.merge.fastq.gz -q 30 -j 7 1> /home/allie/uf_rnaseq/merged/$line.out; done

# test map with bwa (see if same as Alejandro's data)
bwa index GRCh38.p13.genome.fa
bwa mem -t 60 /home/allie/chelsea_data/human_map/GRCh38.p13.genome.fa /home/allie/chelsea_data/filtered/21PD_S7.1.fastq.gz /home/allie/chelsea_data/filtered/21PD_S7.2.fastq.gz > 21PD_test_bwa.sam
samtools view -S -b 21PD_test_bwa.sam > 21PD_test_bwa.bam
# get unmapped reads
samtools view -u -f 4 -F264 21PD_test_bwa.bam > 21PD_test_bwa.bam.temp1
samtools view -u -f 8 -F260 21PD_test_bwa.bam > 21PD_test_bwa.bam.temp2
samtools view -u -f 12 -F256 21PD_test_bwa.bam > 21PD_test_bwa.bam.temp3
# merge together
samtools merge -u 21PD_test_bwa.unmap.bam 21PD_test_bwa.bam.temp1 21PD_test_bwa.bam.temp2 21PD_test_bwa.bam.temp3
# sort
samtools sort -n 21PD_test_bwa.unmap.bam -o 21PD_test_bwa.sort.bam
# bam to fastq
bamToFastq -i 21PD_test_bwa.sort.bam -fq 21PD_test_bwa.sort.1.fq -fq2 21PD_test_bwa.sort.2.fq 

#################
# HUMANN ANALYSIS
#################
# install
conda create --name biobakery3 python=3.7
conda activate biobakery3
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
conda install humann -c biobakery
# extras for bowtie2
sudo apt-get install libtbb2 
sudo apt-get install libtbb2 
# download databases
humann3_databases --download chocophlan full /home/allie/refdb/humann3_db --update-config yes &
humann3_databases --download uniref uniref90_diamond /home/allie/refdb/humann3_db --update-config yes &
humann3_databases --download utility_mapping full /home/allie/refdb/humann3_db --update-config yes 
# update diamond
conda install -c bioconda diamond=0.9.36
# test installation
humann_test
humann3 -i demo.fastq.gz -o demo_humann3
# run humann3
ls *fastq.gz | sed 's/.fastq.gz//' | parallel -j 60 'humann3 -i {}.fastq.gz -o ../humann3/{} --remove-temp-output 1>../humann3/{}.out 2>../humann3/{}.err'
# join tables together
cd ..
humann_join_tables -i humann3/ -s -o humann3/merged.tsv  










############################################
# TESTING SUBSET DENVO VS REF BASED ASSEMBLY
############################################
cd ..
mkdir test_assembly && cd test_assembly
# denovo assembly with trinity
Trinity --seqType fq --SS_lib_type RF --left ../filtered/12PD_S9.2.fastq.gz,../filtered/13PF_S20.2.fastq.gz,../filtered/47PER_S3.2.fastq.gz --right ../filtered/12PD_S9.1.fastq.gz,../filtered/13PF_S20.1.fastq.gz,../filtered/47PER_S3.1.fastq.gz --CPU 32 --max_memory 10G
# get stats on the run
TrinityStats.pl Trinity.fasta
# what proportion of reads map back onto the assembled contigs?
# ################################
# ## Counts of transcripts, etc.
# ################################
# Total trinity 'genes':	802666
# Total trinity transcripts:	1172055
# Percent GC: 50.00

# ########################################
# Stats based on ALL transcript contigs:
# ########################################

# 	Contig N10: 11624
# 	Contig N20: 6087
# 	Contig N30: 3653
# 	Contig N40: 2253
# 	Contig N50: 1387

# 	Median contig length: 349
# 	Average contig: 751.62
# 	Total assembled bases: 880945314


# #####################################################
# ## Stats based on ONLY LONGEST ISOFORM per 'GENE':
# #####################################################

# 	Contig N10: 8249
# 	Contig N20: 4073
# 	Contig N30: 2301
# 	Contig N40: 1367
# 	Contig N50: 855

# 	Median contig length: 324
# 	Average contig: 612.35
# 	Total assembled bases: 491509239

# How many reads map to HOMD vs our denovo assembled data?
wget http://www.homd.org/ftp/genomes/PROKKA/current/tsv/ALL_genomes.tsv
wget http://www.homd.org/ftp/genomes/PROKKA/current/SEQID_info.txt
wget http://www.homd.org/ftp/genomes/PROKKA/current/fna/ALL_genomes.fna
# index bowtie2 databases
bowtie2-build ALL_genomes.fna ALL_genomes.db &
bowtie2-build trinity_out_dir/Trinity.fasta Trinity.db &
# get fasta index files for igv
samtools faidx ALL_genomes.fna
samtools faidx trinity_out_dir/Trinity.fasta
# map to transcriptome and homd -- which perfoms better?
bowtie2 -x Trinity.db -1 filtered/12PD_S9.1.fastq.gz -2 filtered/12PD_S9.2.fastq.gz -S test.sam --end-to-end --no-unal
bowtie2 -x ALL_genomes.db -1 filtered/12PD_S9.1.fastq.gz -2 filtered/12PD_S9.2.fastq.gz -S test_HOMD.sam --end-to-end --no-unal






# Testing humann v3.0
# download databases
humann_databases --download chocophlan full /home/allie/refdb/humann_db/ --update-config yes









