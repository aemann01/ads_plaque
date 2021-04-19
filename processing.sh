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
