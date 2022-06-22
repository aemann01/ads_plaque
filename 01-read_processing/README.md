# Read processing

### 1. Download raw fastq files from /////

```bash
mkdir raw && cd raw
wget -i //////
```

### 2. Quality filter and trim raw reads

First need to generate some basic quality metrics for each file

```bash
ls *fastq.gz | parallel 'fastqc {}'
```

Filter by quality score (Q-30), remove TruSeq adapters, poly G tails. Since there is so much data this will take a long while. On a linux computer with 8 cores and 64G of ram this took about a weekend to complete.

```bash
mkdir cutadapt
ls *_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | parallel 'cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --nextseq-trim=20 -o cutadapt/{}.1.trim.fastq.gz -p cutadapt/{}.2.trim.fastq.gz --trim-n --minimum-length 100 --max-n 0 -q 30,30 {}_R1_001.fastq.gz {}_R2_001.fastq.gz 1>cutadapt/{}.trim.out'
```

NOTE: At this point you might want to remove raw reads to free up space

### 3. Map ribosomal RNA sequences

First need to download a reference database of LSU and SSU 16S and 18S rRNA sequences (here using SILVA v.138) to detect ribosomal RNAs.

```bash
mkdir ../../ribosomalRNA && cd ../../ribosomalRNA
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gzip -d SILVA*
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > silva_rRNA.fa
rm SILVA*
```

Index database

```bash
bowtie2-build silva_rRNA.fa silva_rRNA.db
```

Align sequences with bowtie2

```bash
cd ../raw/cutadapt
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 50 --gnu 'bowtie2 -x ../../ribosomalRNA/silva_rRNA.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal --no-head --no-sq -t -S ../../ribosomalRNA/{}.sam 2>../../ribosomalRNA/{}.out'
```

Get sequence identifiers that map to SILVA, remove sam files

```bash
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
rm *sam 
```

### 5. Map human sequences

```bash
cd ..
mkdir human_map && cd human_map
```

Download the human reference genome

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz
gzip -d GRCh38.p13.genome.fa.gz
```

Index the genome

```bash
bowtie2-build GRCh38.p13.genome.fa GRCh38.p13.genome.db
```

Align your reads

```bash
cd ../raw/cutadapt/
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ../../human_map/GRCh38.p13.genome.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal --no-head --no-sq -t -S ../../human_map/{}.sam 2>../../human_map/{}.out 1>../../human_map/{}.err'
```

Get sequence identifiers that map

```bash
ls *sam | while read line; do awk '{print $1}' $line | sort | uniq > $line.ids; done
rm *sam
```

### 6. Remove rRNA and human sequences

Concatenate ids to filter

```bash
mkdir ../raw/filtered
ls *ids | sed 's/.sam.ids//' | parallel 'cat {}.sam.ids ../ribosomalRNA/{}.sam.ids | sort | uniq > ../raw/filtered/{}.filt.ids'
```

Remove unwanted reads

```bash
cd ../raw/filtered/
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 32 'python ../../../00-scripts/remove_seqs.py -f ../cutadapt/{}.1.trim.fastq.gz -i {}.filt.ids -o {}.1.fastq'
ls *ids | sed 's/.filt.ids//' | parallel --gnu -j 32 'python ../../../00-scripts/remove_seqs.py -f ../cutadapt/{}.2.trim.fastq.gz -i {}.filt.ids -o {}.2.fastq'
```

### 8. OPTIONAL Test BWA mapping

```bash
cd ../03-human_map
bwa index GRCh38.p13.genome.fa
bwa mem -t 60 GRCh38.p13.genome.fa ../04-filtered/21PD_S7.1.fastq.gz ../04-filtered/21PD_S7.2.fastq.gz > 21PD_test_bwa.sam
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
```