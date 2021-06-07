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
mkdir ../cutadapt
ls *_R1_*.fastq.gz | sed 's/_R1_001.fastq.gz//' | parallel 'cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --nextseq-trim=20 -o ../cutadapt/{}.1.trim.fastq.gz -p ../cutadapt/{}.2.trim.fastq.gz --trim-n --minimum-length 100 --max-n 0 -q 30,30 {}_R1_001.fastq.gz {}_R2_001.fastq.gz 1>../cutadapt/{}.trim.out'
```

NOTE: At this point you might want to remove raw reads to free up space

### 3. Remove ribosomal RNA sequences

First need to download a reference database of LSU and SSU 16S and 18S rRNA sequences (here using SILVA v.138) to detect ribosomal RNAs.

```bash
mkdir ribosomalRNA && cd ribosomalRNA
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
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 50 --gnu 'bowtie2 -x ../ribosomalRNA/silva_rRNA.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal --no-head --no-sq -t -S ../ribosomalRNA/{}.sam 2>../ribosomalRNA/{}.out'
```


