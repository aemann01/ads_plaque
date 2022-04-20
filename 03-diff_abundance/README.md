
```R
library(DESeq2)
library(Glimma)
# load data
metadata <- read.table("map.txt", header=T, sep="\t")
# only compare PD to PF
metadata <- metadata[metadata$tooth_type != "PE",]
row.names(metadata) <- metadata$UFsample_id
head(metadata)
# read in gene count matrix
genecounts <- read.delim("raw_counts.matrix.txt", header=T, sep="\t", row.names=1)
# only keep columns found in metadata
genecounts <- genecounts[, colnames(genecounts) %in% row.names(metadata)]
head(genecounts)
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = genecounts, colData = metadata, design = ~tooth_type)
# filter out any genes with fewer than 10 reads total
star_results <- star_results[rowSums(counts(star_results)) >= 10,]
star_results
# set factor level (this determines which direction the comparisions are made -- by default it's by alphabetical order)
star_results$tooth_type <- factor(star_results$tooth_type, levels=c("PD", "PF"))
# run deseq
res <- DESeq(star_results)
res <- results(res, alpha=0.01)
summary(res)
# write results to file
write.table(as.data.frame(res), file="deseq_results.csv", quote=F, sep="\t")
```










Now want to get gene function associated with each gene name

```bash
awk '{print $1}' deseq_results.csv | while read line; do grep -m 1 $line ../03-STAR/star_resuts/kraken_refgenomes/all_genomes.gtf ; done > test

```

HOMD doesn't generate any read count feature files, use featurecount to pull this info from the bam files and then can do deseq 

```bash
 conda install -c bioconda subread 


# remove any transcripts with no gene id
# how many entries have no gene id?
grep -v 'gene_id ""' kraken_refgenomes/all_genomes.gtf > kraken_refgenomes/all_genomes.w.geneid.gtf






featureCounts -a ../kraken_refgenomes/all_genomes.w.geneid.gtf -o readcounts -T 60 ../*Aligned.out.bam
```




















### Install interproscan database


```bash
# download interproscan database
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.55-88.0/interproscan-5.55-88.0-64-bit.tar.gz.md5
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.55-88.0/interproscan-5.55-88.0-64-bit.tar.gz
# checksum
md5sum -c interproscan-5.55-88.0-64-bit.tar.gz.md5
# untar gz
tar xvzf interproscan-5.55-88.0-64-bit.tar.gz
# remove old DB
rm -rf /home/allie/src/miniconda3/envs/uf_rnaseq/share/InterProScan/data/
# move new database
mv interproscan--/data /home/allie/src/miniconda3/envs/uf_rnaseq/share/InterProScan/
```

### Sort bam files by read name

```bash
ls *Aligned.out.bam | while read line; do samtools sort -n -o $line.sort -O bam -@ 60 $line 2>$line.sort.err -m 2G; done



```

Some bam files were truncated when transferring from palmetto with globus -- requeue transfer and checksum

```bash
md5sum 25PER_S1Aligned.out.bam
# b47fde40653893c2568363c8510a9949
md5sum 44PD_S13Aligned.out.bam
# 75db0d6dff43c7c484f8c3a89944f348
```

These are getting re mapped because of truncated bam files -- in the mean time just continue on with pipeline

### Get gene counts using htseq-count


















```bash
# test
grep "CDS" kraken_refgenomes/all_genomes.w.geneid.gtf | grep "gene_id" > kraken_refgenomes/test.gtf
# test feature count
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -M -C -B -O -a kraken_refgenomes/test.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam --extraAttributes gene -t CDS -g gene_id; done






```



grep -v 'gene_id ""' kraken_refgenomes/all_genomes.gtf | grep "gene_id" > kraken_refgenomes/all_genomes.w.geneid.gtf



### Get gene counts using feature counts

```bash
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -M -C -B -O -a kraken_refgenomes/all_genomes.w.geneid.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam --extraAttributes gene; done 
```


### Try FADU

```bash
git clone https://github.com/IGS/FADU.git
cd FADU


```












### 1. Get matrix of gene expression data from all samples
```bash
paste ../03-STAR/star_resuts/*ReadsPerGene.out.tab | grep -v "^N_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > temp
# add headers
sed -e "1igeneName\t$(ls ../03-STAR/star_resuts/*ReadsPerGene.out.tab  | sed 's/ReadsPerGene.out.tab//g' | sed 's/..\/03-STAR\/star_resuts\///' | sed 's/_.*//' | sed 's/^/UF/' |  tr '\n' '\t')" temp > raw_counts.matrix.txt
sed -i 's/\t$//' raw_counts.matrix.txt
rm temp
````

```R
library(DESeq2)
library(Glimma)
# load data
metadata <- read.table("map.txt", header=T, sep="\t")
# only compare PD to PF
metadata <- metadata[metadata$tooth_type != "PE",]
row.names(metadata) <- metadata$UFsample_id
head(metadata)
# read in gene count matrix
genecounts <- read.delim("raw_counts.matrix.txt", header=T, sep="\t", row.names=1)
# only keep columns found in metadata
genecounts <- genecounts[, colnames(genecounts) %in% row.names(metadata)]
head(genecounts)
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = genecounts, colData = metadata, design = ~tooth_type)
# filter out any genes with fewer than 10 reads total
star_results <- star_results[rowSums(counts(star_results)) >= 10,]
star_results
# set factor level (this determines which direction the comparisions are made -- by default it's by alphabetical order)
star_results$tooth_type <- factor(star_results$tooth_type, levels=c("PD", "PF"))
# run deseq
res <- DESeq(star_results)
res <- results(res, alpha=0.01)
summary(res)
# write results to file
write.table(as.data.frame(res), file="deseq_results.csv", quote=F, sep="\t")
```










Now want to get gene function associated with each gene name

```bash
awk '{print $1}' deseq_results.csv | while read line; do grep -m 1 $line ../03-STAR/star_resuts/kraken_refgenomes/all_genomes.gtf ; done > test

```

HOMD doesn't generate any read count feature files, use featurecount to pull this info from the bam files and then can do deseq 

```bash
 conda install -c bioconda subread 


# remove any transcripts with no gene id
# how many entries have no gene id?
grep -v 'gene_id ""' kraken_refgenomes/all_genomes.gtf > kraken_refgenomes/all_genomes.w.geneid.gtf






featureCounts -a ../kraken_refgenomes/all_genomes.w.geneid.gtf -o readcounts -T 60 ../*Aligned.out.bam
```













