### Post mapping processing

#### Re run featurecounts with homd mapped data with transcript_id and CDS as target

```bash
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -M -C -B -O -a ALL_genomes.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t CDS -g transcript_id; done
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













