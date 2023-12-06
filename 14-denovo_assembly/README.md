# Denovo assembly analysis 

## Sample by sample full denovo assembly

Since cannot meet a Trinity assembly checkpoint with full denovo assembly within the walltime constraints on palmetto with all the data merged together, going to go with a tiered approach where we:

(1) assemble sample by sample
(2) merge the assemblies into a single reference
(3) dereplicate full length transcripts to cut down on multimapping from same strain 
(4) annotate to get transcripts with full ADS operon
(5) map individual samples to these references ONLY

First do a sample by sample assembly 

```bash
#PBS -N denovo_SAMP
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/amann3/ads_raw_data/filtered/

# load trinity module
module add trinityrnaseq/2.15.1

# run
   Trinity \
   --seqType fq \
   --max_memory 750G \
   --left SAMP.1.fastq.gz \
   --right SAMP.2.fastq.gz \
   --CPU 40 \
   --output /scratch1/amann3/trinity_denovo_SAMP \
   --full_cleanup
```

From above example pbs script, generate for each sample

```bash
cat sample.list | while read line; do sed "s/SAMP/$line/g" example.pbs > $line.pbs; done
```

After sample wise assembly has completed, merge all into one file, dereplicate with vsearch, run prokka

```bash
#PBS -N prokka-denovo
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l walltime=72:00:00
#PBS -j oe

# move to folder with data
cd /scratch1/amann3/complete_assemblies

# load modules
module add prokka/1.14.5 vsearch/2.21.1

# merge and dereplicate
cat trinity*fasta > all_assemblies.fasta
vsearch --derep_fulllength all_assemblies.fasta \
   --output all_assemblies.uniq.fasta
echo "all transcripts"
grep ">" all_assemblies.fasta
echo "unique transcripts"
grep ">" all_assemblies.uniq.fasta

# run prokka
prokka --prefix all_assemblies  all_assemblies.uniq.fasta

echo "prokka complete"
```

Pull any transcripts with arcABC

```bash
grep -E "eC_number=3.5.3.6|eC_number=2.1.3.3|eC_number=2.7.2.2" all_assemblies_bigmem.gff | awk '{print $1}' | sort | uniq -c > ads_transcript.ids

# how many truncated operons? i.e., only one or two genes included in the same transcript?
sed -i 's/^[ \t]*//g' ads_transcript.ids
# how many with only a single gene 
grep "^1" ads_transcript.ids | wc -l
# 7792
# pull partial or non truncated
grep -v "^1" ads_transcript.ids | awk '{print $2}' > ads_transcripts_nontrunc.ids
module load seqtk/1.3-r106
seqtk subseq all_assemblies_bigmem.fna ads_transcripts_nontrunc.ids > ads_transcripts_nontrunc.fna
```

Rebuild gff file from this and then generate star reference database

```bash
# this complete very quickly so not work messing with filtering the gff
#PBS -N prokka-denovo
#PBS -l select=1:ncpus=20:mem=300gb
#PBS -l walltime=72:00:00
#PBS -j oe

# move to folder with data
cd /scratch1/amann3/complete_assemblies

# load modules
module add prokka/1.14.5

# run prokka
prokka --prefix ads_assemblies \
   --force \
   --cpus 20 \
   --norrna \
   --notrna \
   all_assemblies_bigmem/ads_transcripts_nontrunc.fna

echo "prokka complete"
```

Now convert GFF to GTF for star

working directory on palmetto: /home/amann3/uf_rnaseq/14-denovo_assembly/samp_assembly/star_mapping

```bash
#PBS -N gffread_ads
#PBS -l select=1:ncpus=10:mem=50gb
#PBS -l walltime=3:00:00
#PBS -j oe

# move to folder with data
cd /scratch1/amann3/complete_assemblies/ads_assemblies/

# load module
module add gffread/0.12.7

# run
gffread -T ads_assemblies.gff -o ads_assemblies.gtf

echo "gffread complete"
```

Generate reference database for STAR

```bash
#PBS -N ads-genomeGenerate
#PBS -l select=1:ncpus=10:mem=750gb
#PBS -l walltime=72:00:00
#PBS -q bigmem

# load module
module add star/2.7.5c

# move into scratch
cd /scratch1/amann3/complete_assemblies/ads_assemblies/

# run command
STAR --runMode genomeGenerate \
   --genomeFastaFiles /scratch1/amann3/complete_assemblies/ads_assemblies/ads_assemblies.fna  \
   --runThreadN 10 \
   --limitGenomeGenerateRAM 228697571594 \
   --sjdbGTFfile /scratch1/amann3/complete_assemblies/ads_assemblies/ads_assemblies.gtf \
   --genomeChrBinNbits 10 \
   --limitSjdbInsertNsj 3204994 \
   --genomeSAindexNbases 10\
   --sjdbGTFfeatureExon CDS
```

Map raw reads to contigs (serial -- takes ~ 3hrs per sample)

```bash
#PBS -N STAR-map_ads
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=168:00:00
#PBS -j oe
#PBS -q bigmem

# load module
module add star/2.7.5c

# move into scratch
cd /scratch1/amann3/ads_raw_data/filtered

# run command
ls *.1.*fastq.gz | sed 's/.1.fastq.gz//' | while read line; do STAR \
   --runThreadN 40 \
   --genomeDir /scratch1/amann3/complete_assemblies/ads_assemblies/GenomeDir \
   --readFilesIn <(gunzip -c $line\.1.fastq.gz) <(gunzip -c $line\.2.fastq.gz) \
   --outFileNamePrefix $line.ads_denovo \
   --outSAMtype BAM Unsorted \
   --quantMode TranscriptomeSAM GeneCounts \
   --alignIntronMax 1 \
   --chimOutType SeparateSAMold ; done
```

Sample specific star mapping (test -- if this works, run instead of above)

```bash
cat samps.ids | while read line; do sed "s/SAMP/$line/g" example_star.pbs > $line.star_map.pbs; done
```

Download resulting bam files to hillary, run feature counts on output

```bash
cd ~/ads_plaque/14-denovo_assembly/samp_assembly/star_map_palmetto
scp amann3@login.palmetto.clemson.edu:/scratch1/amann3/ads_raw_data/filtered/\*bam .
scp amann3@login.palmetto.clemson.edu:/scratch1/amann3/complete_assemblies/ads_assemblies/ads_assemblies.gtf .
conda activate 2022-ADS_plaque

mkdir featurecounts
ls *ads_denovoAligned.out.bam | sed 's/.ads_denovoAligned.out.bam//' | while read line; do featureCounts -f -p -C -B -M -a ads_assemblies.gtf \
   -o featurecounts/$line.out \
   -T 16 $line.ads_denovoAligned.out.bam \
   -t CDS -g transcript_id; done

# some of these fail individually for some reason, try to run individually
featureCounts -f -p -C -B -M -a ads_assemblies.gtf \
   -o featurecounts/87PFR_S21.out \
   87PFR_S21.ads_denovoAligned.out.bam \
   -t CDS -g transcript_id

# segmentation fault, the bam files are corrupted, re run on palmetto with a ton of cores to try and finish up...
# once corrupted bam files are reprocessed, continue to combine output

# combine output
cd featurecounts
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt

# clean up sample names
sed 's/.ads_denovoAligned.out.bam//g' read_counts.txt | sed 's/_S..\t/\t/'g | sed 's/_S.\t/\t/g' > temp
mv temp read_counts.txt
```

Now want to run deseq analysis on these results 

```bash
## RUN WITH CONDA ENVIRONMENT ON LOBSANG
cd /Users/mann/github/ads_plaque/14-denovo_assembly
conda activate 2022-ADS_plaque-R
scp allie@hillary.clemson.edu:/home/allie/ads_plaque/14-denovo_assembly/samp_assembly/star_map_palmetto/featurecounts/read_counts.txt .
R
```

```R
library(pheatmap, warn.conflicts = F, quietly = T)
library(ggplot2, warn.conflicts = F, quietly = T)
library(limma, warn.conflicts = F, quietly = T)
library(DESeq2, warn.conflicts = F, quietly = T)
library(edgeR, warn.conflicts = F, quietly = T)
library(Glimma, warn.conflicts = F, quietly = T)
library(apeglm, warn.conflicts = F, quietly = T)

# load data
metadata <- read.table("map.txt", header=T, sep="\t")
row.names(metadata) <- metadata$sample_id
genecounts <- read.table("read_counts.txt", header=T, sep="\t", row.names=1)
colnames(genecounts) <- gsub(x = names(genecounts), pattern = "X", replacement = "UF")  

# format for deseq analysis
# only compare PD to PF
metadata <- metadata[metadata$tooth_type != "PE",]
# only keep columns found in metadata
genecounts <- genecounts[, colnames(genecounts) %in% row.names(metadata)]
metadata <- metadata[row.names(metadata) %in% colnames(genecounts),]

# reorder by metadata rownames
metadata <- metadata[order(colnames(genecounts)),]
# colnames(genecounts)
# rownames(metadata)
# create deseq object
star_results <- DESeqDataSetFromMatrix(countData = genecounts, colData = metadata, design = ~tooth_type)
# filter out any genes with fewer than 10 reads total
star_results <- star_results[rowSums(counts(star_results)) >= 10,]
star_results
# set factor level (this determines which direction the comparisions are made -- by default it's by alphabetical order)
star_results$tooth_type <- factor(star_results$tooth_type, levels=c("PD", "PF"))

# Run deseq2
ptm <- proc.time()
se_star <- DESeq(star_results, parallel=TRUE)
proc.time() - ptm 

# how many genes are up or down regulated at adjusted alpha level of 0.05?
res <- results(se_star, alpha=0.05)
# order by p value
res <- res[order(res$padj),]
paste("number of genes with adjusted p value lower than 0.01: ", sum(res$padj < 0.05, na.rm=TRUE))
summary(res)
# [1] "number of genes with adjusted p value lower than 0.01:  5502"

# out of 12443 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2180, 18%
# LFC < 0 (down)     : 3322, 27%
# outliers [1]       : 0, 0%
# low counts [2]     : 750, 6%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# check directionality of results   
resultsNames(se_star)
# [1] "Intercept"           "tooth_type_PF_vs_PD"

# filter out low count genes with LFC shrinkage
# filter out low count genes
resLFC <- lfcShrink(se_star, coef="tooth_type_PF_vs_PD", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
paste("number of genes with adjusted p value lower than 0.05: ", sum(resLFC$padj < 0.05, na.rm=TRUE))
summary(resLFC)
# write results to file
write.table(resLFC, file="deseq_results_denovoADS.txt", quote=F, sep="\t")

# test volcano plot
plotMA(res)
plotMA(resLFC)
     
# transform for visualizations
vld <- vst(se_star, fitType="local")

#Get 50 top varying genes
topVarGenes <- head(order(rowVars(assay(vld)), decreasing=TRUE), 50)

#make a subset of the log transformed counts for just the top 50 varying genes
top50Counts <- assay(vld)[topVarGenes,]
write.csv(top50Counts, file="top50counts.vld.denovoADS.csv", quote=FALSE)
 
#PLOT PCA
#PCA using top 500 varying genes
pdf("pca_ADSdenovo.pdf")
plotPCA(vld, intgroup=c("tooth_type"), ntop=500) + theme_minimal()
dev.off()
plotPCA(vld, intgroup=c("tooth_type"), ntop=500) + theme_minimal()

# heatmap of the top 50 genes (variance between PD and PF)
df <- as.data.frame(colData(vld)[,c("ads_nmol_min_mg_of_protein","tooth_type")])
df$ads_nmol_min_mg_of_protein <- log2(as.numeric(df$ads_nmol_min_mg_of_protein))
colnames(df) <- c("ads", "tooth_type")

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
library(RColorBrewer)
x <- pheatmap(top50Counts, annotation_col = df, color = brewer.pal(9, "Greys"))
save_pheatmap_pdf(x, "heatmap_ADSdenovo.pdf")

# now we want to generate a volcano plot with all genes
library(ggpubr)
options(ggrepel.max.overlaps = Inf)
# get labels for top up and top down regulated genes
upreg <- res[res$log2FoldChange > 2,]
downreg <- res[res$log2FoldChange < -2,]
# rank by p value, get top 10
upregname <- row.names(upreg[rank(upreg$padj) <= 10,])
# bottom 10
downregname <- row.names(downreg[rank(downreg$padj) <= 10,])
regnames <- c(upregname, downregname)

options(ggrepel.max.overlaps = Inf)
p <- ggmaplot(res, main = expression("PF" %->% "PD"),
   fdr = 0.05, fc = 4, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames = as.vector(row.names(res)),
   legend = "top", top = 0,
   font.label = c("bold", 11),
   font.legend = "bold",
   font.main = "bold",
   label.rectangle = T,
   label.select = regnames,
   ggtheme = ggplot2::theme_minimal())
p

# write output results to file to filter through (want to label only those from same transcript, all three or at least two arc genes)
write.table(p$data, "sig_genes_from_volcanoplot_ADSdenovo.txt", sep = "\t", quote=F, row.names=F)
```

Get annotated gene information from each suspected ADS gene

```bash
# upload results from R to hillary for further processing
scp sig_genes_from_volcanoplot_ADSdenovo.txt allie@hillary.clemson.edu:/home/allie/ads_plaque/14-denovo_assembly
# on hillary
# get predicted proteins from deseq results, merge into one file for supplement
awk '{print $1}' sig_genes_from_volcanoplot_ADSdenovo.txt| grep "^N" | while read line; do grep $line -m1 ads_assemblies_prokka_annot/ads_assemblies.gff; done > sig_genes_from_volcanoplot_ADSdenovo_gene_info.txt
# add header to gene info
sed -i 1i'HEADER' sig_genes_from_volcanoplot_ADSdenovo_gene_info.txt
# paste together DEseq results with gene info
paste sig_genes_from_volcanoplot_ADSdenovo.txt sig_genes_from_volcanoplot_ADSdenovo_gene_info.txt > final_ADSdenovo_sig_genes.txt
```

Which are arc genes?

```bash
# pull out arcABC ids
grep -E "eC_number=3.5.3.6|eC_number=2.1.3.3|eC_number=2.7.2.2" final_ADSdenovo_sig_genes.txt > temp
mv temp final_ADSdenovo_sig_genes.txt
# add back in header
sed -i 1i'name\tmean\tlfc\tpadj\tsig\tHEADER' final_ADSdenovo_sig_genes.txt
# so this will be part of the supplementary table -- need next to get taxonomy for each of these transcripts
# pull out arcABC ids
cat final_ADSdenovo_sig_genes.txt | awk '{print $1}' | grep "^N" > sig_genes_arc.ids
# pull transcripts with these genes
# first get corresponding transcript ids
cat sig_genes_arc.ids| while read line; do grep $line ads_assemblies_prokka_annot/ads_assemblies.gff | awk '{print $1}' ; done > sig_genes_transcript.ids
# filter to only include unique transcript ids
cat sig_genes_transcript.ids| sort | uniq > temp
mv temp sig_genes_transcript.ids 
# pull transcripts
seqtk subseq ads_assemblies_prokka_annot/ads_assemblies.fna sig_genes_transcript.ids > arc_transcript.fa
# taxonomic assignment with standard database
kraken2 --db ~/refdb/kraken2_standard_db --threads 20 --use-names --output arc_transcript.kraken.out --confidence 0.01 arc_transcript.fa
```

Are the "truncated" transcripts really partial ADS operons? Let's look at an alignment of all Strep sanguinis transcripts

```bash
# transcript IDs from supplementary tables (v7)
cat sang_transcript.ids| sort | uniq > temp
mv temp sang_transcript.ids
# pull from assembled transcripts
seqtk subseq ads_assemblies_prokka_annot/ads_assemblies.fna sang_transcript.ids > sang_transcript.fa
# align with complementation enabled
mafft --thread 40 \
   --adjustdirectionaccurately \
   sang_transcript.fa > sang_transcript.align.fa

# try with a cleaned up alignment as well
trimal -in sang_transcript.align.fa -out sang_transcript.align.trim.fa -htmlout sang_transcript.align.html -gt 0.5 -resoverlap 0.5 -seqoverlap 50
```

Do the same with the Lepto 212 transcripts

```bash
cat lepto212_transcript.ids| sort | uniq > temp
mv temp lepto212_transcript.ids
seqtk subseq ads_assemblies_prokka_annot/ads_assemblies.fna lepto212_transcript.ids > lepto212_transcript.fa
mafft --thread 40 \
   --adjustdirectionaccurately \
   lepto212_transcript.fa > lepto212_transcript.align.fa
# try with a cleaned up alignment as well
trimal -in lepto212_transcript.align.fa -out lepto212_transcript.align.trim.fa -htmlout lepto212_transcript.align.trim.html -gt 0.5 -resoverlap 0.5 -seqoverlap 50
```

Plan:
What we want to do is demonstrate the increase in diversity that we can gain from the denovo assembly data -- to do this make a sort of insert following format of figure three for two taxa (1) Streptococcus sanguinis and (2) the Leptotrichia. Only using those that have all three genes so that we can make the tree BUT report both the minimum number of unique transcripts (i.e., those that have all three genes) as well as the maximum number of unique transcripts (all transcripts that are assigned to that particular taxon).

First let's do the Leptotrichia (references are those that we pulled from the HOMD database -- by genome coordinates)
There are no Leptotrichia that are down regulated based on the deseq results

```bash
# first get transcript IDs for any leptos that were in deseq results 
cd /home/allie/ads_plaque/14-denovo_assembly/lepto_transcript_trees # (ONLY transcripts with at least arcABC)
cat lepto_transcript.ids| sort | uniq > temp
mv temp lepto_transcript.ids

# from this pull coordinates for each gene, put into supplementary table
# pull operon by coordinates (make mock bedfile)
bedtools getfasta -fi ../ads_assemblies_prokka_annot/ads_assemblies.fna -bed lepto_transcript.coords -fo lepto_transcript.ADS.fa
# get references from HOMD data as well
bedtools getfasta -fi ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.fna -bed homd_lepto.coords -fo homd_lepto.ADS.fa
# concatenate together and align
cat homd_lepto.ADS.fa lepto_transcript.ADS.fa > all_leptoADS.fa
# removing UniqTranscript_11438308 -- very long branch, blast shows is probably a strep of some sort (not significant expression too)
# align
mafft --thread 40 \
   --adjustdirectionaccurately \
   all_leptoADS.fa > all_leptoADS.align.fa
# clean up alignment
trimal -in all_leptoADS.align.fa -out all_leptoADS.align.trim.fa -htmlout all_leptoADS.align.trim.html -gt 0.5 -resoverlap 0.5 -seqoverlap 50
# build a tree
rm *tre
raxmlHPC-PTHREADS-SSE3 -T 40 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s all_leptoADS.align.trim.fa
# clean up tree in fig tree, get order of tips in archaeoptryx to make heatmap of average log fold change and base mean
```

Do the same thing for strep sanguinis and parasanguinis

```bash
cat strep_transcript.ids| sort | uniq > temp
mv temp strep_transcript.ids
# from this pull coordinates for each gene, put into supplementary table
# pull operon by coordinates (make mock bedfile)
bedtools getfasta -fi ../ads_assemblies_prokka_annot/ads_assemblies.fna -bed strep_transcript.coords -fo strep_transcript.ADS.fa
# get references from HOMD data as well
bedtools getfasta -fi ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.fna -bed homd_strep.coords -fo homd_strep.ADS.fa
# concatenate together and align
cat homd_strep.ADS.fa strep_transcript.ADS.fa > all_strepADS.fa
# align
mafft --thread 40 \
   --adjustdirectionaccurately \
   all_strepADS.fa > all_strepADS.align.fa
# clean up alignment
trimal -in all_strepADS.align.fa -out all_strepADS.align.trim.fa -htmlout all_strepADS.align.trim.html -gt 0.5 -resoverlap 0.5 -seqoverlap 50
# build a tree
rm *tre
raxmlHPC-PTHREADS-SSE3 -T 40 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s all_strepADS.align.trim.fa
# clean up tree in fig tree, get order of tips in archaeoptryx to make heatmap of average log fold change and base mean
```





