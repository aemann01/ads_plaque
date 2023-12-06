Revision figures

```R
# install.packages("remotes")
# remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
# BiocManager::install("phyloseq")
# install.packages("phytools")
# install.packages("devtools")
# devtools::install_github(repo = "malucalle/selbal")
library(selbal)
library(dplyr)
library(phyloseq)
library(phytools)
library(phyloseq.extended)
library(ggplot2)

# load data running on hillary
seqtab <- read.table("../05-rpoC_processing/sequence_table.merged.txt", header=T, row.names=1)
tax <- read.table("../05-rpoC_processing/taxonomy_bac.txt", header=F, row.names=1, sep="\t")
tree <- read.tree("../05-rpoC_processing/RAxML_bestTree.ref.tre")
tree.root <- midpoint.root(tree)
map <- read.table("../map.txt", sep="\t", header=T, row.names=1)
# create phyloseq object
ps.dat <- phyloseq(otu_table(seqtab, taxa_are_rows=F), sample_data(map), tax_table(as.matrix(tax)), tree.root)
ps.dat
head(sample_data(ps.dat))
```

```R
# rarefaction curve analysis
p <- ggrare(ps.dat, step = 1000, se = TRUE)
p + theme_minimal()
pdf("rarefaction_plots.pdf")
p + theme_minimal()
dev.off()
```

Identify taxa associated with response variable (citrulline expression)

```R
# collapse data to roughly species level to minimize high sparsity
glom <- tax_glom(ps.dat, taxrank=rank_names(ps.dat)[9])
# remove any taxa with fewer than 100 counts and in at least 20% of samples post merging
glom <- filter_taxa(glom, function(x) sum(x > 100) > (0.20*length(x)), TRUE)
# pull data
dat <- as.data.frame(otu_table(glom))
map <- sample_data(glom)

# merge metadata with asv table so response variable in same order 
dat <- merge(dat, map, by="row.names")
# fix row names
rownames(dat) <- dat$Row.names

# define data and response variable
dif <- dim(dat)[2] - dim(map)[2]
x <- dat[,1:dif]
# make sure only numeric data
x <- select_if(x, is.numeric)
# response variable
dif2 <- dim(dat)[2] - 3
y <- dat[,dif2]
z <- data.frame(tooth_type = as.factor(dat$tooth_type))
# using tooth type as the covariate
bal <- selbal(x = x, y = y, covar = z, zero.rep = "bayes")
# run selbal.cv as well
cv.bal <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 100, covar = z, zero.rep = "bayes", seed = 4597)
cv.bal$global.balance
# plot correlation
pdf("global_plot.pdf")
grid.draw(bal$global.plot)
dev.off()
pdf("global_plot.cv.pdf")
grid.draw(cv.bal$global.plot)
dev.off()
pdf("var_barplot.cv.pdf")
cv.bal$var.barplot
dev.off()
pdf("plot_tab.cv.pdf")
plot.tab(cv.bal$cv.tab)
dev.off()
```

What if we don't control for tooth type?

```R
# using tooth type as the covariate
bal <- selbal(x = x, y = y, zero.rep = "bayes")
# run selbal.cv as well
cv.bal <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 100, zero.rep = "bayes", seed = 4597)
cv.bal$global.balance
# plot correlation
pdf("global_plot.nocovar.pdf")
grid.draw(bal$global.plot)
dev.off()
pdf("global_plot.cv.nocovar.pdf")
grid.draw(cv.bal$global.plot)
dev.off()
pdf("var_barplot.cv.nocovar.pdf")
cv.bal$var.barplot
dev.off()
pdf("plot_tab.cv.nocovar.pdf")
plot.tab(cv.bal$cv.tab)
dev.off()
```

Identify taxa associated with response variable (tooth health)

```R
y <- factor(dat$tooth_type)
# run
cv.bal <- selbal.cv(x = x, y = y, n.fold = 5, n.iter = 10, covar = NULL, zero.rep = "bayes")
```

```R
# plot correlation
pdf("global_plot.tooth_health.pdf")
grid.draw(cv.bal$global.plot)
dev.off()
pdf("var_barplot.tooth_health.pdf")
cv.bal$var.barplot
dev.off()
pdf("plot_tab.tooth_health.pdf")
plot.tab(cv.bal$cv.tab)
dev.off()
```

What is the average number of transcripts that aligned to P. acidifaciens?

```bash
cd /home/allie/ads_plaque/03-STAR/homd_mapped_strict/featurecounts
grep "SEQF1851\|SEQF2583" read_counts.txt > read_counts.pacidifaciens.txt
# checked this in excel
```

Want to pull coordinates for all strep in HOMD and then map to each one

```bash
grep "SEQF" ~/ads_plaque/12-arcABC_trees/ads_operon.strep.fa | awk '{print $1}' | sed 's/>//' > strep.ids
sed 's/^/https:\/\/www.homd.org\/ftp\/genomes\/PROKKA\/current\/gff\//' strep.ids | sed 's/$/.gff/' > query
wget -i query
# pull arc genes
ls SEQF*gff | while read line; do grep 'eC_number=3.5.3.6' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done > arcA.coordinates
ls SEQF*gff | while read line; do grep 'eC_number=2.1.3.3' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done > arcB.coordinates
ls SEQF*gff | while read line; do grep 'eC_number=2.7.2.2' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done > arcC.coordinates
# add headers
sed -i '1i SEQID\tACCESSION\tarcA_begin\tarcA_end' arcA.coordinates
sed -i '1i SEQID\tACCESSION\tarcB_begin\tarcB_end' arcB.coordinates
sed -i '1i SEQID\tACCESSION\tarcC_begin\tarcC_end' arcC.coordinates
# seqf2012 has a second arcC downstream but pretty far away, remove from file manually
# SEQF2012	GL872447.1	1218622	1218852
# merge 
paste arcA.coordinates arcB.coordinates arcC.coordinates > arc_operon.coordinates.txt
# clean up in excel (make sure sample ids are the same, genes in right orientation)
scp allie@hillary.clemson.edu:/home/allie/ads_plaque/13-coverage_plots/homd_gff/arc_operon\* .
# pull relevant sequences
awk -F"\t" '{print $1 "_" $2}' arc_operon.coordinates.txt > fasta.ids
seqtk subseq ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.fna fasta.ids > full.fa
# now using coordinates, pull operon sequences
awk -F"\t" '{print $1 "_" $2 "\t" $10 "\t" $11}' arc_operon.coordinates.txt | grep -v "SEQID" | while read line; do bedtools getfasta -fi ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.fna -bed -; done > arc_operon.fna
```

Now can map to the full arc operon strep data

```bash
bowtie2-build arc_operon.fna arc_operon.db
cd ~/ads_plaque/01-read_processing/raw/cutadapt/
# map to database
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/homd_gff/arc_operon.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/homd_gff/{}.sam 2>~/ads_plaque/13-coverage_plots/homd_gff/{}.out 1>~/ads_plaque/13-coverage_plots/homd_gff/{}.err'

# convert to bam and sort by coordinate
ls *sam | sed 's/.sam//' | while read line; do samtools sort $line.sam -o $line.sort.bam; done

# merge bam files by health status
samtools merge all_PD.sort.bam *PD*.bam
samtools merge all_PF.sort.bam *PF*.bam

# get depth by position for each
ls all* | while read line; do samtools depth $line > $line.depth; done
# split by strain
cat strep.ids | while read line; do grep $line all_PD.sort.bam.depth > $line.all_PD.sort.bam.depth ; done
cat strep.ids | while read line; do grep $line all_PF.sort.bam.depth > $line.all_PF.sort.bam.depth ; done

```

Plot in R 

```R
filelist.PD <- list.files(pattern="*PD.sort.bam.depth")
filelist.PF <- list.files(pattern="*PF.sort.bam.depth")

for(i in 1:length(filelist.PD)){
	pd.dat <- read.table(filelist.PD[i], header=F, sep="\t")
	pf.dat <- read.table(filelist.PF[i], header=F, sep="\t")
	pdf(paste(filelist.PD[i], "_cov.plot.pdf", sep=""))
	plot(pd.dat$V2, log10(pd.dat$V3), type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3, ylim=c(1.0,4.0), xlim=c(1,3360))
lines(pf.dat$V2, log10(pf.dat$V3), col="#44CCD0", lwd=3)
dev.off()
}
# most of these are too low coverage, but gordonii actually looks pretty interesting, map gene coordinates for this species
pd.gord <- read.table("SEQF1066.all_PD.sort.bam.depth", header=F, sep="\t")
pf.gord <- read.table("SEQF1066.all_PF.sort.bam.depth", header=F, sep="\t")
# plot
pdf("seqf1066_gordonii_cov.plot.pdf")
plot(pd.gord$V2, log10(pd.gord$V3), type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3, ylim=c(1.0,5.0), xlim=c(1,3600))
lines(pf.gord$V2, log10(pf.gord$V3), col="#44CCD0", lwd=3)
abline(v=947, col="black")
abline(v=1123, col="black")
abline(v=2139, col="black")
abline(v=2189, col="black")
dev.off()
```

Are there specific taxa that are negatively correlated with S. mutans (based on ASV abundace data) -- looking at CLR transformed data

```R
install.packages("remotes")
remotes::install_github("microbiome/microbiome")
install.packages("tidyverse")
library(microbiome)
library(tidyverse)
# rel abund transformation
rel.abund <- transform(ps.dat, "compositional") %>% tax_glom(taxrank=rank_names(ps.dat)[10]) %>% psmelt()
head(rel.abund)
dim(rel.abund)

# first need to generate a dataframe with relative abundance per species of interest grouped by sample name and tooth type
# NOTE: Peptostreptococcus nodatum not in ASV data so not included in this analysis

temp <- rel.abund %>% group_by(Sample, tooth_type) %>%
        reframe(Streptococcus_mutans = Abundance[V9 == "Streptococcus_mutans"],
        Parvimonas_micra = Abundance[V9 == "Parvimonas_micra"],
        Actinomyces_naeslundii = Abundance[V9 == "Actinomyces_naeslundii"],
        Lactobacillus_ultunensis = Abundance[V9 == "Lactobacillus_ultunensis"],
        Olsenella_uli = Abundance[V9 == "Olsenella_uli"],
        Fusobacterium = Abundance[V9 == "Fusobacterium_sp._oral_taxon_370"],
        Bulleidia_extructa = Abundance[V9 == "Bulleidia_extructa"],
        Treponema_vincentii = Abundance[V8 == "Treponema_vincentii"],
        Streptococcus_constellatus = Abundance[V10 == "Streptococcus_constellatus"]
        )
# remove duplicate rows after grouping to avoid duplicate points
temp <- unique(temp)


# run hmisc to get correlation matrix
install.packages("Hmisc")
library(Hmisc)

# correlation matrix (spearman correlation)
taxcor <- rcorr(as.matrix(select_if(temp, is.numeric)))

# correlation plot for supplement
install.packages("corrplot")
library(corrplot)

png("corrplot_smutans.png")
corrplot(as.matrix(taxcor$r), type="upper", order="hclust")
dev.off()

# correlation between abundance of S mutans vs other taxa
cor(temp$Streptococcus_mutans, temp$Parvimonas_micra, method="pearson")
# [1] 0.2802064
cor(temp$Streptococcus_mutans, temp$Actinomyces_naeslundii, method="pearson")
# [1] 0.08489966
cor(temp$Streptococcus_mutans, temp$Lactobacillus_ultunensis, method="pearson")
# [1] -0.05544475
cor(temp$Streptococcus_mutans, temp$Olsenella_uli, method="pearson")
# [1] -0.008819661
cor(temp$Streptococcus_mutans, temp$Fusobacterium, method="pearson")
# [1] -0.07822344
cor(temp$Streptococcus_mutans, temp$Bulleidia_extructa, method="pearson")
# [1] 0.06879705
cor(temp$Streptococcus_mutans, temp$Treponema_vincentii, method="pearson")
# [1] -0.04900711
cor(temp$Streptococcus_mutans, temp$Streptococcus_constellatus, method="pearson")
# [1] -0.05544475

# significance


# scatter plot for each taxa v mutans
png("test.png")
ggplot(temp, aes(x=Parvimonas_micra, y=Streptococcus_mutans)) + geom_point()
dev.off()



```
