### Coverage plots for Strep ADS operon genes

```bash
# first need to pull our reference regions
# comparing S. sanguinis (arcCBA) and S. parasanguinis (arcABC)

bowtie2-build sanguinis_arcCBA.fa sanguinis_arcCBA.db
bowtie2-build parasanguinis_arcABC.fa parasanguinis_arcABC.db
cd ~/ads_plaque/01-read_processing/raw/cutadapt/

# map to sanguinis
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/sanguinis_arcCBA.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/{}.sanguinis.sam 2>~/ads_plaque/13-coverage_plots/{}.sanguinis.out 1>~/ads_plaque/13-coverage_plots/{}.sanguinis.err'
# map to parasanguinis
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/parasanguinis_arcABC.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/{}.parasanguinis.sam 2>~/ads_plaque/13-coverage_plots/{}.parasanguinis.out 1>~/ads_plaque/13-coverage_plots/{}.parasanguinis.err'

# convert to bam and sort by coordinate
ls *sam | sed 's/.sam//' | while read line; do samtools sort $line.sam -o $line.sort.bam; done
# get coverage (quick visualization)
samtools coverage *.sanguinis.sort.bam -A
# these look like pretty even coverage? -- what about read counts though?
# merge bam files by health status
samtools merge all_PD.sanguinis.sort.bam *PD*.sanguinis.sort.bam
samtools merge all_PF.sanguinis.sort.bam *PF*.sanguinis.sort.bam
samtools merge all_PD.parasanguinis.sort.bam *PD*.parasanguinis.sort.bam
samtools merge all_PF.parasanguinis.sort.bam *PF*.parasanguinis.sort.bam
# get depth by position for each
samtools depth all_PD.parasanguinis.sort.bam > all_PD.parasanguinis.depth
samtools depth all_PF.parasanguinis.sort.bam > all_PF.parasanguinis.depth
samtools depth all_PD.sanguinis.sort.bam > all_PD.sanguinis.depth
samtools depth all_PF.sanguinis.sort.bam > all_PF.sanguinis.depth
```

Plot in R

```R
pd.para <- read.table("all_PD.parasanguinis.depth", header=F, sep="\t")

pf.para <- read.table("all_PF.parasanguinis.depth", header=F, sep="\t")

plot(pd.para$V2, pd.para$V3, type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3)
lines(pf.para$V2, pf.para$V3, col="#44CCD0", lwd=3)
abline(v=1229, col="black")
abline(v=2330, col="black")
abline(v=3360, col="black")

pd.san <- read.table("all_PD.sanguinis.depth", header=F, sep="\t")
pf.san <- read.table("all_PF.sanguinis.depth", header=F, sep="\t")
plot(pd.san$V2, log10(pd.san$V3), type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3)
lines(pf.san$V2, log10(pf.san$V3), col="#44CCD0", lwd=3)
abline(v=947, col="black")
abline(v=2130, col="black")
abline(v=3410, col="black")
```

Redo mapping to HOMD species

```bash
bowtie2-build Smitis_arcABC.fa Smitis_arcABC.db
bowtie2-build Ssanguinis_arcABC.fa Ssanguinis_arcABC.db
bowtie2-build Sparasanguinis_arcABC.fa Sparasanguinis_arcABC.db
cd ~/ads_plaque/01-read_processing/raw/cutadapt/
# map to each strep reference
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/Smitis_arcABC.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/{}.smitis.sam 2>~/ads_plaque/13-coverage_plots/{}.smitis.out 1>~/ads_plaque/13-coverage_plots/{}.smitis.err'
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/Ssanguinis_arcABC.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/{}.ssang.sam 2>~/ads_plaque/13-coverage_plots/{}.ssang.out 1>~/ads_plaque/13-coverage_plots/{}.ssang.err'
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/Sparasanguinis_arcABC.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/{}.sparasang.sam 2>~/ads_plaque/13-coverage_plots/{}.sparasang.out 1>~/ads_plaque/13-coverage_plots/{}.sparasang.err'

# convert to bam and sort by coordinate
ls *sam | sed 's/.sam//' | while read line; do samtools sort $line.sam -o $line.sort.bam; done

# merge bam files by health status
samtools merge all_PD.smitis.sort.bam *PD*smitis*.bam
samtools merge all_PF.smitis.sort.bam *PF*smitis*.bam

samtools merge all_PD.sang.sort.bam *PD*sang*.bam
samtools merge all_PF.sang.sort.bam *PF*sang*.bam

samtools merge all_PD.parasang.sort.bam *PD*parasang*.bam
samtools merge all_PF.parasang.sort.bam *PF*parasang*.bam

# get depth by position for each
ls all* | while read line; do samtools depth $line > $line.depth; done
```

Plot in R 

```R
# parasanguinis
pd <- read.table("all_PD.parasang.sort.bam.depth", header=F, sep="\t")

pf <- read.table("all_PF.parasang.sort.bam.depth", header=F, sep="\t")

pdf("all.parasang.pdf")
plot(pd$V2, log10(pd$V3), type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3, ylim=c(1.0,5.0), xlim=c(1,3600))
lines(pf$V2, log10(pf$V3), col="#44CCD0", lwd=3)
abline(v=1229, col="black")
abline(v=1314, col="black")
abline(v=2330, col="black")
abline(v=2416, col="black")
dev.off()

# sanguinis
pd <- read.table("all_PD.sang.sort.bam.depth", header=F, sep="\t")

pf <- read.table("all_PF.sang.sort.bam.depth", header=F, sep="\t")

pdf("all.sang.pdf")
plot(pd$V2, log10(pd$V3), type="l", xlab="base position", ylab="log10 coverage", col="#FA918B", lwd=3, ylim=c(1.0,5.0), xlim=c(1,3600))
lines(pf$V2, log10(pf$V3), col="#44CCD0", lwd=3)
abline(v=1229, col="black")
abline(v=1280, col="black")
abline(v=2296, col="black")
abline(v=2463, col="black")
dev.off()

# mitis
pd <- read.table("all_PD.smitis.sort.bam.depth", header=F, sep="\t")

pf <- read.table("all_PF.smitis.sort.bam.depth", header=F, sep="\t")

pdf("all.smitis.pdf")
plot(pf$V2, log10(pf$V3), type="l", xlab="base position", ylab="log10 coverage", col="#44CCD0", lwd=3, ylim=c(1.0,5.0), xlim=c(1,3600))
lines(pd$V2, log10(pd$V3), col="#FA918B", lwd=3)
abline(v=947, col="black")
abline(v=1260, col="black")
abline(v=2276, col="black")
abline(v=2354, col="black")
dev.off()
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

Again with the two leptotrichia oral taxa

```bash
# pull gff files for both
curl https://www.homd.org/ftp/genomes/PROKKA/current/gff/SEQF2444.gff > SEQF2444.gff
curl https://www.homd.org/ftp/genomes/PROKKA/current/gff/SEQF2748.gff > SEQF2748.gff
# pull arc genes
ls SEQF*gff | while read line; do grep 'eC_number=3.5.3.6' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done 
ls SEQF*gff | while read line; do grep 'eC_number=2.1.3.3' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done 
ls SEQF*gff | while read line; do grep 'eC_number=2.7.2.2' $line | sed 's/;/\t/g' | awk -F"\t" '{print $1 "\t" $4 "\t" $5}' | sed 's/|/\t/'; done 
# clean up in excel (make sure sample ids are the same, genes in right orientation)

# pull relevant sequences
awk -F"\t" '{print $1 "_" $2}' lepto_arc_operon_coordinates.txt > fasta.ids
seqtk subseq ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.fna fasta.ids > full.fa
# now using coordinates, pull operon sequences
bedtools getfasta -fi full.fa -bed bed.txt > lepto_arc_operon.fna
```

Now can map to the full arc operon lepto data

```bash
bowtie2-build lepto_arc_operon.fna lepto_arc_operon.db
cd ~/ads_plaque/01-read_processing/raw/cutadapt/
# map to database
ls *.1.trim.fastq.gz | sed 's/.1.trim.fastq.gz//' | parallel -j 40 --gnu 'bowtie2 -x ~/ads_plaque/13-coverage_plots/homd_gff/lepto/lepto_arc_operon.db -1 {}.1.trim.fastq.gz -2 {}.2.trim.fastq.gz --end-to-end  --qc-filter --no-unal -t -S ~/ads_plaque/13-coverage_plots/homd_gff/lepto/{}.sam 2>~/ads_plaque/13-coverage_plots/homd_gff/lepto/{}.out 1>~/ads_plaque/13-coverage_plots/homd_gff/lepto/{}.err'
# convert to bam and sort by coordinate
ls *sam | sed 's/.sam//' | while read line; do samtools sort $line.sam -o $line.sort.bam; done

# merge bam files by health status
samtools merge all_PD.sort.bam *PD*.bam
samtools merge all_PF.sort.bam *PF*.bam

# get depth for all pooled together
ls all* | while read line; do samtools depth $line > $line.depth; done

# get depth by position for each lepto species
grep "SEQF2444" all_PD.sort.bam.depth > SEQF2444_PD.depth
grep "SEQF2444" all_PF.sort.bam.depth > SEQF2444_PF.depth
grep "SEQF2748" all_PD.sort.bam.depth > SEQF2748_PD.depth
grep "SEQF2748" all_PF.sort.bam.depth > SEQF2748_PF.depth
```

Plot in R 

```R
pd <- read.table("SEQF2444_PD.depth", header=F, sep="\t")
pf <- read.table("SEQF2444_PF.depth", header=F, sep="\t")

pdf("SEQF2444_cov.plot.pdf")
plot(pf$V2, log10(pf$V3), type="l", xlab="base position", ylab="log10 coverage", col="#44CCD0", lwd=3, xlim=c(1,3600), ylim=c(1.0, 5.0))
lines(pd$V2, log10(pd$V3), col="#FA918B", lwd=3)
abline(v=1227, col="black")
abline(v=1359, col="black")
abline(v=2351, col="black")
abline(v=2354, col="black")
abline(v=3289, col="black")
dev.off()

pd <- read.table("SEQF2748_PD.depth", header=F, sep="\t")
pf <- read.table("SEQF2748_PF.depth", header=F, sep="\t")

pdf("SEQF2748_cov.plot.pdf")
plot(pf$V2, log10(pf$V3), type="l", xlab="base position", ylab="log10 coverage", col="#44CCD0", lwd=3, xlim=c(1,3600), ylim=c(1.0, 5.0))
lines(pd$V2, log10(pd$V3), col="#FA918B", lwd=3)
abline(v=935, col="black")
abline(v=937, col="black")
abline(v=1929, col="black")
abline(v=2075, col="black")
abline(v=3292, col="black")
dev.off()
```