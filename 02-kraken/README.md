# Reference based mapping with KRAKEN

```bash
mkdir /home/allie/uf_rnaseq/kraken
ls *1.derep.fa.gz | sed 's/.1.derep.fa.gz//' | while read line; do kraken2 --db /home/allie/refdb/kraken_nt/ --threads 8 --use-names --gzip-compressed --output /home/allie/uf_rnaseq/kraken/$line.out $line.1.derep.fa.gz; done
```
Get taxonomic ids from kraken results
```bash
cat *out | awk -F"\t" '{print $3}' | awk -F"(" '{print $2}' | sed 's/taxid //' | sed 's/)//' | sort | uniq > taxids
```

Get the reference genome assembly data for all ncbi genomes

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
```

Query using our taxid list (note: not all will have reference genomes) and pull assembly 

```bash
grep -E '^[0-9]' taxids > temp
mv temp taxids
# only want representative genomes
grep "representative genome" assembly_summary_genbank.txt > rep_genomes.txt
cat taxids | while read line; do grep -w -m 1 $line rep_genomes.txt | awk -F"\t" '{print $20}' | sed 's/$/\/*genomic.fna.gz/' ; done > query.ftp
wget -i query.ftp
```

Get all gtf files for your genomes

```bash
sed 's/*genomic.fna.gz/*genomic.gtf.gz/' query.ftp > query.gtf
wget -i query.gtf
```

Concatenate all together and map using STAR

```bash
# first clean up
rm *rna* *cds* *.1

```





###########
# ALIGNMENT
###########
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
	--readFilesIn filtered.fastq \
	--outFileNamePrefix star-results \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold 

#################
# GET GENE COUNTS
#################
# convert bam to sam formatted file
samtools view -h -o \
	star-resultsAligned.toTranscriptome.out.sam \
	star-resultsAligned.toTranscriptome.out.bam
# what genes were identified?
wget http://www.homd.org/ftp/HOMD_prokka_genomes/tsv/ALL_genomes.tsv
# pull from locus id
grep -v "^@" star-resultsAligned.toTranscriptome.out.sam | awk -F"\t" '{print $3}' | while read line; do grep -w -m 1 $line ALL_genomes.tsv ; done > test.txt
