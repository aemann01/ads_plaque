# Reference based mapping with KRAKEN

### 1. Download kraken2 nt database

```bash
kraken2-build --download-taxonomy --db kraken_nt
kraken2-build --download-library nt --db kraken_nt
kraken2-build --build --db kraken_nt --kmer-len 45 --threads 8
```

### 2. Get taxonomic assignments with kraken2

```bash
cd //////
ls *1.derep.fa.gz | sed 's/.1.derep.fa.gz//' | while read line; do kraken2 --db /home/allie/refdb/kraken_nt/ --threads 8 --use-names --gzip-compressed --output /home/allie/uf_rnaseq/kraken/$line.out $line.1.derep.fa.gz; done
```

Pull taxonomic IDs

```bash
cat *out | awk -F"\t" '{print $3}' | awk -F"(" '{print $2}' | sed 's/taxid //' | sed 's/)//' | sort | uniq > taxids
ls *out | parallel 'gzip {}'
```

Download the full assembly report from NCBI

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
```

Query using our taxid list for any genomes with gtf files (not all will have this)

```bash
grep -E '^[0-9]' taxids > temp
mv temp taxids
# only want representative genomes and those with gtf files
grep "representative genome" assembly_summary_genbank.txt > rep_genomes.txt
# note: in the command below you'll get a warning where "cyano" can't be found, still works
cat taxids | while read line; do grep -w -m 1 $line rep_genomes.txt | awk -F"\t" '{print $20}' | sed 's/$/\/*genomic.gtf.gz/' ; done > query.gtf
wget -i query.gtf
```

Get associated genomes from gtf file

```bash
# only keep those records that had a corresponding gtf file




sed 's/*genomic.gtf.gz/*genomic.fna.gz/' query.ftp > query.genomes
wget -i query.genomes
```

Concatenate all together 

```bash
# first clean up
rm *rna* *cds* *.1

```

### 3. Generate reference index

NOTE: this step takes a long time but only has to be run once

```bash
cd 04-STAR
mkdir kraken_map && cd kraken_map
STAR --runMode genomeGenerate \
	--genomeFastaFiles ../02-kraken/refgenomes/all_genomes.fna  \
	--runThreadN 8 \
	--limitGenomeGenerateRAM 66959267424 \
	--sjdbGTFfile ../02-kraken/refgenomes/all_genomes.gtf \
	--genomeChrBinNbits 15
```

### 4. Map to reference genome

```bash
STAR --runThreadN 4 \
	--genomeDir GenomeDir \
	--readFilesIn filtered.fastq \
	--outFileNamePrefix star-results \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold 
```

### 5. Get gene counts

Convert bam to sam

```bash
samtools view -h -o \
	star-resultsAligned.toTranscriptome.out.sam \
	star-resultsAligned.toTranscriptome.out.bam
```



