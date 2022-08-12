## Mapping to ADS pathway genes

Pull genes annotated for the active enzymes in the ADS pathway (arcABDC) and annotate predicted enzyme

```bash
grep "Arginine deiminase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcA.ids
grep "Delta(1)-pyrroline-2-carboxylate reductase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcB.ids
grep "Carbamate kinase 1" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcC.ids
grep "Arginine/ornithine antiporter" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcD.ids
# filter from annotated fasta
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn arcA.ids > arcA.ffn &
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn arcB.ids > arcB.ffn &
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn arcC.ids > arcC.ffn &
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn arcD.ids > arcD.ffn
# annotate with predicted gene
sed -i 's/^>/>arcA_/' arcA.ffn
sed -i 's/^>/>arcB_/' arcB.ffn
sed -i 's/^>/>arcC_/' arcC.ffn
sed -i 's/^>/>arcD_/' arcD.ffn
# merge together into single database
cat arcA.ffn arcB.ffn arcC.ffn arcD.ffn > arcABCD.ffn
```

Use these sequences as database to map filtered reads

```bash
bowtie2-build arcABCD.ffn arcABCD.db
cd ../01-read_processing/raw/filtered/
ls *1.fastq.gz | sed 's/.1.fastq.gz//' | parallel -j 55 'bowtie2 -x ~/ads_plaque/09-arcABCD/arcABCD.db -1 {}.1.fastq.gz -2 {}.2.fastq.gz --end-to-end --qc-filter --no-unal -S ~/ads_plaque/09-arcABCD/{}.ads.sam 1>~/ads_plaque/09-arcABCD/{}.ads.out 2>~/ads_plaque/09-arcABCD/{}.ads.err'
```

Alignment rate (proportion of reads) that mapped for each sample?

```bash
grep "overall alignment rate" *err
# 102PFR_S21.ads.err:0.00% overall alignment rate
# 106PFR_S20.ads.err:0.01% overall alignment rate
# 10PD_S8.ads.err:0.01% overall alignment rate
# 115PD_S36.ads.err:0.01% overall alignment rate
# 12PD_S9.ads.err:0.01% overall alignment rate
# 13PF_S20.ads.err:0.01% overall alignment rate
# 14PF_S15.ads.err:0.01% overall alignment rate
# 15PD_S10.ads.err:0.00% overall alignment rate
# 17PF_S29.ads.err:0.01% overall alignment rate
# 20PD_S6.ads.err:0.01% overall alignment rate
# 21PD_S7.ads.err:0.00% overall alignment rate
# 22PF_S26.ads.err:0.02% overall alignment rate
# 23PF_S25.ads.err:0.01% overall alignment rate
# 24PF_S31.ads.err:0.01% overall alignment rate
# 25PD_S8.ads.err:0.00% overall alignment rate
# 25PER_S1.ads.err:0.01% overall alignment rate
# 26PD_S9.ads.err:0.00% overall alignment rate
# 26PF_S22.ads.err:0.01% overall alignment rate
# 27PF_S23.ads.err:0.01% overall alignment rate
# 30PF_S28.ads.err:0.01% overall alignment rate
# 31PFR_S24.ads.err:0.01% overall alignment rate
# 32PF_S25.ads.err:0.01% overall alignment rate
# 33PF_S34.ads.err:0.01% overall alignment rate
# 35PFR_S19.ads.err:0.01% overall alignment rate
# 36PER_S7.ads.err:0.01% overall alignment rate
# 36PF_S19.ads.err:0.01% overall alignment rate
# 39PD_S10.ads.err:0.00% overall alignment rate
# 40PFR_S35.ads.err:0.01% overall alignment rate
# 42PD_S11.ads.err:0.01% overall alignment rate
# 42PER_S2.ads.err:0.03% overall alignment rate
# 43PD_S12.ads.err:0.01% overall alignment rate
# 44PD_S13.ads.err:0.00% overall alignment rate
# 46PD_S14.ads.err:0.02% overall alignment rate
# 47PD_S5.ads.err:0.02% overall alignment rate
# 47PER_S3.ads.err:0.01% overall alignment rate
# 48PF_S26.ads.err:0.01% overall alignment rate
# 49PDR_S32.ads.err:0.02% overall alignment rate
# 49PD_S13.ads.err:0.01% overall alignment rate
# 49PER_S4.ads.err:0.01% overall alignment rate
# 49PFR_S27.ads.err:0.01% overall alignment rate
# 50PER_S5.ads.err:0.01% overall alignment rate
# 51PD_S6.ads.err:0.03% overall alignment rate
# 52PDR_S35.ads.err:0.02% overall alignment rate
# 52PD_S4.ads.err:0.01% overall alignment rate
# 52PER_S1.ads.err:0.03% overall alignment rate
# 53PDR_S34.ads.err:0.02% overall alignment rate
# 53PD_S15.ads.err:0.01% overall alignment rate
# 53PER_S2.ads.err:0.02% overall alignment rate
# 53PFR_S28.ads.err:0.03% overall alignment rate
# 54PD_S16.ads.err:0.01% overall alignment rate
# 54PER_S3.ads.err:0.01% overall alignment rate
# 55PF_S27.ads.err:0.03% overall alignment rate
# 56PFR_S30.ads.err:0.01% overall alignment rate
# 58PD_S12.ads.err:0.03% overall alignment rate
# 59PD_S17.ads.err:0.01% overall alignment rate
# 60PFR_S33.ads.err:0.01% overall alignment rate
# 61PFR_S23.ads.err:0.01% overall alignment rate
# 64PFR_S18.ads.err:0.01% overall alignment rate
# 67PDR_S11.ads.err:0.01% overall alignment rate
# 68PD_S33.ads.err:0.01% overall alignment rate
# 69PD_S14.ads.err:0.02% overall alignment rate
# 6PF_S16.ads.err:0.01% overall alignment rate
# 70PF_S22.ads.err:0.01% overall alignment rate
# 73PF_S24.ads.err:0.01% overall alignment rate
# 76PDR_S18.ads.err:0.03% overall alignment rate
# 77PF_S31.ads.err:0.01% overall alignment rate
# 79PF_S29.ads.err:0.02% overall alignment rate
# 7PF_S17.ads.err:0.01% overall alignment rate
# 82PFR_S32.ads.err:0.01% overall alignment rate
# 87PFR_S21.ads.err:0.01% overall alignment rate
# 8PF_S30.ads.err:0.02% overall alignment rate
```

Really low across the board -- what about read counts themselves?

```bash
ls *sam | while read line; do grep -v "^@" $line | awk '{print $1}' | sort | uniq | wc -l; done
```

Pull concordant pairs for each sample

```bash
ls *sam | sed 's/.sam//' | parallel 'samtools view -hf 0x2 {}.sam > {}.concordant.sam'
# convert to fastq
ls *concordant.sam | sed 's/.sam//' | while read line; do samtools fastq $line.sam -1 $line.1.fq -2 $line.2.fq; done
```

Merge paired end reads with pear

```bash
ls *1.fq | sed 's/.1.fq//' | parallel 'pear -f {}.1.fq -r {}.2.fq -o {} 1>{}.pear.out'
```

Fastq to fasta

```bash
ls *.assembled.fastq | sed 's/.fastq//' | while read line; do seqtk seq -a $line.fastq > $line.fasta; done
```

Get nr database for diamond, format for diamond

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr
```

Classify sequences with diamond

```bash
ls *fasta | sed 's/.fasta//' | parallel 'diamond blastx -d ~/refdb/nr_diamond/nr -q {}.fasta -o {}.diamond.m8 --sensitive --evalue 1e-10 --query-cover 0.95'
# re running samples that did not complete
cat diamond_redo.ids | sed 's/.fasta//' | while read line; do diamond blastx -d ~/refdb/nr_diamond/nr -q $line.fasta -o $line.diamond.m8 --sensitive --evalue 1e-10 --query-cover 0.95 -p 55; done
```

Merge together PF and PD diamond results

```bash
cat *PF*m8 > PF.diamond.m8
cat *PD*m8 > PD.diamond.m8
```

Now can import these into megan to get taxonomic and functional (EC) differences between groups








