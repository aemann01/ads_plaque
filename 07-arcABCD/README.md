






## trying to get read counts from original star mapping to get differential abundance of just those genes that are essential for ADS pathway

```bash
cat ads.ids | awk '{print $1}' | while read line; do grep -m 1 -w $line ~/ads_plaque/03-STAR/homd_mapped_strict/featurecounts/read_counts.txt ; done > ads.read_counts.txt
head -n 1 ~/ads_plaque/03-STAR/homd_mapped_strict/featurecounts/read_counts.txt > header
cat header ads.read_counts.txt > temp
mv temp ads.read_counts.txt
```

Do DESeq2 analysis on this (differential_abundance_ADSgenes.ipynb)




## Mapping to ADS pathway genes

Pull genes from HOMD annotated with functions involved in the ADS pathway

```bash
grep "Arginine/ornithine antiporter" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv > arcD.ids
grep "Ornithine carbamoyltransferase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv > arcB.ids
grep "Carbamate kinase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv > arcC.ids
grep "Arginine deiminase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv > arcA.ids
cat arcA.ids arcB.ids arcC.ids arcD.ids > ads.ids
# filter from annotated fasta
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn ads.ids > ads.ffn
```

Use these sequences as database to map filtered reads (this cuts down on time takes to run diamond)

```bash
bowtie2-build ads.ffn ads.db
cd ../01-read_processing/raw/filtered/
ls *1.fastq.gz | sed 's/.1.fastq.gz//' | parallel -j 55 'bowtie2 -x ~/ads_plaque/09-arcABCD/ads.db -1 {}.1.fastq.gz -2 {}.2.fastq.gz --end-to-end --qc-filter --no-unal -S ~/ads_plaque/09-arcABCD/{}.ads.sam 1>~/ads_plaque/09-arcABCD/{}.ads.out 2>~/ads_plaque/09-arcABCD/{}.ads.err'
```

Pull concordant pairs for each sample and merge

```bash
cd ~/ads_plaque/09-arcABCD
ls *sam | sed 's/.sam//' | parallel 'samtools view -hf 0x2 {}.sam > {}.concordant.sam'
# convert to fastq
ls *concordant.sam | sed 's/.sam//' | while read line; do samtools fastq $line.sam -1 $line.1.fq -2 $line.2.fq; done
# merge paired end reads
ls *1.fq | sed 's/.1.fq//' | parallel 'pear -f {}.1.fq -r {}.2.fq -o {} 1>{}.pear.out'
```

Fastq to fasta

```bash
ls *.assembled.fastq | sed 's/.fastq//' | while read line; do seqtk seq -a $line.fastq > $line.fasta; done
```

AnnoTree protein database for diamond alignment

```bash
wget https://software-ab.informatik.uni-tuebingen.de/download/megan-annotree/annotree.fasta.gz
# build database
diamond makedb --in annotree.fasta.gz -d annotree
```

Classify sequences with diamond

```bash
ls *fasta | sed 's/.fasta//' | parallel 'diamond blastx -d annotree -q {}.fasta -o {}.daa --sensitive --evalue 1e-10 --query-cover 0.95 -f 100'
```

Merge together PF and PD diamond results

```bash
cat *PD*daa > PD.daa
cat *PF*daa > PF.daa
```

Now can import these into megan to get taxonomic and functional (EC) differences between groups



ARC tree

seqtk subseq ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.ffn seq.ids > seq.ffn
mafft seq.ffn > seq.align.fa
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s seq.align.fa





