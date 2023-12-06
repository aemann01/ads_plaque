## Homologous gene clustering for ADS pathway genes

Pull files

```bash
# download the protein entries for all homd genomes
wget https://www.homd.org/ftp/genomes/PROKKA/current/faa/ALL_genomes.faa
# remove line wrap 
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ALL_genomes.faa > temp
mv temp ALL_genomes.faa
```

Pull genes annotated for the active enzymes in the ADS pathway (arcABDC) and annotate predicted enzyme

```bash
grep "Arginine deiminase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcA.ids
grep "Ornithine carbamoyltransferase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcB.ids
grep "Carbamate kinase" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcC.ids
grep "Arginine/ornithine antiporter" ../03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' | sort | uniq > arcD.ids
# filter from protein fasta
seqtk subseq ALL_genomes.faa arcA.ids > arcA.faa
seqtk subseq ALL_genomes.faa arcB.ids > arcB.faa
seqtk subseq ALL_genomes.faa arcC.ids > arcC.faa
seqtk subseq ALL_genomes.faa arcD.ids > arcD.faa
# annotate them with their prediced enzymatic function
sed -i 's/^>/>arcA_/' arcA.faa
sed -i 's/^>/>arcB_/' arcB.faa
sed -i 's/^>/>arcC_/' arcC.faa
sed -i 's/^>/>arcD_/' arcD.faa
# merge together into single database
cat arcA.faa arcB.faa arcC.faa arcD.faa > arcABCD.faa
```

Now can use MCL blastline to cluster genes that were not annotated as one of the ADS pathway genes

```bash
# format fasta headers 
grep -v "^>" ALL_genomes.faa > seqs
grep "^>" ALL_genomes.faa | sed 's/ .*//' | sed 's/>//' | awk '{print $0 "|" $0}' | sed 's/_.*|/|/' | sed 's/^/>/' > headers
paste -d '\n' headers seqs > homd_for_mcl.faa
```

```bash
# make blast database
makeblastdb -in arcABCD.faa -dbtype prot -out arcABCD.blast.db
# now run blast 
cat homd_for_mcl.faa | parallel --block 100k --recstart '>' --pipe blastp -db arcABCD.blast.db -evalue 1e-10 -outfmt 6 -query - > blast.out
```

Add in cluster numbers to mclblast output

```bash
seq 1 122 | sed 's/^/Cluster/' | paste - dump.out.blast.out.I18 > cluster_map.txt
```

Now need to get cluster number for each gene in read counts -- only include those that have cluster number, add as column in read count file, collapse by cluster number then do differential abundance

```bash
awk '{print $1}' ../03-STAR/homd_mapped_strict/featurecounts/read_counts.txt | grep -v "Geneid" > locus.tags
# get cluster ids from locus tags
cat locus.tags | while read line; do grep -w -m 1 $line cluster_map.txt || echo "NA"  ; done > clusters.ids





```


```bash
# run mclblastline on output
mclblastline --mcl-I=1.8 --blast-m9 blast.out 1>mclblastline.out 2>mclblastline.err
```

Create cluster map by adding cluster numbers

```bash
seq 1 122 | sed 's/^/Cluster/' | paste - dump.out.blast.out.I18 > cluster_map.txt
# # clean up names and convert to tab format
# sed -i 's/|SEQF[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9]//g' cluster_map.txt
```

Get list of genes to use as reference to map sequences back onto

```bash
sed 's/\t/\n/g' dump.out.blast.out.I18 | sed 's/^.*|//' | grep "^SEQF" | sort | uniq > cluster.ids
seqtk subseq ../03-STAR/homd_mapped_strict/ALL_genomes.ffn cluster.ids > all_ads_genes.ffn
```







Pull clusters that have assigned arcA gene ids

```bash
grep "arcA" ../03-STAR/homd_mapped/ALL_genomes.tsv | awk '{print $1}' > arcA.ids
cat arcA.ids | while read line; do grep $line dump.out.blast.out.I18 ; done | sort | uniq > temp
# how many clusters?
wc -l temp
# 40 gene clusters
sed 's/\t/\n/g' temp | awk -F"|" '{print $2}' > temp.ids
# filter these from the fasta file
seqtk subseq ~/ads_plaque/03-STAR/homd_mapped/ALL_genomes.ffn temp.ids > temp.fa
# dereplicate for testing
vsearch --derep_fulllength temp.fa --output temp.uniq.fa
# align
mafft temp.uniq.fa > temp.align.fa

```











#### 1. Filter reference database to only include those genes that were actually detected

```R
# load in read counts file
dat <- read.table("../03-STAR/homd_mapped/featurecounts/read_counts.txt", header=T, sep="\t", row.names = 1)
# make sure there aren't any columns with NA
dat <- dat[,colSums(is.na(dat)) == 0]
# filter out any rows with fewer than 100 mapped reads 
filtdat <- dat[rowSums(dat)>100,]
# write to file
write.table(filtdat, "filtered.counts.txt", sep="\t")
```

```bash
# use file to filter ffn 
awk '{print $1}' filtered.counts.txt | sed 's/"//g' > filtered.query
seqtk subseq /home/allie/ads_plaque/03-STAR/homd_mapped/ALL_genomes.ffn  filtered.query > ALL_genomes.filt.ffn
grep ">" ALL_genomes.filt.ffn -c
```

#### 2. Cluster reference database at 50% sequence identity to get gene clusters

```bash
# dereplicate full length sequences
vsearch --derep_fulllength ALL_genomes.filt.ffn --strand both --sizeout --threads 50 --output ALL_genomes.uniq.ffn --minseqlength 1
# sort by size
vsearch --sortbylength ALL_genomes.uniq.ffn --output ALL_genomes.sort.ffn --sizeout --threads 50
# cluster
vsearch --cluster_smallmem ALL_genomes.sort.ffn --id 0.50 --sizeout --threads 61 --minseqlength 100 --clusters cluster --consout cluster.fa


ls clusters* | parallel 'size=$(grep -c "^>" {}); if [[ ${size} == 1 ]]; then rm {}; fi'


# evaluate how clusters did based on read grouping
cat arcA.ids | while read line; do find . -type f | xargs grep $line; done > test_cluster.results




```

Homologous gene clustering with cblast

```bash
cd /home/allie/ads_plaque/07-homologous_gene_clust/cblaster_test
wget https://www.homd.org/ftp/genomes/PROKKA/current/gbk/ALL_genomes.gbk
# fix problems with genbank file (name and length colline in the LOCUS linerna )
grep "LOCUS" ALL_genomes.gbk | grep -v "\.1 " | awk '{print $2}' > broken
grep "LOCUS" ALL_genomes.gbk | grep -v "\.1 " | awk '{print $2}' | sed 's/\.1/\.1 /' > fixed
paste broken fixed | sed "s/^/sed -i \'s\//" | sed 's/\t/\//' | sed "s/$/\/\' ALL_genomes.gbk/" | sed 's/\.1/\\.1/g' > fix_gbk.sh




# make database
cblaster makedb ALL_genomes.gbk homd.db
# get protein fasta of genes of interest



```









#### 3. Assign GO terms to clusters

```bash
# first need to download interproscan databases
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.9-50.0/interproscan-5.9-50.0-64-bit.tar.gz
```





#### 2. Pull clusters of interest

```bash
# get list of genes that were significant annotated as arcA -- what clusters do these belong to?
grep "arcA" deseq_results_sig.genes.txt | awk '{print $1}' > ../05-homologous_gene_clust/arcA.ids
cat arcA.ids | while read line; do grep -l $line cluster*; done
# arcA is cluster0
# get list of significant genes annotated as arcB 
grep "arcB" deseq_results_sig.genes.txt | awk '{print $1}' > ../05-homologous_gene_clust/arcB.ids
cat arcB.ids | while read line; do grep -l $line cluster*; done
# all of these were also in cluster 0
grep "arcC" deseq_results_sig.genes.txt | awk '{print $1}' > ../05-homologous_gene_clust/arcC.ids
cat arcC.ids | while read line; do grep -l $line cluster*; done
# arcD gene ids
grep "arcD" deseq_results_sig.genes.txt | awk '{print $1}' > ../05-homologous_gene_clust/arcD.ids
# no sig gene ids for arcT or arcR
# queA gene ids
grep "queA" deseq_results_sig.genes.txt | awk '{print $1}' > ../05-homologous_gene_clust/queA.ids


```
