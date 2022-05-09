### *rpoC* amplicon processing

#### 1. Setup

Download raw data from ////

```bash
cd $HOME/ads_plaque/04-rpoC_processing
mkdir raw && cd raw
wget ////
cd ..
```

#### 2. Run quality filtering, read merging, ASV generation, chimeria removal with DADA2

```bash
jupyter-notebook DADA2_processing.ipynb
```

#### 3. Taxonomic assignment with Kraken2

```bash
kraken2 --db ~/refdb/kraken_rpoc/ --threads 15 --use-names --output rep_set.kraken.out --unclassified-out rep_set.unclassified.kraken.out --confidence 0.01 rep_set.fa

awk -F"\t" '{print $3}' rep_set.kraken.out | sed 's/^.*(//' | sed 's/taxid //' | sed 's/)//' > query

cat query | while read line; do grep -w -m 1 ^$line ~/refdb/ncbi_taxonomy/fullnamelineage.dmp | awk -F"|" '{print $3, $2}' | sed 's/\t//g' | sed 's/  / /g' | sed 's/cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages

awk '{print $2}' rep_set.kraken.out > asvids
paste asvids lineages > taxonomy.txt

grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" | awk '{print $1}' > wanted.ids
seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa
grep "Bacteria" taxonomy.txt | grep -v "Bacteria$" > taxonomy_bac.txt
python fix_taxonomy.py taxonomy_bac.txt > temp
mv temp taxonomy_bac.txt
sed -i 's/;/\t/g' taxonomy_bac.txt

mafft --thread 55 rep_set.filt.fa > rep_set.align.fa
raxmlHPC-PTHREADS-SSE3 -T 25 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n ref.tre -s rep_set.align.fa
```

#### 3a. Build database for dada2 taxonomic assignment

```bash
grep ">" rpoc_ref.fa | awk -F"|" '{print $1}' | sed 's/>//' > ids
cat ids | while read line; do esearch -db assembly -query $line < /dev/null | elink -target taxonomy | efetch -format native -mode xml | grep "ScientificName" | awk -F ">|<" 'BEGIN{ORS=", ";}{print $3;}' | awk -F"," 'BEGIN { OFS=";"}{first = $1; $1=""; print $0, first}' | sed 's/; ;/;/' | sed 's/; cellular organisms; //' | sed 's/; /;/g' | sed 's/ /_/g'; done > lineages

```
