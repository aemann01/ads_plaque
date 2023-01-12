neighbor joining trees of strep and leptotrichia arcA B and C

```bash
grep "Arginine deiminase" ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' > arcA.ids

grep "Ornithine carbamoyltransferase" ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' > arcB.ids

grep "Carbamate kinase" ~/ads_plaque/03-STAR/homd_mapped_strict/ALL_genomes.tsv | awk '{print $1}' > arcC.ids

grep "Streptococcus" ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt | awk '{print $1}' > strep.ids

cat strep.ids | while read line; do grep $line arcA.ids ; done > strep_arcA.ids
cat strep.ids | while read line; do grep $line arcB.ids ; done > strep_arcB.ids
cat strep.ids | while read line; do grep $line arcC.ids ; done > strep_arcC.ids

seqtk subseq arcA.fa strep_arcA.ids > strep_arcA.fa
seqtk subseq arcB.fa strep_arcB.ids > strep_arcB.fa
seqtk subseq arcC.fa strep_arcC.ids > strep_arcC.fa

mafft strep_arcA.fa > strep_arcA.align.fa
mafft strep_arcB.fa > strep_arcB.align.fa
mafft strep_arcC.fa > strep_arcC.align.fa

fasttree -nt strep_arcA.align.fa > strep_arcA.tre
fasttree -nt strep_arcB.align.fa > strep_arcB.tre
fasttree -nt strep_arcC.align.fa > strep_arcC.tre

cat strep_arcA.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste strep_arcA.ids temp > strep_arcA.annotations.txt
cat strep_arcB.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste strep_arcB.ids temp > strep_arcB.annotations.txt
cat strep_arcC.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste strep_arcC.ids temp > strep_arcC.annotations.txt

####### Now do the same with leptotrichia
grep "Leptotrichia" ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt | awk '{print $1}' > lepto.ids

cat lepto.ids | while read line; do grep $line arcA.ids ; done > lepto_arcA.ids
cat lepto.ids | while read line; do grep $line arcB.ids ; done > lepto_arcB.ids
cat lepto.ids | while read line; do grep $line arcC.ids ; done > lepto_arcC.ids

seqtk subseq arcA.fa lepto_arcA.ids > lepto_arcA.fa
seqtk subseq arcB.fa lepto_arcB.ids > lepto_arcB.fa
seqtk subseq arcC.fa lepto_arcC.ids > lepto_arcC.fa

mafft lepto_arcA.fa > lepto_arcA.align.fa
mafft lepto_arcB.fa > lepto_arcB.align.fa
mafft lepto_arcC.fa > lepto_arcC.align.fa

fasttree -nt lepto_arcA.align.fa > strep_arcA.tre
fasttree -nt lepto_arcB.align.fa > strep_arcB.tre
fasttree -nt strep_arcC.align.fa > strep_arcC.tre

cat lepto_arcA.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste lepto_arcA.ids temp > lepto_arcA.annotations.txt
cat lepto_arcB.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste lepto_arcB.ids temp > lepto_arcB.annotations.txt
cat lepto_arcC.ids | sed 's/_.*//' | while read line; do grep $line ~/ads_plaque/03-STAR/homd_mapped_strict/SEQID_info.txt ; done > temp
paste lepto_arcC.ids temp > lepto_arcC.annotations.txt
```

Now concatenate sequences together

```bash
seqtk subseq strep_arcA.fa strep_synt_genes.ids > strep_arcA.synt.fa
seqtk subseq strep_arcB.fa strep_synt_genes.ids > strep_arcB.synt.fa
seqtk subseq strep_arcC.fa strep_synt_genes.ids > strep_arcC.synt.fa
sed -i 's/_/ /' strep_arcA.synt.fa strep_arcB.synt.fa strep_arcC.synt.fa
seqkit concat strep_arcA.synt.fa strep_arcB.synt.fa strep_arcC.synt.fa > ads_operon.strep.fa
# remove duplicate entries (more than one arcC sequence)
```

Build concatenated operon tree

```bash
mafft --thread 55 ads_operon.strep.fa > ads_operon.strep.align.fa
raxmlHPC-PTHREADS-SSE3 -T 55 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s ads_operon.strep.align.fa
```

Get expanded streptococcus tree

```bash
# first need to get a list of all the strep in HOMD, merge with the taxa that we've already pulled and then remake tree
# only missing lineages of sanguinis, cristatus, pneumoniae, oralis, and infantis
# wd: /Users/mann/github/ads_plaque/12-arcABC_trees/expanded_strep_tre
grep "Streptococcus" ../../02-STAR/homd_mapped_strict/SEQID_info.txt | grep "sanguinis" | awk '{print $1}' > sang.ids
grep "Streptococcus" ../../02-STAR/homd_mapped_strict/SEQID_info.txt | grep "cristatus" | awk '{print $1}' > crist.ids
grep "Streptococcus" ../../02-STAR/homd_mapped_strict/SEQID_info.txt | grep "pneumoniae" | awk '{print $1}' > pneu.ids
grep "Streptococcus" ../../02-STAR/homd_mapped_strict/SEQID_info.txt | grep "oralis" | awk '{print $1}' > oral.ids
grep "Streptococcus" ../../02-STAR/homd_mapped_strict/SEQID_info.txt | grep "infantis" | awk '{print $1}' > infant.ids
# need to pull arc genes, make sure they are syntenous 
cat crist.ids infant.ids oral.ids pneu.ids sang.ids > temp
cat temp | while read line; do grep $line ../arcA.ids ; done > temp.arcA
cat temp | while read line; do grep $line ../arcB.ids ; done > temp.arcB
cat temp | while read line; do grep $line ../arcC.ids ; done > temp.arcC
# placed all syntenous genes in a file, pull from fasta files
seqtk subseq ../arcA.fa synt_all_strep.ids > all_strep.arcA.fa
seqtk subseq ../arcB.fa synt_all_strep.ids > all_strep.arcB.fa
seqtk subseq ../arcC.fa synt_all_strep.ids > all_strep.arcC.fa
# modify sequence headers
sed -i 's/_/ /' all_strep.arcA.fa all_strep.arcB.fa all_strep.arcC.fa
# concatenate them together
seqkit concat all_strep.arcA.fa all_strep.arcB.fa all_strep.arcC.fa > all_strep.operon.fa
# build tree
mafft --thread 55 all_strep.operon.fa > all_strep.operon.align.fa
raxmlHPC-PTHREADS-SSE3 -T 55 -m GTRCAT -c 25 -e 0.001 -p 31514 -f a -N 100 -x 02938 -n tre -s all_strep.operon.align.fa
```




