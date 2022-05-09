### Differential abundance analysis

Install and activate the environment

```bash
conda env create -f environment.yml
conda activate 2022-ADS_plaque_diffAbund
```

Add R kernel to Jupyter

```bash
R -e 'IRkernel::installspec()'
```

To turn off the environment

```bash
conda deactivate
```

#### After differential abundance analysis with R

```bash
awk -F"\t" '$7 <= 0.01' deseq_results.txt > deseq_results_sig.txt
awk '{print $1}' deseq_results_sig.txt | while read line; do grep -w -m 1 $line ../02-STAR/homd_mapped/ALL_genomes.tsv ; done > deseq_results_sig.genes.txt
# get gene ref file
awk -F"\t" '{print $1, "\t", $3}' ../02-STAR/homd_mapped/ALL_genomes.tsv > gene_ref



```


```bash
# dereplicate gene tsv file from homd
cat gene_ref | sort | uniq > temp
mv temp gene_ref
```