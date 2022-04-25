## STAR mapping and post mapping processing

All samples were mapped to the Human Oral Microbiome Database using STAR v 2.7.5c on the Palmetto cluster at Clemson University. PBS scripts for each sample as well as genome indexing instructions can be found in the pbs_scripts folder

### Post mapping processing

#### 1. Re run featurecounts with homd mapped data with transcript_id and CDS as target

```bash
mkdir featurecounts
ls *Aligned.out.bam | sed 's/Aligned.out.bam//' | while read line; do featureCounts -f -p -C -B -a ALL_genomes.gtf -o featurecounts/$line.out -T 60 $line\Aligned.out.bam -t transcript -g transcript_id; done
```

#### 2. Get matrix of gene expression data from all samples

```bash
paste *out | grep -v "^#" | awk '{printf "%s\t", $1}{for (i=7;i<=NF;i+=7) printf "%s\t", $i; printf "\n"}' > read_counts.txt
```
