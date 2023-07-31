## Denovo assembly analysis (Streptococcus subset)

First trying to assemble just Strep to handle resource limitations on palmetto. If this works, will scale up to full dataset.

Merge all forward and reverse into two files

```bash
# on hillary -- these are the reads post quality/human/rRNA filtering
cd /home/allie/ads_plaque/01-read_processing/raw/filtered
cat *.2.fastq.gz > all.2.fastq.gz
cat *.1.fastq.gz > all.1.fastq.gz
# once these are complete, move to palmetto scratch1 directory
```

Download standard database for kraken2 (pbs script)

```bash
#PBS -N kraken2_ntdb
#PBS -l select=1:ncpus=10:mem=50gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe

# load kraken
module add kraken2/2.1.2

# move into working directory
cd /home/amann3/uf_rnaseq/14-denovo_assembly

kraken2-build --standard -db kraken2_standard_db
```

Start denovo assembly of full dataset -- time this to see if practical for future analyses

```bash
#PBS -N denovo_all
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -m abe
#PBS -q bigmem

# move to folder with data
cd /scratch1/amann3

# load trinity module
module add trinityrnaseq/2.15.1

# run
   Trinity \
   --seqType fq \
   --max_memory 750G \
   --left all.1.fastq.gz \
   --right all.2.fastq.gz \
   --CPU 40 \
   --output /scratch1/amann3/trinity_denovo_all \
   --full_cleanup
```

```bash
# pull strep read ids
qsub -I -l walltime=3:00:00
cd /scratch1/amann3/
grep "Streptococcus" kraken.out | awk -F"\t" '{print $2}' > strep.ids




wc -l strep.ids
cat strep.ids | sort | uniq | wc -l
mkdir strep_assembly
module add seqtk/1.3-r106
seqtk subseq strep.ids all.1.fastq.gz > strep_assembly/strep.1.fastq.gz
seqtk subseq strep.ids all.2.fastq.gz > strep_assembly/strep.2.fastq.gz



```


Pull any transcript that aligns to Streptococcus and Leptotrichia, assemble with Trinity. 



kraken2/2.1.2
trinotate/4.0.1




### Exercise 1: Preprocessing and Assembly

   ```
   # Make a denovo and assembly directory for Trinity outputs
   mkdir -p /output/denovo/1_assembly/Trinity
   
   # Run Trinity on raw fastq files
   Trinity \
   --seqType fq \
   --max_memory 14G \
   --left /input/denovo/raw_data/MM_1.fastq \
   --right /input/denovo/raw_data/MM_2.fastq \
   --CPU 4 \
   --trimmomatic \
   --output /output/denovo/1_assembly/Trinity/trinity_bat_assembly \
   --quality_trimming_params "ILLUMINACLIP:$TRINITY_HOME/trinity-plugins/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25" \
   --full_cleanup
   ```

### Exercise 2: Trinity Assembly Statistics

   ```
   # Make a statistics directory for TrinityStats output
   mkdir -p /output/denovo/1_assembly/Trinity/statistics
   
   # Run TrinityStats on the 'trinity_assembly.Trinity.fasta' file
   /$TRINITY_HOME/util/TrinityStats.pl \
   /output/denovo/1_assembly/Trinity/trinity_bat_assembly.Trinity.fasta \
   > /output/denovo/1_assembly/Trinity/statistics/TrinityStats.txt
   
   # Exit from the Trinity Singularity image
   exit
   ```

### Exercise 3: TransDecoder

   ```
   singularity shell \
   --bind /data2/workshop/rna_2023/input:/input \
   --bind /data2/workshop/rna_2023/$USER:/output \
   /data/gsl/bsl/workshop/3_rna/singularity/rnaseq_tools_denovo_v2.sif
   ```

   ```
   # Make a TransDecoder directory within 1_assembly
   mkdir -p /output/denovo/1_assembly/TransDecoder/final
   
   # Run first step of TransDecoder with the output directory flag
   TransDecoder.LongOrfs \
   -t /output/denovo/1_assembly/Trinity/trinity_bat_assembly.Trinity.fasta \
   --gene_trans_map /output/denovo/1_assembly/Trinity/trinity_bat_assembly.Trinity.fasta.gene_trans_map \
   --output_dir /output/denovo/1_assembly/TransDecoder
   
   blastp -query /output/denovo/1_assembly/TransDecoder/longest_orfs.pep \
   -db /db/uniprot_sprot.fasta \
   -max_target_seqs 1 \
   -outfmt 6 \
   -evalue 1e-5 \
   -mt_mode 1 \
   -num_threads 4 \
   > /output/denovo/1_assembly/TransDecoder/blastp.outfmt6
   ```
