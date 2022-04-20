gffread ALL_genomes.gff -T -o ALL_genomes.gtf








#PBS -N STAR-map-102PFR_S21
#PBS -l select=1:ncpus=40:mem=750gb
#PBS -l walltime=336:00:00
#PBS -j oe
#PBS -m abe
#PBS -q cugi

# load module
module add star/2.7.5c

# move into scratch
cd /scratch1/amann3/uf_rnaseq/

# run command
STAR --runThreadN 40 \
	--genomeDir kraken_refgenomes/GenomeDir \
	--readFilesIn <(gunzip -c 102PFR_S21_R1_001.fastq.gz) <(gunzip -c 102PFR_S21_R2_001.fastq.gz) \
	--outFileNamePrefix 102PFR_S21 \
	--outSAMtype BAM Unsorted \
	--outReadsUnmapped 102PFR_S21.unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold


ls *P*pbs | while read line; do sed -i 's/STAR-map-//' $line ; done

STAR-map-