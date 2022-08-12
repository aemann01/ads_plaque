## GO and KEGG number differential abundance analysis

We first want to get the GO/KEGG numbers associated with our genes and add them to the gene count table. Use prokka2kegg to get kegg pathway annotations for each Uniprot predicted protein in the HOMD genbank file.

```bash
# download genbank file from HOMD and prokka2kegg files
wget https://www.homd.org/ftp/genomes/PROKKA/current/gbk/ALL_genomes.gbk
wget https://github.com/SilentGene/Bio-py/raw/master/prokka2kegg/idmapping_KO.tab.gz
wget https://github.com/SilentGene/Bio-py/raw/master/prokka2kegg/prokka2kegg.py
chmod +x prokka2kegg.py
# run prokka2kegg on genbank file
./prokka2kegg.py -i ALL_genomes.gbk -d idmapping_KO.tab.gz -o ALL_genomes_KO.txt
```

Now get corresponding uniprot ids that map to locus tags from the genbank file

```bash
python3.8 locus_tag_to_uniprot.py
```

Get refseq IDs for each Uniprot ID, then gene ID from each accession

```bash
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz
# unzip
gzip -d gene_refseq_uniprotkb_collab.gz
gzip -d gene2accession.gz
# only want to include specific columns from gene2accession to reduce memory load
awk -F "\t" '{print $2, "\t", $6}' gene2accession > gene_protein_acc
# remove version annotation
sed -i 's/\..*//' gene_protein_acc
# remove spaces
sed -i 's/ //g' gene_protein_acc
# merge files together 
python3 merge_tables.py 
# clean up header
sed -i 's/#//g' read_counts_goKegg.txt
```

Now can open jupyter notebook to run stat analyses






