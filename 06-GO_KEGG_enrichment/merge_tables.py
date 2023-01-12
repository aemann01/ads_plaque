import pandas as pd

# read in files
read_counts = pd.read_csv("../02-STAR/homd_mapped_strict/featurecounts/read_counts.txt", delimiter="\t")
# clean up empty column
read_counts = read_counts.dropna(axis=1)
unprot = pd.read_csv("locus2uniprot.txt", delimiter="\t")
uni2ref = pd.read_csv("gene_refseq_uniprotkb_collab", delimiter="\t")
gene2acc = pd.read_csv("gene_protein_acc", delimiter="\t")
kegg = pd.read_csv("ALL_genomes_KO.txt", delimiter="\t")
goterms = pd.read_csv("gene2go", delimiter="\t")

# merge files
merged = pd.merge(unprot, read_counts, right_on="Geneid", left_on="Locus_tag")
merge2 = pd.merge(uni2ref, merged, right_on="Uniprot", left_on="UniProtKB_protein_accession")
merge3 = pd.merge(gene2acc, merge2, right_on="#NCBI_protein_accession", left_on="protein_accession")
merge4 = pd.merge(goterms, merge3, left_on="GeneID", right_on="GeneID")
merge5 = pd.merge(kegg, merge4, right_on="Locus_tag", left_on="SeqID")

# drop duplicate entries
merge5 = merge5.drop_duplicates("Locus_tag")

# write merged data to file
with open("read_counts_goKegg.txt", "w") as outfile:
	merge5.to_csv(outfile, sep="\t", index=False)
