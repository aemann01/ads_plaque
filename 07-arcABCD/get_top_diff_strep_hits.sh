awk '$5 >= 0.001' deseq_results_all.txt | awk '{print $1}' > top_hits_padj0.001
