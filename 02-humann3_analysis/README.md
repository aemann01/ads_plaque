# HUMAnN 3.0 Processing

[Full Documentation](https://github.com/biobakery/humann#humann_join_tables)

### 1. Setup environment

```bash
conda create --name biobakery3 python=3.7
conda activate biobakery3
```

Next add channels needed to run HUMAnN 3.0

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
```

Install HUMAnN

```bash
conda install humann -c biobakery
```

If you have problems installing or there are bowtie2 errors you may need to install the libraries below (on Ubuntu)

```bash
sudo apt-get install libtbb2 
sudo apt-get install libtbb2 
```

Download HUMAnN databases

```bash
mkdir refdb
humann3_databases --download chocophlan full refdb/humann3_db --update-config yes &
humann3_databases --download uniref uniref90_diamond refdb/humann3_db --update-config yes &
humann3_databases --download utility_mapping full refdb/humann3_db --update-config yes 
```

If you get a diamond version error:

```bash
conda install -c bioconda diamond=0.9.36
```

Test installation

```bash
humann_test
humann3 -i demo.fastq.gz -o demo_humann3
```

Finally, run HUMAnN

```bash
ls *fastq.gz | sed 's/.fastq.gz//' | parallel -j 60 'humann3 -i {}.fastq.gz -o ../humann3/{} --remove-temp-output 1>../humann3/{}.out 2>../humann3/{}.err'
```

Join tables together

```bash
cd ..
humann_join_tables -i humann3/ -s -o humann3/merged.tsv  
```