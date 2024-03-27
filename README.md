# extract-COI

Collection of small scripts to identify and extract candidate COI regions from tblastn result tables formatted in a specific way. 

## Dependencies
You need the following external software to run this workflow:
- bedtools   (https://github.com/arq5x/bedtools2)
- ncbi-blast (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

You also need these python packages:
- pandas
- biopython

## How to run

First, clone the repository and go inside the cloned directory:
```sh
git clone https://github.com/etkayapar/extract-COI
cd extract-COI
```

Originally the Bombyx mori COI sequence queried against input genome assemblies (it is included in this repo for you to easily test). To change the query COI sequence, replace the `COI_query.fasta` file with the query sequence of your choice.

To run the entire workflow, you need to have a plaintext file with absolute paths to your assemblies (one assembly for each line in the file).

Then run the workflow:

```sh
scripts/main.sh <GENOME_PATHS_FILE> <NUMBER_OF_THREADS>
```

After a successful run you will find your nucleotide and amino-acid COI sequences extracted in the project root directory named as:
```
all_COI.fa
all_COI.faa
```

