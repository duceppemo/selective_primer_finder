# selective_primer_finder
Finds candidate kmers with high variability between an inclusion and an exclusion group to design qPCR assays. This tool works as following:
1. Extract shared kmers from inclusion genomes
2. Extract all kmers from exclusion genomes
3. Subtract exclusion kmers to inclusion kmers
4. Dump inclusion-specific kmers and convert to fasta
5. Assemble kmers into contigs
6. Map contigs to an exclusion genome to identify the variants
7. Filter contigs to only the keep the ones that have at least three variants within 21 bp (the average PCR primer length)
8. Make sure that the remaining contigs are present in all inclusion genomes
9. Further filter contigs to only keep the ones with variants conserved in all exclusion genomes (compensate for the fact that the variants are determined using a single exclusion genome)


## Installation
Create and activate a python 3.6 conda virtual environment
```
conda create -n selective_primer_finder python=3.6
conda activate selective_primer_finder
```
Clone github repo
```
git clone https://github.com/duceppemo/selective_primer_finder
```
Install dependencies
```
cd selective_primer_finder
conda install -y --file requirements.txt
```
## Usage
```
usage: primer_finder.py [-h] -i /inclusion_folder/ -e /exclusion_folder/ -o
                        /output_folder/ [-t 48] [-m 459] [-k 99] [-d 1]

Find group specific qPCR assays from fasta files

optional arguments:
  -h, --help            show this help message and exit
  -i /inclusion_folder/, --inclusion /inclusion_folder/
                        Folder that contains the genomes they assay should
                        amplify
  -e /exclusion_folder/, --exclusion /exclusion_folder/
                        Folder that contains the genomes the assay should not
                        amplify
  -o /output_folder/, --output /output_folder/
                        Folder to hold the result files
  -t 48, --threads 48   Number of CPU. Default is maximum CPU available(48)
  -m 459, --memory 459  Memory in GB. Default is 85% of total memory (459)
  -k 99, --kmer_size 99
                        kmer size to use for KMC (1-256). Default 99.
  -d 1, --duplication 1
                        Maximum number of times a kmer can be found in each
                        inclusion genome. Default is 1, meaning repeated
                        regions are discarded.
```
## Output
The output is a fasta file with inclusion-specific bases in lowercase. A cigar string is in the heard to help locate the variants.

## TODO
I would like the output to be a the qPCR assay _per se_. To do so, contigs in the current output file sould be sorted in a way that the sequence from which it would be the easiest to make an inclusion-specific qPCR assay should be first (those are the large insertions). I'm planing to run Primer3 to get the assays and validate them using _in silico_ PCR.
