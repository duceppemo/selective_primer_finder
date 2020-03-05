# selective_primer_finder
Finds candidate kmers with high variability between an inclusion and an exclusion group to design qPCR assays.

## Installation
Create and activate a python 3.6 conda virtual environment
```
conda create -n selectiver_primer_finder python=3.6
conda activate selectiver_primer_finder
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
                        
