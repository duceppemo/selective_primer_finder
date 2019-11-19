# selective_primer_finder
Finds candidate kmers with high variability between an inclusion and an exclusion group to design qPCR assays.

## Usage
```
usage: primer_finder.py [-h] -i /inclusion_folder/ [-e /exclusion_folder/] -o
                        /output_folder/ [-t 48] [-k KMER_SIZE]
                        [-s SKETCH_SIZE]

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
  -t 48, --threads 48   Number of CPU.Default is maximum CPU available(48)
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        kmer size to use for Mash. Default 99.
  -s SKETCH_SIZE, --sketch_size SKETCH_SIZE
  
                        Sketch size to use for Mash. Default 1000.
```
                        
