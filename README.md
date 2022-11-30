# selective_primer_finder
Finds candidate kmers with high variability between an inclusion and an exclusion group to design qPCR assays. This tool works as following:
1. Extract shared kmers from inclusion genomes (kmc)
2. Extract all kmers from exclusion genomes (kmc)
3. Subtract exclusion kmers to inclusion kmers (kmc)
4. Dump inclusion-specific kmers and convert to fasta (kmc)
5. Assemble kmers into contigs (skesa or spades)
6. Map contigs to an exclusion genome to identify the variants (minimap2)
7. Filter contigs to only the keep the ones that have at least two variants within 21 bp (the average PCR primer length) (samtools)
8. Make sure that the remaining contigs are present in all inclusion genomes (blast)
9. Further filter contigs to only keep the ones with variants conserved in all exclusion genomes (compensate for the fact that the variants are determined using a single exclusion genome) (blast)

*Note* that `selective_primer_finder` is programmed in such a way that it will only keep perfect matches. In other words, a kmer needs to be present in ALL inclusion genomes, with no mismatches, and absent from ALL exclusion genomes. It is thus extremely sensitive the quality of the genomes used as input and proper assignement of genome to the inclusiosn and exclusion groups. Input genomes must be carefully currated to ensure proper results. A tool such as `genome_comparator` (https://github.com/duceppemo/genome_comparator) can be really helpful for that purpose.

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
  -a skesa, --assembler skesa
                        Choose assembler: "skesa" or "spades".
```

## Input
Here's an example about how to create the inclusion and exclusion folders in a semi-automated way, if you don't already have your inclusion and exclusion genomes separated:
1. You need a sorted list of all the genomes you're planing to use (inclusion and exclusion):
```
# If all genomes are in a single folder:
ls /all_genomes | sort > all.list

# If you need to get just the ".fasta" files recursively from a "master folder":
find -L  /all_genomes -type f -name "*.fasta" | sort > all.list  # "-L" is to grab the files present as symbolic links
```
2. You need a sorted text file with only your inclusion genomes:
```
# Make sure the list is sorted
sort inclusion.list > inclusion.list1  # sort the file
mv inclusion.list1 inclusion.list  # overwrite the original file with the sorted one

# Create a sorted exclusion list
# comm -3 big_file small_file > difference_file
comm -3 all.list inclusion.list | sort > exclusion.list
```
3. Create symbolic links of your inclusion and exclusion genomes in their respective folder
```
# Create folders to be used as the inclusion and exclusion inputs for selective_primer_finder and populate it:
[ -d inclusion/ ] || mkdir -p inclusion/
[ -d exclusion/ ] || mkdir -p exclusion/

# Inclusion
rm inclusion/*
for i in $(cat inclusion.list); do
    find -L /all_genomes -type f -name "*"${i}"*" \
    -exec ln -s {} inclusion/ \;
done

# Exclusion
rm exclusion/*
for i in $(cat exclusion.list); do
    find -L /all_genomes -type f -name "*"${i}"*" \
    -exec ln -s {} exclusion/ \;
done
```

## Output
The output is a fasta file with inclusion-specific bases in lowercase. Positions of the variants are in the header of the sequences. The criteria to keep a sequence are the following:
* Must have a least 2 bases difference withing 21 bp.
* Must be absent from all samples in exclusion group.
* Must be present in all samples in inclusion group.

## TODO
* Add the possibility to keep kmers that are not in 100% of the inclusion genomes. User could set a threshold.
* I would like the output to be a the qPCR assay _per se_. To do so, contigs in the current output file sould be sorted in a way that the sequence from which it would be the easiest to make an inclusion-specific qPCR assay should be first (those are the large insertions). I'm planing to run Primer3 to get the assays and validate them using _in silico_ PCR.
