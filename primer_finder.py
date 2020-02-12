#!/usr/local/env python3

# Dependencies
# conda install bbmap
# conda install pigz
# pip install pysam
# conda install blast
# conda install bowtie2
# conda install samtools
# conda install kmc
# conda install jellyfish
# conda install psutil
# conda install biopython


from argparse import ArgumentParser
from multiprocessing import cpu_count
from primer_finder_methods import Methods
import os
import sys
from glob import glob
from shutil import rmtree
from random import randint
from psutil import virtual_memory
from collections import OrderedDict


class PrimerFinder(object):
    def __init__(self, args):
        # Command line arguments
        self.inclusion = args.inclusion
        self.exclusion = args.exclusion
        self.out_folder = args.output
        self.cpu = args.threads
        self.mem = args.memory
        self.kmer_size = args.kmer_size
        self.ref = args.reference

        # Data
        self.inclusion_fasta_list = list()
        self.exclusion_fasta_list = list()
        self.seq_dict = dict()
        self.kmer_dict = dict([])

        # Run
        self.run()

    def run(self):
        self.run_checks()
        fasta_kmer_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers.fasta'
        assembly_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers_assembly.fasta'
        # self.run_kmc(fasta_kmer_file)
        # self.run_assembly(fasta_kmer_file, assembly_file)
        ordered_dict = self.filter_mapping(assembly_file)
        self.filter_blast(ordered_dict)

        print('Done!')

    def run_checks(self):
        print('Checking a few things...')
        # Check if number of CPU and memory requested are valid
        self.cpu = Methods.check_cpus(self.cpu)
        self.mem = Methods.check_mem(self.mem)

        # Check if folders are not empty
        result = Methods.check_folder_empty(self.inclusion)
        if result == 0:
            raise Exception('Inclusion folder does not contain files with accepted file extensions: {}'.format(
                Methods.accepted_extensions))
        result = Methods.check_folder_empty(self.exclusion)
        if result == 0:
            raise Exception('Exclusion folder does not contain files with accepted file extensions: {}'.format(
                Methods.accepted_extensions))

        # Check if fasta files are valid
        self.inclusion_fasta_list = Methods.get_files(self.inclusion, Methods.accepted_extensions)
        self.exclusion_fasta_list = Methods.get_files(self.exclusion, Methods.accepted_extensions)
        # for my_list in [self.inclusion_fasta_list, self.exclusion_fasta_list]:
        #     for fasta in my_list:
        #         result = Methods.check_fasta(fasta)
        #         if result == 0:
        #             raise Exception('{} is not a valid fasta file'.format(fasta))

        # Check reference exclusion genome
        if self.ref:
            result = Methods.check_ref(self.ref, self.exclusion_fasta_list)
            if result == 0:
                raise Exception('Please provide the absolute path of a file selected from the exclusion folder')

    def run_kmc(self, fasta_kmer_file):
        # Using kmc
        Methods.create_output_folders(self.out_folder)
        Methods.create_output_folders(self.out_folder + '/kmc')

        print('Running KMC on inclusion group...')
        Methods.list_to_file(self.inclusion_fasta_list, self.out_folder + '/inclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/inclusion_list.txt',
                        self.out_folder + '/inclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem,
                        len(self.inclusion_fasta_list), len(self.inclusion_fasta_list))

        print('Running KMC on exclusion group...')
        Methods.list_to_file(self.exclusion_fasta_list, self.out_folder + '/exclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/exclusion_list.txt',
                        self.out_folder + '/exclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem, 1, '1e9')

        print('Finding inclusion-specific kmers...')
        Methods.kmc_subtract(self.out_folder + '/inclusion',
                             self.out_folder + '/exclusion',
                             self.out_folder + '/inclusion-specific',
                             self.cpu)

        # Dump kmers
        dump_file = self.out_folder + '/dump.txt'
        Methods.kmc_transform(self.out_folder + '/inclusion-specific', dump_file, self.cpu)
        # convert dump file to fasta
        Methods.dump_to_fasta(dump_file, fasta_kmer_file)

    def run_assembly(self, fasta_kmer_file, assembly_file):
        # Assemble kmers to reduce the number of sequences to map later, since this is the longest part
        print('Assembling {} kmers...'.format(Methods.count_fasta_entries(fasta_kmer_file)))
        # Methods.assemble_tadpole(fasta_kmer_file, assembly_file, self.mem, self.cpu)
        Methods.assemble_skesa(fasta_kmer_file, assembly_file, self.mem, self.cpu)

        # Check if assembly file empty
        if os.stat(assembly_file).st_size == 0:
            raise Exception('kmers could not be assembled.')

    def filter_mapping(self, assembly_file):
        # Use mapping to get the unique and conserved bases from the CIGAR string
        # Pick one genome randomly from the exclusion group
        # all kmers will be tested against this genome to find differences (insertion, deletion and mismatches)
        if self.ref:
            excl_ref = self.ref
        else:
            excl_ref = self.exclusion_fasta_list[randint(0, len(self.exclusion_fasta_list) - 1)]
        index_prefix = '.'.join(excl_ref.split('.')[:-1])

        print('Indexing exclusion genome (picked one randomly: {})...'.format(os.path.basename(index_prefix)))
        Methods.index_bowtie2(excl_ref, index_prefix, self.cpu)

        print('Mapping {} assembled kmers...'.format(Methods.count_fasta_entries(assembly_file)))
        # parse assembly file to dictionary
        assembly_dict = Methods.fasta_2_dict(assembly_file)
        Methods.parallel_map_bowtie2(self.out_folder + '/bam', index_prefix, assembly_dict, self.cpu)

        print('Filtering assembled kmers by looking at cigar strings...')
        bam_list = glob(self.out_folder + '/bam' + '/*.bam')
        filtered_dict = Methods.filter_kmers(assembly_dict, bam_list)

        # Sort best kmers
        # Want the ones with the most difference first in the fasta file.
        # Just using length of cigar string to do the sorting
        # TODO -> come up with some more advanced way of sorting
        print('Sorting {} filtered kmers...'.format(len(filtered_dict.keys())))
        ordered_dict = Methods.sort_best_kmers(filtered_dict)

        # Write filtered kmers to file
        with open(self.out_folder + '/best_kmers.fasta', 'wt') as best_fh:
            # for ident, info in filtered_dict.items():
            for ident, info in ordered_dict.items():
                best_fh.write('>{} {}\n{}\n'.format(info.ident, info.desc, info.seq))

        rmtree(self.out_folder + '/bam')
        bt2_files = glob(self.exclusion + '/*.bt2')
        for f in bt2_files:
            os.remove(f)

        return ordered_dict

    def filter_blast(self, ordered_dict):
        # TODO -> Find an automated way to pick the best candidate kmers
        #         not too sure what to prioritize? Insertions, deletions or mismatches?
        # TODO -> Automatocally create qPCR assays with Primer3
        #         Have to make sure that the regions with insertions, deletions or mismatches
        #         are in the primers (3' end ideally) or in the probe. I think it's better if the specificity
        #         is comming from the primers than the probe, but have no evidence for that!
        # TODO -> Test the qPCR assays. Maybe with the in silico PCR?
        #         Score the assays and sort the with the best first

        # Further filter kmers using blast
        print('Filtering kmers using blast...')
        # Merge all exclusion genomes into a single file and make a blast database with it
        blast_out = self.out_folder + '/blast'
        Methods.create_output_folders(blast_out)

        # First: check if assembled inclusion kmers are in all inclusion genomes
        # Only keep the shared ones
        # TODO ->  This part I'm working on right now.
        print('\tIndexing blast database of all inclusion genomes...')
        blast_handle_list = list()
        for incl_genome in self.inclusion_fasta_list:
            Methods.makeblastdb(incl_genome)
            blast_handle_incl = Methods.run_blastn(incl_genome,
                                                   self.out_folder + '/best_kmers.fasta',
                                                   self.cpu,
                                                   len(self.inclusion_fasta_list))
            blast_handle_list.append(blast_handle_incl)

        in_all_incl_list = list()
        for blast_handle in blast_handle_list:
            found = Methods.is_positive_hit(blast_handle)
            in_all_incl_list.append(found[0])

            # Check if the in_all_incl_list returns all True
            # if not, drop the assembled kmer

        # Second: check if different from all exclusion genomes
        print('\tIndexing blast database of all exclusion genomes...')
        Methods.cat_files(self.exclusion_fasta_list, blast_out + '/exclusion_merged.fasta')
        Methods.makeblastdb(blast_out + '/exclusion_merged.fasta')

        # Run the blast
        print('\tRunning blast...')
        blast_handle = Methods.run_blastn(blast_out + '/exclusion_merged.fasta',
                                          self.out_folder + '/best_kmers.fasta',
                                          self.cpu,
                                          len(self.exclusion_fasta_list))
        print('\tFiltering blast results...')
        Methods.filter_blast(blast_handle, self.out_folder + '/final_kmers.fasta', ordered_dict)


if __name__ == "__main__":
    max_cpu = cpu_count()
    max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB

    parser = ArgumentParser(description='Find group specific qPCR assays from fasta files')
    parser.add_argument('-i', '--inclusion', metavar='/inclusion_folder/',
                        required=True,
                        help='Folder that contains the genomes they assay should amplify')
    parser.add_argument('-e', '--exclusion', metavar='/exclusion_folder/',
                        required=False,
                        help='Folder that contains the genomes the assay should not amplify')
    parser.add_argument('-o', '--output', metavar='/output_folder/',
                        required=True,
                        help='Folder to hold the result files')
    parser.add_argument('-t', '--threads', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU. Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-m', '--memory', metavar=str(max_cpu),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85% of total memory ({})'.format(max_mem))
    parser.add_argument('-k', '--kmer_size',
                        required=False,
                        default='99',
                        help='kmer size to use for Mash. Default 99.')
    parser.add_argument('-r', '--reference', metavar='/exclusion_folder/exclusion_genome.fasta',
                        required=False,
                        help='Force a specific genome from exclusion group to be used for variant identification. '
                             'Default: pick on randomly. Must use absolute path.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    PrimerFinder(arguments)
