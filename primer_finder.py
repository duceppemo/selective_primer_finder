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

        # Using kmc
        Methods.create_output_folders(self.out_folder)
        Methods.create_output_folders(self.out_folder + '/kmc')

        print('Running KMC on inclusion group...')
        Methods.list_to_file(self.inclusion_fasta_list, self.out_folder + '/inclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/inclusion_list.txt',
                        self.out_folder + '/inclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem, len(self.inclusion_fasta_list))

        print('Running KMC on exclusion group...')
        Methods.list_to_file(self.exclusion_fasta_list, self.out_folder + '/exclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/exclusion_list.txt',
                        self.out_folder + '/exclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem, 1)

        print('Finding inclusion-specific kmers...')
        Methods.kmc_subtract(self.out_folder + '/inclusion',
                             self.out_folder + '/exclusion',
                             self.out_folder + '/inclusion-specific',
                             self.cpu)

        # Dump kmers
        dump_file = self.out_folder + '/dump.txt'
        fasta_kmer_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers.fasta'
        Methods.kmc_transform(self.out_folder + '/inclusion-specific', dump_file, self.cpu)
        # convert dump file to fasta
        Methods.dump_to_fasta(dump_file, fasta_kmer_file)

        # Using Jellyfish
        # Methods.create_output_folders(self.out_folder)
        # Methods.create_output_folders(self.out_folder + '/jellyfish')
        # print('Running Jellyfish on inclusion group...')
        # Only keep the kemrs that are in all genomes of inclusion group
        # comparing kmer counts to number of inclusion for that purpose
        # Not perfect since a kmer can be found multiple times in a single genome
        # Methods.jellyfish_count_parallel(self.inclusion_fasta_list, self.out_folder + '/jellyfish/inclusion/',
        #                                  self.cpu, self.kmer_size, 1)
        # jf_list = glob(self.out_folder + '/jellyfish/inclusion' + '/*.jf')
        # Methods.jellyfish_merge(jf_list, self.out_folder + '/jellyfish/inclusion_merged.jf',
        #                         len(self.inclusion_fasta_list))
        # Methods.jellyfish_count(self.out_folder + '/tmp.fasta',
        #                         self.out_folder + '/jellyfish/inclusion_merged.jf',
        #                         self.cpu, self.kmer_size, len(self.inclusion_fasta_list))  # Not perfect
        # os.remove(self.out_folder + '/tmp.fasta')
        #
        # # Any kmer present at least 1 time in the exclusion group will we remove from the inclusion group
        # Methods.cat_files(self.exclusion_fasta_list, self.out_folder + '/tmp.fasta')
        # print('Running Jellyfish on exclusion group...')
        # Methods.jellyfish_count(self.out_folder + '/tmp.fasta',
        #                         self.out_folder + '/jellyfish/exclusion_merged.jf',
        #                         self.cpu, self.kmer_size, 1)
        # os.remove(self.out_folder + '/tmp.fasta')
        #
        # # Dump kmers to dictionary
        # print('Dumping inclusion kmers to dictionary...')
        # incl_dict = Methods.jellyfish_dump(self.out_folder + '/jellyfish/inclusion_merged.jf')
        # print('Dumping exclusion kmers to dictionary...')
        # Methods.jellyfish_dump_to_file(self.out_folder + '/jellyfish/exclusion_merged.jf',
        #                                self.out_folder + '/jellyfish/exclusion_merged.cnt')
        # excl_dict = Methods.jellyfish_dump(self.out_folder + '/jellyfish/exclusion_merged.jf')
        #
        # # Delete jellyfish database files
        # rmtree(self.out_folder + '/jellyfish')
        #
        # # Subtract exclusion kmers to the inclusion kmers
        # print('Subtracting kmers...')
        # self.seq_dict = {str(i): key for i, key in enumerate(incl_dict.keys()) if key not in excl_dict.keys()}
        #
        # # Write kmer to file in fasta format
        # fasta_kmer_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers.fasta'
        # if self.seq_dict:
        #     Methods.dict_to_fasta(self.seq_dict, fasta_kmer_file)
        # else:
        #     raise Exception('No unique kmers found...')

        # Assemble kmers to reduce the number of sequences to map later, since this is the longest part
        print('Assembling {} kmers with tadpole...'.format(Methods.count_fasta_entries(fasta_kmer_file)))
        assembly_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers_assembly.fasta'
        Methods.assemble(fasta_kmer_file, assembly_file, self.mem, self.cpu)

        # Check if assembly file empty
        if os.stat(assembly_file).st_size == 0:
            raise Exception('kmers could not be assembled.')

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

        # TODO -> Find an automated way to pick the best candidate kmers
        #         not too sure what to prioritize? Insertions, deletions or mismatches?
        # TODO -> Automatocally create qPCR assays with Primer3
        #         Have to make sure that the regions with insertions, deletions or mismatches
        #         are in the primers (3' end ideally) or in the probe. I think it's better if the specificity
        #         is comming from the primers than the probe, but have no evidence for that!
        # TODO -> Test the qPCR assays. Maybe with the in silico PCR?
        #         Score the assays and sort the with the best first

        print('Done!')


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
