#!/usr/local/env python3

# Dependencies
# python=3.7
# conda install pigz=2.4
# conda install pysam=0.15.4
# conda install blast=2.9.0
# conda install bowtie2=2.3.5
# conda install samtools=1.9
# conda install kmc=3.1.1
# conda install psutil=5.7.0
# conda install skesa=2.3.0
# conda install minimap2=2.17


from argparse import ArgumentParser
from multiprocessing import cpu_count
from primer_finder_methods import Methods
from glob import glob
from shutil import rmtree
from random import randint
from psutil import virtual_memory
import os


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
        self.dup = args.duplication

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
        self.run_kmc(fasta_kmer_file)
        self.run_assembly(fasta_kmer_file, assembly_file)
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

        # Check if kmer size is between 1 and 256
        if 1 > int(self.kmer_size) > 260:
            raise Exception('Valid kmer sizes are integer from 1 to 256')

    def run_kmc(self, fasta_kmer_file):
        # Using kmc
        Methods.create_output_folders(self.out_folder)
        Methods.create_output_folders(self.out_folder + '/kmc')

        print('\nRunning KMC on inclusion group...')
        Methods.list_to_file(self.inclusion_fasta_list, self.out_folder + '/inclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/inclusion_list.txt',
                        self.out_folder + '/inclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem,
                        len(self.inclusion_fasta_list), len(self.inclusion_fasta_list) * self.dup)

        print('\nRunning KMC on exclusion group...')
        Methods.list_to_file(self.exclusion_fasta_list, self.out_folder + '/exclusion_list.txt')
        Methods.run_kmc(self.out_folder + '/exclusion_list.txt',
                        self.out_folder + '/exclusion',
                        self.out_folder + '/kmc/',
                        self.kmer_size, self.cpu, self.mem, 1, '1e9')

        print('\nFinding inclusion-specific kmers (subtracting exclusion from inclusion kmers)...')
        Methods.kmc_subtract(self.out_folder + '/inclusion',
                             self.out_folder + '/exclusion',
                             self.out_folder + '/inclusion-specific',
                             self.cpu)

        # Dump kmers
        print('\nDumping kmers to file...')
        dump_file = self.out_folder + '/dump.txt'
        Methods.kmc_transform(self.out_folder + '/inclusion-specific', dump_file, self.cpu)

        # Check if dump file is empty
        if os.stat(dump_file).st_size == 0:
            raise Exception('Could not find inclusion-specific kmers!')

        # Count line number in dump file
        dmp_lines = 0
        with open(dump_file, 'r', encoding='utf-8', errors='ignore') as f:
            dmp_lines = sum(bl.count("\n") for bl in Methods.count_lines_in_file(f))

        # Only keep kmers present in all inclusion genomes
        print('\n{} kmers were dumped'.format(dmp_lines))

        # convert dump file to fasta
        print('Converting dump file to fasta...')
        Methods.dump_to_fasta(dump_file, fasta_kmer_file)

        # Clean up
        rmtree(self.out_folder + '/kmc')

    def run_assembly(self, fasta_kmer_file, assembly_file):
        # Assemble kmers to reduce the number of sequences to map in a later step
        print('\nAssembling {} kmers...'.format(Methods.count_fasta_entries(fasta_kmer_file)))
        Methods.assemble_skesa(fasta_kmer_file, assembly_file, self.mem, self.cpu)

        # Check if assembly file empty
        if os.stat(assembly_file).st_size == 0:
            raise Exception('kmers could not be assembled.')

    def filter_mapping(self, assembly_file):
        """
        To know what makes the assembled kmer specific to the inclusion group:
          - is it just one SNP?
          - is the whole assembled kmer unique to the inclusion group?
          - is it a multiple SNPs and indels?
          - How far apart there variants are from each other?
          - Can we get a stretch of at least 3 variants within 21 nucleotides to design selective primers?

        The one drawback of this method is that it only uses one exclution genome to check for the variants.
        Those variants may not be present in all exclusion genomes

        :param assembly_file: string; fasta file of assembled kmers
        :return: Ordered dictionary; contains the assembled kmers which are passing retention criteria for primer design
        """

        # Use mapping to get the unique and conserved bases from the CIGAR string
        # all kmers will be tested against this genome to find differences (insertion, deletion and mismatches)
        if self.ref:
            excl_ref = self.ref
        else:
            excl_ref = self.exclusion_fasta_list[randint(0, len(self.exclusion_fasta_list) - 1)]
        index_prefix = '.'.join(excl_ref.split('.')[:-1])

        if self.ref:
            print('Indexing selected exclusion genome: {}...'.format(os.path.basename(index_prefix)))
        else:
            print('Indexing exclusion genome (picked one randomly: {})...'.format(os.path.basename(index_prefix)))
        Methods.index_bowtie2(excl_ref, index_prefix, self.cpu)

        # split assembly file
        fasta_folder = self.out_folder + '/fasta'
        Methods.split_fasta(assembly_file, fasta_folder)

        print('Mapping {} assembled kmers...'.format(Methods.count_fasta_entries(assembly_file)))
        # parse assembly file to dictionary
        assembly_dict = Methods.fasta_2_dict(assembly_file)

        bam_folder = self.out_folder + '/bam'
        # Methods.parallel_map_bowtie2(bam_folder, index_prefix, fasta_folder, self.cpu)
        # output_folder, ref, fasta_folder, cpu
        Methods.parallel_map_minimap2(bam_folder, excl_ref, fasta_folder, self.cpu)

        # Cleanup fasta files
        rmtree(self.out_folder + '/fasta')

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
        """
        This is an attempt to consider all inclusion and all exclusion genomes to keep only the assembled
        kmers that are the most likely to generate inclusion specific PCR assays.

        :param ordered_dict: Ordered dictionary; pre-filtered assembled kmers
        :return: None
        """
        # Further filter kmers using blast
        print('Filtering kmers using blast...')
        # Merge all exclusion genomes into a single file and make a blast database with it
        blast_out = self.out_folder + '/blast'
        Methods.create_output_folders(blast_out)

        # First: check if assembled inclusion kmers are in all inclusion genomes
        # Only keep the shared ones
        print('\tChecking assembled kmer presence in all inclusion genomes...')
        blast_incl_dict = dict()
        # TODO -> make parallel
        for incl_genome in self.inclusion_fasta_list:
            # print('\t\t{}'.format(os.path.basename(incl_genome)))
            Methods.makeblastdb(incl_genome)
            blast_handle_incl = Methods.run_blastn(incl_genome,
                                                   self.out_folder + '/best_kmers.fasta',
                                                   self.cpu,
                                                   len(self.inclusion_fasta_list))
            if blast_handle_incl:
                # Which genomes harbour which assembled kmers?
                blast_incl_dict[os.path.basename(incl_genome)] = Methods.is_positive_hit(blast_handle_incl)
            else:
                # None of the assembled kmers are found in that specific inclusion genome!!!
                blast_incl_dict[os.path.basename(incl_genome)] = ''
            Methods.clean_blast_index_files(incl_genome)

        # Check if the in_all_incl_list returns all True
        # if not, drop the assembled kmer
        incl_dict = dict()
        for incl, tupple_list in blast_incl_dict.items():
            for t in tupple_list:
                contig = t[0]
                present = t[1]
                try:
                    incl_dict[contig].append(present)
                except KeyError:
                    incl_dict[contig] = [present]

        all_incl_out = self.out_folder + '/all_inclusion_contigs.fasta'
        with open(all_incl_out, 'w') as f:
            for contig, present_list in incl_dict.items():
                # print('{}: {}'.format(contig, present_list))
                if all(present_list):
                    # Write to file
                    # use that file for following steps
                    f.write('>{}\n{}\n'.format(contig, ordered_dict[contig].seq))

        # Second: check if different from all exclusion genomes
        print('\tIndexing blast database of all exclusion genomes...')
        Methods.cat_files(self.exclusion_fasta_list, blast_out + '/exclusion_merged.fasta')
        Methods.makeblastdb(blast_out + '/exclusion_merged.fasta')

        # Run the blast
        print('\tRunning blast...')
        blast_handle = Methods.run_blastn(blast_out + '/exclusion_merged.fasta',
                                          # self.out_folder + '/best_kmers.fasta',
                                          all_incl_out,
                                          self.cpu,
                                          len(self.exclusion_fasta_list))
        Methods.clean_blast_index_files(blast_out + '/exclusion_merged.fasta')
        print('\tFiltering blast results...')
        Methods.filter_blast(blast_handle, self.out_folder + '/final_kmers.fasta', ordered_dict)
        print(
            'Final number of contigs: {}'.format(Methods.count_fasta_entries(self.out_folder + '/final_kmers.fasta')))

        # removed merged exlusion genomes file
        rmtree(blast_out)
        # os.remove(blast_out + '/exclusion_merged.fasta')


    # TODO -> Find an automated way to pick the best candidate kmers
    #         not too sure what to prioritize? Insertions, deletions or mismatches?
    # TODO -> Automatocally create qPCR assays with Primer3
    #         Have to make sure that the regions with insertions, deletions or mismatches
    #         are in the primers (3' end ideally) or in the probe. I think it's better if the specificity
    #         is comming from the primers than the probe, but have no evidence for that!
    # TODO -> Test the qPCR assays. Maybe with the in silico PCR?
    #         Score the assays and sort the with the best first


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
    parser.add_argument('-m', '--memory', metavar=str(max_mem),
                        required=False,
                        type=int, default=max_mem,
                        help='Memory in GB. Default is 85%% of total memory ({})'.format(max_mem))
    parser.add_argument('-k', '--kmer_size',
                        required=False,
                        default='99',
                        help='kmer size to use for Mash. Default 99.')
    parser.add_argument('-d', '--duplication', metavar='',
                        required=False,
                        type=int, default=1,
                        help='Maximum number of times a kmer can be found in each inclusion genome. '
                             'Default is 1 time, meaning repeated regions are discarded.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    PrimerFinder(arguments)
