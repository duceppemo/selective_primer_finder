#!/usr/local/env python3

from argparse import ArgumentParser
from multiprocessing import cpu_count
from primer_finder_methods import Methods
from concurrent import futures
import os


class PrimerFinder(object):
    def __init__(self, args):
        # Command line arguments
        self.inclusion = args.inclusion
        self.exclusion = args.exclusion
        self.out_folder = args.output
        self.cpu = args.threads
        self.sketch_size = args.sketch_size
        self.kmer_size = args.kmer_size

        # Data
        self.inclusion_fasta_list = list()
        self.exclusion_fasta_list = list()
        self.inclusion_sig_list = list()
        self.exclusion_sig_list = list()

        self.seq_dict = dict()
        self.kmer_dict = dict([])

        # Run
        self.run()

    def run(self):

        print('Checking a few things...')
        # Check if number of CPU requested id valid
        self.cpu = Methods.check_cpus(self.cpu)

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
        for my_list in [self.inclusion_fasta_list, self.exclusion_fasta_list]:
            for fasta in my_list:
                result = Methods.check_fasta(fasta)
                if result == 0:
                    raise Exception('{} is not a valid fasta file'.format(fasta))

        # print('Computing signatures for inclusion...')
        # Methods.parallel_sourmash_compute(self.inclusion_fasta_list,
        #                                   os.path.join(self.out_folder, 'inclusion/'),
        #                                   self.cpu,
        #                                   self.kmer_size,
        #                                   self.sketch_size)
        # print('Computing signatures for exclusion...')
        # Methods.parallel_sourmash_compute(self.exclusion_fasta_list,
        #                                   os.path.join(self.out_folder, 'exclusion/'),
        #                                   self.cpu,
        #                                   self.kmer_size,
        #                                   self.sketch_size)

        # put signature files into list
        self.inclusion_sig_list = Methods.get_files(os.path.join(self.out_folder, 'inclusion/'), ['.sig'])
        self.exclusion_sig_list = Methods.get_files(os.path.join(self.out_folder, 'exclusion/'), ['.sig'])

        # print('Finding common signatures of inclusion group...')
        # Methods.sourmash_signature_intersect(self.inclusion_sig_list,
        #                                      self.out_folder + '/common_to_inclusion.sig')
        # print('Subtracting signatures from exclusion group...')
        # Methods.sourmash_signature_subtract(self.out_folder + '/common_to_inclusion.sig',
        #                                     self.exclusion_sig_list,
        #                                     self.out_folder + '/specific_to_inclusion.sig')

        # Use the first genome of inclusion list to get the kmer sequences.
        # Any would work because sourmash reports only perfect matches,
        # thus present in all genomes from the inclusion group
        csv_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers.csv'

        # print('Fetching kmer sequences...')
        # Methods.sourmash_signature_to_kmer(csv_file,
        #                                    self.out_folder + '/specific_to_inclusion.sig',
        #                                    self.inclusion_fasta_list[0])

        # Parse csv and convert to fasta
        # check if file empty
        # TODO -> only keep the one that applies for 'no kmer'
        if os.stat(csv_file).st_size == 0:
            raise Exception('No unique sequences found.')
        elif sum(1 for line in open(csv_file, 'r')) < 2:  # header and a least one kmer
            raise Exception('No unique sequences found.')
        else:
            Methods.parse_cvs(csv_file, self.seq_dict)

        fasta_kmer_file = self.out_folder + '/inclusion-specific_' + self.kmer_size + '-mers.fasta'
        # if self.seq_dict:
        #     Methods.dict_to_fasta(self.seq_dict, fasta_kmer_file)

        # Check which part of the kmer is really unique
        print('Finding unique regions of kmers...')
        # count kmer in exclusion group
        # merge all exclusion genomes in a single file
        # assume that they are all gzipped or all not gzipped
        # Methods.cat_files(self.exclusion_fasta_list, self.out_folder + '/tmp.fasta')
        # Methods.jellyfish_count(self.out_folder + '/tmp.fasta',
        #                         self.out_folder + '/jellyfish/exclusion_merged.jf',
        #                         self.cpu)
        # os.remove(self.out_folder + '/tmp.fasta')
        '''
        Do not work properly
        
        Methods.parallel_jellyfish_count(self.exclusion_fasta_list, self.out_folder + '/jellyfish/', self.cpu)
        Methods.jellyfish_merge(Methods.get_files(self.out_folder + '/jellyfish/', ['.jf']),
                                self.out_folder + '/jellyfish/')
        '''

        # check all kmers for each unique sequence
        kmer_dict = dict([])
        for ident, seq in self.seq_dict.items():
            cnt_list = list()
            kmers = Methods.sliding_window(seq, 21, step=1)
            for kmer in kmers:
                cnt = Methods.jellyfish_query(kmer, self.out_folder + '/jellyfish/exclusion_merged.jf')
                cnt_list.append(int(cnt))
            kmer_dict[ident] = cnt_list

        # Create a pattern to return a condensed representation of the coverage (0 or >0)
        # Something like the CIRGAR string in the SAM format
        # 21U51C would mean 21 bases with zero coverage in the exclusion group (Unique to inclusion)
        # and 51 bases with coverage in the exclusion group (Common to exclusion)
        # that would tell that using 21-mers there is 1 bp difference at position 21.
        # it not perfect because if there are more than 1 difference within 21 bp, we can't know how many exactly
        coverage_dict = dict()
        coverage_string = ''
        for ident, cnts in kmer_dict.items():
            # report the counts of consecutive zero and non-zero values
            coverage_list = []
            counter = 0
            zero = 0  # track if the previous values was a zero (1) or not (0)
            for cnt in cnts:
                if cnt == 0 and zero == 0:  # hit a new zero value
                    if counter > 0:
                        # coverage_string = coverage_string + str(counter) + 'C'
                        coverage_list.append((counter, 'C'))
                    counter = 1  # reset counter
                    zero = 1
                elif cnt == 0 and zero == 1:  # hit a consecutive zero value
                    counter += 1
                elif cnt > 0 and zero == 1:  # hit a new non-zero value
                    if counter > 0:
                        # coverage_string = coverage_string + str(counter) + 'U'
                        coverage_list.append((counter, 'U'))
                    counter = 1
                    zero = 0
                elif cnt > 0 and zero == 0:  # hit a consecutive non-zero value
                    counter += 1
            if counter > 0 and zero == 1:
                # coverage_string = coverage_string + str(counter) + 'U'
                coverage_list.append((counter, 'U'))
            elif counter > 0 and zero == 0:
                # coverage_string = coverage_string + str(counter) + 'C'
                coverage_list.append((counter, 'C'))
            coverage_dict[ident] = coverage_list

        # find More interesting kmers
        print('Filtering kmers...')
        with open(self.out_folder + '/best_kmers.fasta', 'w') as f:
            best_kmers_dict = dict()
            for ident, cnt_tuples in coverage_dict.items():
                if len(cnt_tuples) > 3:
                    best_kmers_dict[ident] = cnt_tuples
                    f.write('>{} {}\n{}\n'.format(ident, cnt_tuples, self.seq_dict[ident]))
                    continue
                # or if the "U-values" are > 2
                else:
                    for cnt, cov_type in cnt_tuples:
                        if cov_type == 'U' and cnt > 21:
                            best_kmers_dict[ident] = cnt_tuples
                            f.write('>{} {}\n{}\n'.format(ident, cnt_tuples, self.seq_dict[ident]))

        # output best kmer
        # print('Writing best (most diverse) inclusion-specific kmers...')
        # Methods.dict_to_fasta(best_kmers_dict,
        #                       self.out_folder + '/best_kmers.fasta')

        # Create jelly database of exclusion genomes
        # Methods.find_unique_regions(self.seq_dict)

        # find stretches of consecutive bases
        # Will be used to force a primer or probe location
        # Since Sourmash uses only perfect matches, retained kmers can only have 1 bp difference
        # Thus hard to get selective in a real PCR
        print('Done!')


if __name__ == "__main__":
    max_cpu = cpu_count()

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
                        help='Number of CPU.'
                             'Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-k', '--kmer_size',
                        required=False,
                        default='99',
                        help='kmer size to use for Mash. Default 99.')
    parser.add_argument('-s', '--sketch_size',
                        required=False,
                        default='1000',
                        help='Sketch size to use for Mash. Default 1000.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    PrimerFinder(arguments)
