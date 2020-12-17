#!/usr/local/env python3

from glob import glob
import subprocess
from multiprocessing import cpu_count
from psutil import virtual_memory
import sys
import gzip
import os
import pathlib
from concurrent import futures
from shutil import copyfileobj
import pysam
from collections import OrderedDict
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from io import StringIO
import ahocorasick


class KmerObject(object):
    def __init__(self, ident, desc, seq):
        # Create seq object with its attributes
        self.ident = ident
        self.desc = desc
        self.seq = seq


class Methods(object):
    accepted_extensions = ['.fa', '.fasta', '.fna',
                           '.fa.gz', '.fasta.gz', '.fna.gz']

    @staticmethod
    def get_files(in_folder, accepted_ext):
        """
        Get all files in folder.
        Must supports symbolic links.
        :param in_folder: path of input folder
        :return: list of files in folder
        """
        file_list = list()
        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(in_folder):
            for filename in filenames:
                if filename.endswith(tuple(accepted_ext)):  # accept a tuple or string
                    file_path = os.path.join(root, filename)
                    file_list.append(file_path)
        # return [f for f in glob(in_folder + "**/* _result.txt", recursive=False)]
        return file_list

    @staticmethod
    def check_folder_empty(folder):
        status = 0
        # List content of folder
        dir_content = os.listdir(folder)
        # if folder is not empty and all files have the accepted extensions
        test_file_ext = [x.endswith(tuple(Methods.accepted_extensions)) for x in dir_content]
        if dir_content and any(test_file_ext):
            status = 1
        return status

    @staticmethod
    def check_fasta(fasta):
        status = 0

        with gzip.open(fasta, 'rt') if fasta.endswith('.gz') else open(fasta, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('>'):
                status = 1

        return status

    @staticmethod
    def check_ref(fasta, exclusion_list):
        status = 0

        # Check if used absolute path
        file_path = os.path.dirname(fasta)
        if not file_path:
            status = 0

        # check if in exclusion list
        if fasta in exclusion_list:
            status = 1
        return status

    @staticmethod
    def check_cpus(requested_cpu):
        total_cpu = cpu_count()
        
        if requested_cpu:
            if requested_cpu > total_cpu:
                requested_cpu = total_cpu
                sys.stderr.write("Number of threads was set higher than available CPUs ({})".format(total_cpu))
                sys.stderr.write("Number of threads was set to {}".format(requested_cpu))
        else:
            requested_cpu = total_cpu

        return requested_cpu

    @staticmethod
    def check_mem(requested_mem):
        max_mem = int(virtual_memory().total * 0.85 / 1000000000)  # in GB
        if requested_mem:
            if requested_mem > max_mem:
                requested_mem = max_mem
                sys.stderr.write("Requested memory was set higher than available system memory ({})".format(max_mem))
                sys.stderr.write("Memory was set to {}".format(requested_mem))
        else:
            requested_mem = max_mem

        return requested_mem

    @staticmethod
    def create_output_folders(output):
        """
        Create the output folder and subfolders
        :param output: string; absolute path to output folder
        :return:
        """
        # Create output folder is it does not exist
        pathlib.Path(output).mkdir(parents=True, exist_ok=True)

    @staticmethod
    def dict_to_fasta(seq_dict, fasta_file):
        with open(fasta_file, 'w') as f:
            for k, v in seq_dict.items():
                f.write('>{}\n{}\n'.format(k, v))

    @staticmethod
    def cat_files(fasta_list, output_file):
        with open(output_file, 'w') as f:
            for fasta in fasta_list:
                copyfileobj(open(fasta, 'r'), f)

    @staticmethod
    def assemble_skesa(fasta_file, output_file, mem, cpu):
        cmd = ['skesa', '--cores', str(cpu), '--mem', str(mem), '--fasta', fasta_file, '--contigs_out', output_file]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def assemble_spades(fasta_file, output_dir, mem, cpu):
        cmd = ['spades.py', '--s', '1', fasta_file, '--isolate', '--threads', str(cpu), '-o', output_dir]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def fasta_2_dict(fasta_file):
        my_dict = dict()
        with open(fasta_file, 'r') as f:
            header = None
            seq = list()
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>') and not seq:  # first entry
                    header = ''.join(line[1:]).split(',')[0]  # tadpole specific
                if line.startswith('>') and seq:  # all other enties, except last one
                    my_dict[header] = ''.join(seq)
                    header = ''.join(line[1:]).split(',')[0]  # tadpole specific
                    seq = list()
                if not line.startswith('>'):
                    seq.append(line)
            my_dict[header] = ''.join(seq)  # For the last entry
        return my_dict

    @staticmethod
    def count_fasta_entries(fasta_file):
        cnt = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    cnt += 1
        return cnt

    @staticmethod
    def count_lines_in_file(my_file, size=65536):
        while True:
            b = my_file.read(size)
            if not b:
                break
            yield b

    @staticmethod
    def index_bowtie2(ref, prefix, cpu):
        cmd = ['bowtie2-build', '--threads', str(cpu), ref, prefix]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def map_bowtie2(ref_index, query, cpu, output_bam):
        """
        Can run into an error with bowtie2 if the contigs are too long.
        TODO -> write sequences to fasta files instead of passing them as strings

        :param ref_index:
        :param query:
        :param cpu:
        :param output_bam:
        :return:
        """

        bowtie2_align_cmd = ['bowtie2',
                             '--mm',
                             '--very-sensitive',
                             '--xeq ',
                             '-x', ref_index,
                             '-f', query,
                             '--threads', str(cpu)]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-b', '-h', '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        # samtools_index_cmd = ['samtools', 'index',
        #                       output_bam]
        p1 = subprocess.Popen(bowtie2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
        p2.stdout.close()
        p3.communicate()

    @staticmethod
    def parallel_map_bowtie2(output_folder, ref, fasta_folder, cpu):
        Methods.create_output_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((ref, fasta_file, 1, output_folder + '/' + '.'.join(os.path.basename(fasta_file).split('.')[0]) + '.bam')
                    for fasta_file in glob(fasta_folder + '/*.fasta'))
            # ref_index, query, cpu, output_bam
            for results in executor.map(lambda x: Methods.map_bowtie2(*x), args):
                pass

    @staticmethod
    def map_minimap2(ref, query, cpu, output_bam):
        minimap2_align_cmd = ['minimap2',
                              '-t', str(cpu),
                              '-a', '--eqx',
                              ref, query]
        samtools_view_cmd = ['samtools', 'view',
                             '-@', str(cpu),
                             '-b', '-h', '-']
        samtools_sort_cmd = ['samtools', 'sort',
                             '-@', str(cpu),
                             '-o', output_bam,
                             '-']
        p1 = subprocess.Popen(minimap2_align_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(samtools_sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
        p2.stdout.close()
        p3.communicate()

    @staticmethod
    def parallel_map_minimap2(output_folder, ref, fasta_folder, cpu):
        Methods.create_output_folders(output_folder)
        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((ref, fasta_file, 1, output_folder + '/' + '.'.join(os.path.basename(fasta_file).split('.')[:-1]) + '.bam')
                    for fasta_file in glob(fasta_folder + '/*.fasta'))
            # ref, query, cpu, output_bam
            for results in executor.map(lambda x: Methods.map_minimap2(*x), args):
                pass

    @staticmethod
    def split_fasta(input_fasta, output_folder):
        Methods.create_output_folders(output_folder)
        with open(input_fasta, 'r') as ifh:
            header = ''
            seq = list()
            for line in ifh:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>') and not seq:  # first entry
                    header = line.split()[0].replace('>', '')
                elif line.startswith('>') and seq:  # just getting to a new entry
                    outfile = output_folder + '/' + header + '.fasta'
                    with open(outfile, 'w') as ofh:
                        ofh.write('>{}\n{}\n'.format(header, ''.join(seq)))
                    header = line.split()[0].replace('>', '')
                    seq = list()
                else:
                    seq.append(line)
            else:  # last entry
                outfile = output_folder + '/' + header + '.fasta'
                with open(outfile, 'w') as ofh:
                    ofh.write('>{}\n{}\n'.format(header, ''.join(seq)))

    @staticmethod
    def lower_diff(tuple_list, seq):
        index = 0
        new_seq = ''
        for match_type, length in tuple_list:
            # 1 -> insertion
            # 2 -> deletion
            # 7 -> match
            # 8 -> mismatch
            if match_type == 8 or match_type == 1 or match_type == 2:
                new_seq += ''.join(seq[index:index+length]).lower()
            else:
                new_seq += ''.join(seq[index:index+length])
            index += length
        return new_seq

    @staticmethod
    def tuple_to_cirgar(tuple_list):
        convert = {
            '1': 'I',  # insertion
            '2': 'D',  # deletion
            '3': 'N',  # skipped
            '4': 'S',  # soft clipped
            '5': 'H',  # hard clipped
            '6': 'P',  # padding
            '7': '=',  # match
            '8': 'X'   # mismatch
        }
        cigar = list()
        for match_type, length in tuple_list:
            cigar.extend([str(length), convert[str(match_type)]])
        return str(''.join(cigar))

    @staticmethod
    def filter_kmers(assembly_dict, bam_list):
        filtered_dict = dict()
        for bam in bam_list:
            ident = '.'.join(os.path.basename(bam).split('.')[:-1])
            sam_file = pysam.AlignmentFile(bam, 'rb')

            for s in sam_file:
                cigar_list = s.cigar
                """
                M   0   alignment match (can be a sequence match or mismatch)
                I   1   insertion to the reference
                D   2   deletion from the reference
                N   3   skipped region from the reference
                S   4   soft clipping (clipped sequences present in SEQ)
                H   5   hard clipping (clipped sequences NOT present in SEQ)
                P   6   padding (silent deletion from padded reference)
                =   7   sequence match
                X   8   sequence mismatch
                """
                if not cigar_list:  # All kmer is unique!
                    filtered_dict[ident] = KmerObject(ident, str(s.query_length) + 'I',
                                                      assembly_dict[ident].lower())
                elif len([x for x, y in cigar_list if x == 8]) > 1:
                    # Assuming that a tuple starting with 8 is followed by a tuple starting with 7 and them a 8
                    for i, c in enumerate(cigar_list):
                        if c[0] == 8:
                            # of two mismatched regions are less than 21 bp apart
                            if i < len(cigar_list) - 3\
                                    and cigar_list[i + 2][0] in (8, 1, 2)\
                                    and cigar_list[i + 1][1] < 21:
                                filtered_dict[ident] = KmerObject(ident, Methods.tuple_to_cirgar(cigar_list),
                                                                  Methods.lower_diff(cigar_list, assembly_dict[ident]))
                                break
                # elif len(cigar_list) == 1 and cigar_list[0][0] in (1, 2, 8):  # All kmer is unique!
                #     filtered_dict[ident] = KmerObject(ident, Methods.tuple_to_cirgar(cigar_list),
                #                                       Methods.lower_diff(cigar_list, assembly_dict[ident]))
                else:  # Mixture of matches and mismatches
                    for c in cigar_list:
                        if c[0] == 8 and c[1] > 1:  # Two or more consecutive mismatches
                            filtered_dict[ident] = KmerObject(ident, Methods.tuple_to_cirgar(cigar_list),
                                                              Methods.lower_diff(cigar_list, assembly_dict[ident]))
                            break  # no need to check others if present. This is already good
                        elif c[0] == 1 and c[1] > 1:  # insertion > 1bp
                            filtered_dict[ident] = KmerObject(ident, Methods.tuple_to_cirgar(cigar_list),
                                                              Methods.lower_diff(cigar_list, assembly_dict[ident]))
                            break  # no need to check others if present. This is already good
                        elif c[0] == 2 and c[1] > 1:  # deletion > 1bp
                            filtered_dict[ident] = KmerObject(ident, Methods.tuple_to_cirgar(cigar_list),
                                                              Methods.lower_diff(cigar_list, assembly_dict[ident]))
                            break  # no need to check others if present. This is already good
        return filtered_dict

    @staticmethod
    def sort_best_kmers(filtered_dict):
        # Just using length of cigar string to do the sorting
        ordered_dict = OrderedDict(sorted(filtered_dict.items(), key=lambda x: len(x[1].desc), reverse=True))
        return ordered_dict

    @staticmethod
    def create_primer3_setting_file(setting_file, seq_id, seq):
        # Potential interesting tags

        # SEQUENCE_INCLUDED_REGION
        # SEQUENCE_EXCLUDED_REGION

        field_dict = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_TARGET': '37,21',
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': '1',
            'PRIMER_PICK_INTERNAL_OLIGO': '1',
            'PRIMER_PICK_RIGHT_PRIMER': '1',
            'PRIMER_OPT_SIZE': '18',
            'PRIMER_MIN_SIZE': '15',
            'PRIMER_MAX_SIZE': '21',
            'PRIMER_MAX_NS_ACCEPTED': '1',
            'PRIMER_PRODUCT_SIZE_RANGE': '75-100',
            'P3_FILE_FLAG': '1',
            'SEQUENCE_INTERNAL_EXCLUDED_REGION': '37,21',
            'PRIMER_EXPLAIN_FLAG': '1'
        }

        with open(setting_file, 'w') as f:
            for k, v in field_dict.items():
                f.write('{}={}\n'.format(k, v))
            f.write('=\n')  # terminates the record

    @staticmethod
    def run_primer3(setting_file, out_file):
        """
        http://primer3.org/manual.html#invokingPrimer3
        :param setting_file:
        :param out_file:
        :return:
        """
        cmd = ['primer3_core', '--default_version=2', '-format_output',
               '--p3_settings_file=' + setting_file,
               '--output=file_path', out_file, ]
        subprocess.run(cmd)

    @staticmethod
    def list_to_file(my_list, output_file):
        with open(output_file, 'wt') as f:
            for l in my_list:
                f.write('{}\n'.format(l))

    @staticmethod
    def run_kmc(list_file, output_file, work_dir, kmer_size, cpu, mem, min_cnt, max_cnt):
        cmd = ['kmc', '-k{}'.format(kmer_size),
               '-t{}'.format(cpu),
               '-m{}'.format(mem),
               '-fm',
               '-ci{}'.format(min_cnt),
               '-cx{}'.format(max_cnt),
               '@{}'.format(list_file),
               output_file, work_dir]
        subprocess.run(cmd)  # ), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def kmc_subtract(db1, db2, out_db, cpu):
        cmd = ['kmc_tools', '-t{}'.format(cpu), 'simple',
               db1, db2,
               'kmers_subtract',
               out_db]
        subprocess.run(cmd)  # , stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def kmc_transform(db, output_fasta, cpu):
        cmd = ['kmc_tools', '-t{}'.format(cpu), 'transform',
               db,
               'dump',
               output_fasta]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def dump_to_fasta(dmp, fasta):
        cnt = 0
        with open(fasta, 'w') as o:
            with open(dmp, 'r') as i:
                for line in i:
                    line = line.rstrip()
                    if line == '':
                        continue
                    o.write('>{}\n{}\n'.format('kmer_' + str(cnt), line.split()[0]))
                    cnt += 1

    @staticmethod
    def makeblastdb(fasta):
        cmd = ['makeblastdb', '-in', fasta, '-dbtype', 'nucl']
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    @staticmethod
    def clean_blast_index_files(fasta):
        remove_list = list()
        remove_list += glob(fasta + '*.nhr')
        remove_list += glob(fasta + '*.nin')
        remove_list += glob(fasta + '*.nsq')
        remove_list += glob(fasta + '*.nal')

        for f in remove_list:
            os.remove(f)

    @staticmethod
    def run_blastn(ref_db, query, cpu, max_targets):
        """
        Perform blastn using biopython
        :param ref_db: A fasta file for which "makeblastdb' was already run
        :param query: Protein fasta file
        :param cpu: number of threads
        :param max_targets: maximum targets returned by blast
        :return: blast handle
        """
        # if max_targets > 20:
        #     max_targets = 20  # limit number of sequences to compare

        blastn = NcbiblastnCommandline(db=ref_db, query=query, evalue='1e-10',
                                       outfmt=5, max_target_seqs=max_targets,
                                       num_threads=cpu)
        (stdout, stderr) = blastn()
        if stderr and not 'Warning' in stderr:
            raise Exception('There was a problem with the blast:\n{}'.format(stderr))
        # blast_handle = None
        # if stdout.find('Hsp') != -1:
        #     # Convert stdout (string; blastp output in xml format) to IO object
        #     blast_handle = StringIO(stdout)
        blast_handle = StringIO(stdout)
        return blast_handle

    @staticmethod
    def filter_by_blast(blast_handle, output_file, ordered_dict):
        from collections import defaultdict
        records_dict = SearchIO.to_dict(SearchIO.parse(blast_handle, 'blast-xml'))

        with open(output_file, 'w') as f:
            for seq_id, qresult in records_dict.items():
                # query_len = qresult.seq_len
                similarity_dict = defaultdict(list)
                # if inclusion contig has not hit in the exclusion genomes
                if not qresult.hsps:
                    seq = ordered_dict[seq_id].seq
                    desc = ordered_dict[seq_id].desc
                    f.write('>{} {}\n{}\n'.format(seq_id, desc, seq))  # write to to final file
                    continue
                # If hits were found, save the alignments in a dictionary
                # There's one hsp per exclusion genome
                for h in qresult.hsps:
                    # TODO -> Add check for alignment length?
                    similarity_dict[seq_id].append(h.aln_annotation['similarity'])
                    if h.evalue > 0.0000000001:
                        print('>{} {}\n{}\n'.format(seq_id, desc, seq))

                # Have a deeper look at the alignments
                index_list = list()
                common_variants = list()
                for seq_id, sim_string_list in similarity_dict.items():
                    # Grab the position of all the misalignments between the inclusion contig and the exclusion genome
                    for s in sim_string_list:
                        # index of mismatches
                        idx = [i for i, char in enumerate(s) if char != '|']
                        index_list.append(idx)
                    # Only keep the mismatch position common for all exclusion genomes
                    for i in index_list[0]:
                        if all([i in sublist for sublist in index_list]):
                            common_variants.append(i)
                    # A contig needs more than 1 variant nucleotide to be kept
                    if len(common_variants) > 1:
                        # Check the distance between the variants
                        for i, p in enumerate(common_variants):
                            # If there are more than 2 variants within 21 nucleotides, the inclusion contig is kept
                            if i < len(common_variants) - 3\
                                    and common_variants[i + 1] - common_variants[i] < 21:
                                seq = ordered_dict[seq_id].seq
                                # Variants positions are added to the fasta header
                                # Variants are converted to lower in sequence. Non variable nucleotides are upper case
                                f.write('>{} {}\n{}\n'.format(seq_id, common_variants,
                                                              Methods.lower_indexes(seq.upper(), common_variants)))
                                break

    @staticmethod
    def is_positive_hit(blast_handle):
        found_list = list()
        records_dict = SearchIO.to_dict(SearchIO.parse(blast_handle, 'blast-xml'))
        for seq_id, qresult in records_dict.items():
            found = False
            if qresult.hsps:
                for qr in qresult.hsps:
                    if qr.evalue < 0.0000000001:  # 1e-10
                        found = True
                        break
                found_list.append((seq_id, found))
            else:
                found_list.append((seq_id, found))
                # Maybe do some additional checks about the quality and/or length of the match
        return found_list

    @staticmethod
    def lower_indexes(string, index_list):
        out_string = ''
        for i, char in enumerate(string):
            if i in index_list:
                out_string += char.lower()
            else:
                out_string += char
        return out_string

    # Needle and hay stack
    @staticmethod
    def screen_kmc_kmers(f):
        pass

    @staticmethod
    def reverse_complement(seq):
        """
        Generates the reverse complement of a given DNA sequence.
        :param seq:
        :return:
        """
        complement_dict = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A'}

        return ''.join([complement_dict[base] for base in seq[::-1]])

    @staticmethod
    def create_haystack(fasta_file):
        """

        :param fasta_list:
        :return:
        """
        seq = []
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    line = line.rstrip()
                    seq.append(line)
        return ''.join(seq)

    @staticmethod
    def create_needles(kmer_file):
        aho = ahocorasick.Automaton()
        i = 0
        with open(kmer_file, 'r') as f:
            for line in f:
                kmer = line.split()[0]
                kmer_rc = Methods.reverse_complement(kmer)  # A.add_word(key, (idx, key))
                aho.add_word(kmer, (i, kmer))
                i += 1
                aho.add_word(kmer_rc, (i, kmer_rc))
                i += 1
        aho.make_automaton()
        return aho

    @staticmethod
    def find_needles(kmer_file, fasta_list, fasta_kmer_file):
        # To track if kmers are present in all inclusion genomes
        kmer_dict = dict()
        with open(kmer_file, 'r') as kf:
            for line in kf:
                line = line.rstrip()
                (seq, cnt) = line.split()
                kmer_dict[seq] = 0

        aho = Methods.create_needles(kmer_file)
        for fasta in fasta_list:
            haystack = Methods.create_haystack(fasta)
            for end_index, (insert_order, original_value) in aho.iter(haystack):
                try:
                    kmer_dict[original_value] += 1
                except KeyError:
                    try:
                        kmer_dict[Methods.reverse_complement(original_value)] += 1
                    except KeyError:
                        pass

        with open(fasta_kmer_file, 'w') as f:
            i = 0
            for seq, cnt in kmer_dict.items():
                if cnt == len(fasta_list):
                    f.write('>{}\n{}\n'.format(i, seq))
                    i += 1
