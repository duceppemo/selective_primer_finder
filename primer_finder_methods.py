#!/usr/local/env python3

from glob import glob
import subprocess
from multiprocessing import cpu_count
import sys
import gzip
import os
import pathlib
from concurrent import futures
from shutil import copyfileobj


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
        """
        Check if fasta file has a valid format
        :param fasta:
        :return:
        """
        status = 0

        with gzip.open(fasta, 'rt') if fasta.endswith('.gz') else open(fasta, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('>'):
                status = 1

        return status

    @staticmethod
    def check_cpus(requested_cpu):
        """

        :param requested_cpu:
        :return:
        """
        total_cpu = cpu_count()
        
        if requested_cpu:
            if requested_cpu > total_cpu:
                requested_cpu = total_cpu
                sys.stderr.write("Number of threads was set higher than available CPUs ({})".format(total_cpu))
                sys.stderr.write("Number of threads was set to {}".format(total_cpu))
        else:
            requested_cpu = total_cpu

        return requested_cpu

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
    def sourmash_compute(num_hashes, ksizes, output, cpu, input_fasta):
        cmd = ['sourmash', 'compute', '-n', str(num_hashes), '-p', str(cpu),
               '-k', str(ksizes), '-o', output, input_fasta]
        subprocess.run(cmd)

    @staticmethod
    def parallel_sourmash_compute(fasta_list, output_folder, cpu, kmer, sketch):
        # Create folder
        Methods.create_output_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((sketch, kmer, output_folder + '.'.join(os.path.basename(genome).split('.')[:-1]) + '.sig',
                   cpu, genome) for genome in fasta_list)
            # sourmash_compute(num_hashes, ksizes, output, cpu, input_fasta)
            for results in executor.map(lambda x: Methods.sourmash_compute(*x), args):
                pass

    @staticmethod
    def sourmash_signature_intersect(incl_sig_list, out_sig):
        """
        TODO -> log output
        :param incl_sig_list: list of signature files made from the inclusion group
        :param out_sig: file containing the common signatures
        :return:
        """

        cmd = ['sourmash', 'signature', 'intersect', '-o', out_sig] + incl_sig_list
        subprocess.run(cmd)

    @staticmethod
    def sourmash_signature_subtract(incl_sig, excl_sig_list, out_sig):
        """
        TODO -> log output
        :param incl_sig: file containing the common signatures
        :param excl_sig_list: list of signature files made from the exclusion group
        :param out_sig: file containing the inclusion-specific signatures
        :return:
        """
        cmd = ['sourmash', 'signature', 'subtract', '-o', out_sig, incl_sig] + excl_sig_list
        subprocess.run(cmd)

    @staticmethod
    def sourmash_signature_to_kmer(out_csv, in_sig, genome):
        """
        Fetch the sequences linked to the signatures
        :param out_csv:
        :param in_sig:
        :param genome:
        :return:
        """
        # Where sourmash is installed
        # install_path = which('sourmash')
        # sig2kmer = install_path + '/utils/signature-to-kmers.py'
        sig2kmer = '/home/bioinfo/prog/sourmash/utils/signature-to-kmers.py'
        cmd = ['python3', sig2kmer, '--output-kmers', out_csv, in_sig, genome]
        subprocess.run(cmd)

    @staticmethod
    def parse_cvs(csv_file, seq_dict):
        with open(csv_file, 'r') as f:
            f.readline()  # skip header
            for line in f:
                line = line.rstrip()
                field_list = line.split(',')
                seq_dict[field_list[1]] = field_list[0]

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
    def jellyfish_count(seq_file, count_file, cpu):
        kmer_size = 21
        hash_size = '100M'

        cmd = ['jellyfish', 'count', '-m', str(kmer_size), '-s', hash_size, '-t', str(cpu), '-C',
               '-o', count_file, seq_file]

        subprocess.run(cmd)

    @staticmethod
    def parallel_jellyfish_count(fasta_list, output_folder, cpu):
        # Create folder
        Methods.create_output_folders(output_folder)

        with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
            args = ((genome, output_folder + '.'.join(os.path.basename(genome).split('.')[:-1]) + '.jf',
                     1) for genome in fasta_list)
            for results in executor.map(lambda x: Methods.jellyfish_count(*x), args):
                pass

    @staticmethod
    def jellyfish_merge(jf_list, folder):
        cmd = ['jellyfish', 'merge', '-o', folder + '/exclusion_merged.jf'] + jf_list
        subprocess.run(cmd)

        # Delete individial jellyfish database files
        for f in jf_list:
            os.remove(f)

    # @staticmethod
    # def jellyfish_query(ident, kmer, jf_db, folder):
    #     cmd = ['jellyfish', 'query', '-o', folder + '/' + ident + '.cnt', jf_db, kmer]
    #     results = subprocess.run(cmd, stdout=subprocess.PIPE)
    #
    #     return tuple(results.stdout)

    @staticmethod
    def jellyfish_query(kmer, jf_db):
        cmd = ['jellyfish', 'query', jf_db, kmer]
        results = subprocess.run(cmd, stdout=subprocess.PIPE)
        cnt = results.stdout.rstrip().decode('ascii').split(' ')[1]
        return cnt

    # @staticmethod
    # def parallel_jellyfish_query(kmer_list, output_folder, cpu):
    #     with futures.ThreadPoolExecutor(max_workers=cpu) as executor:
    #         args = ((ident, kmer, output_folder + '.'.join(os.path.basename(genome).split('.')[:-1]) + '.jf')
    #                 for ident, kemr in fasta_list)
    #         for results in executor.map(lambda x: Methods.jellyfish_count(*x), args):
    #             pass

    @staticmethod
    def sliding_window(sequence, win_size, step=1):
        """
        https://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/

        Returns a generator that will iterate through
        the defined chunks of input sequence.  Input sequence
        must be iterable.
        """

        # Verify the inputs
        try:
            it = iter(sequence)
        except TypeError:
            raise Exception("**ERROR** sequence must be iterable.")
        if not ((type(win_size) == type(0)) and (type(step) == type(0))):
            raise Exception("**ERROR** type(winSize) and type(step) must be int.")
        if step > win_size:
            raise Exception("**ERROR** step must not be larger than winSize.")
        if win_size > len(sequence):
            raise Exception("**ERROR** winSize must not be larger than sequence length.")

        # Pre-compute number of chunks to emit
        num_chunks = int(((len(sequence) - win_size) / step) + 1)

        # Do the work
        for i in range(0, num_chunks * step, step):
            yield sequence[i:i + win_size]

    @staticmethod
    def find_unique_regions(seq_dict):
        import difflib
        for ident, seq in seq_dict.items():
            pass

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
            'SEQUENCE_INTERNAL_EXCLUDED_REGION':'37,21',
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
