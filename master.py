import sys
import os
import gzip
from multiprocessing import Pool
from functools import partial

from conf import *

from trimmer import trimmer
import composer

def fastq_reader(project_dir):
    fastq_list, gz, pairs = [], [], {}
    for filename in os.listdir(project_dir):
        if filename == R1_barcodes:
            pass
        elif filename == R2_barcodes:
            pass
        else:
            try:
                fastq_test, pairs = file_type.is_fq(project_dir + '/' + filename, pairs)
                gz.append(0)
            except UnicodeDecodeError:
                fastq_test, pairs = file_type.is_gz(project_dir + '/' + filename, pairs)
                gz.append(1)
            if fastq_test is None:
                raise TypeError
            fastq_list.append(project_dir + '/' + filename)
    return fastq_list, gz, pairs

def barcode_reader(project_dir, barcodes):
    barcode_list = []
    try:
        with open(project_dir + '/' + barcodes) as f:
            for line in f:
                barcode_list.append(line.rstrip())
            return barcode_list
    except FileNotFoundError:
        sys.exit('''
        based on your configuration file your project
        directory must contain a newline-separated list
        of barcodes named ''' + barcodes
        )

class file_type:
    def is_fq(filename, pairs):
        with open(filename) as f:
            for i, line in enumerate(f):
                if i == 1 and line[0] != '@':
                    return
                else:
                    pairs = file_type.is_paired(filename, line, pairs)
                if i == 3 and line[0] != '+':
                    return
                if i == 5 and line[0] != '@':
                    return
                else:
                    return True, pairs

    def is_gz(filename, pairs):
        with gzip.open(filename, 'rt') as f:
            for i, line in enumerate(f):
                if i == 1 and line[0] != '@':
                    return
                else:
                    pairs = file_type.is_paired(filename, line, pairs)
                if i == 3 and line[0] != '+':
                    return
                if i == 5 and line[0] != '@':
                    return
                else:
                    sys.exit('''
                sorry, gzipped functionality is not currently supported
                ''')
#                    return True, pairs

    def is_paired(filename, line, pairs):
        for i, x in enumerate(line):
            if x == ' ':
                space_pos = i
        header = line[:space_pos]
        if header in pairs:
            pairs[header].append(filename)
        else:
            pairs[header] = [filename]
        return pairs

#def pipe_writer_forward1(input1, input2, cutoff, chunk, barcodes, project_dir):
#    output1 = input1
#    output2 = input2
#    input1 = project_dir + '/trimmer_' + input1
#    input2 = project_dir + '/trimmer_' + input2
#    composer.line_writer(input1, input2, cutoff, chunk, barcodes,
#                        project_dir, output1, output2)

if __name__ == '__main__':
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('''
                project directory not found
                ''')
    fastq_list, gz, pairs = fastq_reader(project_dir)
    if R1_barcodes:
        R1_barcodes = barcode_reader(project_dir, R1_barcodes)
    if R2_barcodes:
        R2_barcodes = barcode_reader(project_dir, R2_barcodes)
    print(fastq_list)
    print(gz)
    print(pairs)
    print(R1_barcodes)


#    if trim == True:
#        for x, input_file in enumerate(fastq_list):
#            output_file =  project_dir + '/trimmer_' + os.path.basename(input_file)
#            trimmer.test(input_file, b, e, output_file, gz[x])


#    input1_list, input2_list = [], []
#    for values in pairs.values():
#        if len(values) == 2:
#            for filename in values:
#                with open(filename) as f:
#                    header = f.readline()
#                    for pos, x in enumerate(header):
#                        if x == ' ':
#                            space_pos = pos
#                    end = header[space_pos+1]
#                    if int(end) == 1:
#                        input1_list.append(os.path.basename(filename))
#                    if int(end) == 2:
#                        input2_list.append(os.path.basename(filename))
#        else:
#            print('is this a single end library?')


#    print(input1_list)
#    print(input2_list)
#    if R1_barcodes == True:
#        for x, input1 in enumerate(input1_list):
#            input1 = input1
#            input2 = input2_list[x]
#            pipe_writer_forward1(input1, input2, cutoff, 200000, R1_barcodes, project_dir)
