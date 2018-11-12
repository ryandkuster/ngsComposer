import sys
import os
import gzip
from multiprocessing import Pool
from functools import partial
from conf import *
from trimmer import *
from composer import *


def fastq_reader(project_dir):
    fastq_list, gz, pairs_list = [], [], {}
    for filename in os.listdir(project_dir):
        if filename == R1_barcodes:
            pass
        elif filename == R2_barcodes:
            pass
        else:
            try:
                fastq_test, pairs_list = is_fq(project_dir + '/' + filename, pairs_list)
                gz.append(0)
            except UnicodeDecodeError:
                fastq_test, pairs_list = is_gz(project_dir + '/' + filename, pairs_list)
                gz.append(1)
            if fastq_test is None:
                raise TypeError
            fastq_list.append(project_dir + '/' + filename)
    return fastq_list, gz, pairs_list


def is_fq(filename, pairs_list):
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            else:
                pairs_list = is_paired(filename, line, pairs_list)
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                return True, pairs_list


def is_gz(filename, pairs_list):
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            else:
                pairs_list = is_paired(filename, line, pairs_list)
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                sys.exit('''
sorry, gzipped functionality is not currently supported
            ''')
#                    return True, pairs_list


def is_paired(filename, line, pairs_list):
    for i, x in enumerate(line):
        if x == ' ':
            space_pos = i
    header = line[:space_pos]
    if header in pairs_list:
        pairs_list[header].append(filename)
    else:
        pairs_list[header] = [filename]
    return pairs_list


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


def input_sort(paired, pairs_list):
    input1_list, input2_list, ignore = [], [], False
    for values in pairs_list.values():
        if paired == True:
            if len(values) == 2:
                for filename in values:
                    with open(filename) as f:
                        header = f.readline()
                        for i, x in enumerate(header):
                            if x == ' ':
                                space_pos = i
                        end = header[space_pos+1]
                        if int(end) == 1:
                            input1_list.append(filename)
                        if int(end) == 2:
                            input2_list.append(filename)
            else:
                sys.exit('''
paired forward and end reads don't match expected number
                    ''')
        if paired == False:
            if len(values) == 1:
                for filename in values:
                    input1_list.append(filename)
            elif ignore == True:
                for filename in values:
                    input1_list.append(filename)
            else:
                print('''
unexpected paired libraries found
                    ''')
                answer = input('continue treating all files as single-end libraries?\n')
                ignore = True if answer in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()
                for filename in values:
                    input1_list.append(filename)
    return input1_list, input2_list


if __name__ == '__main__':
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('''
project directory not found
                ''')
    fastq_list, gz, pairs_list = fastq_reader(project_dir)
    input1_list, input2_list = input_sort(paired, pairs_list)
    if R1_barcodes:
        R1_barcodes = barcode_reader(project_dir, R1_barcodes)
    if R2_barcodes:
        R2_barcodes = barcode_reader(project_dir, R2_barcodes)
    if front_trim > 0:
        trim_part = partial(trimmer, front_trim, back_trim, project_dir)
        pool = Pool(threads)        
        pool.map(trim_part, fastq_list)
        pool.close()
        for i, filename in enumerate(input1_list):
            input1_list[i] = project_dir + '/trimmed_' + os.path.basename(filename)
        for i, filename in enumerate(input2_list):
            input2_list[i] = project_dir + '/trimmed_' + os.path.basename(filename)
    if R1_barcodes:
#        comp_part = partial(comp_piper, input1_list, input2_list, mismatch, R1_barcodes, project_dir)
#        pool = Pool(threads)        
#        pool.map(comp_part, input1_list)
#        pool.close()
        tmp1, tmp2 = [], []
        for x in range(len(R1_barcodes)):
            for i, filename in enumerate(input1_list):
                tmp1.append(project_dir + '/' + str(x + 1) + '_' + os.path.basename(filename))
            for i, filename in enumerate(input2_list):
                tmp2.append(project_dir + '/' + str(x + 1) + '_' + os.path.basename(filename))
        input1_list, input2_list = tmp1, tmp2

