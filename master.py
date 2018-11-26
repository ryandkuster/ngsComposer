import sys
import os
import gzip
from multiprocessing import Pool
from functools import partial
from conf import *
from trimmer import *
from composer import *
from overhang import overhang

def fastq_reader(project_dir):
    fastq_list, gz, pairs_list = [], [], {}
    for filename in os.listdir(project_dir):
        if filename == barcodes_file:
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


def barcode_reader(project_dir, barcodes_file):
    header = []
    try:
        with open(project_dir + '/' +  barcodes_file) as f:
            for i, line in enumerate(f):
                if i == 1:
                    barcodes_matrix = array_maker(header)
                for j, item in enumerate(line.split()):
                    if i == 0:
                        header.append(item)
                    else:
                        barcodes_matrix[j].append(item)
        return barcodes_matrix
    except FileNotFoundError:
        sys.exit('''
based on your configuration file your project
directory must contain a newline-separated list
of barcodes named ''' + barcodes_file
        )


def array_maker(header):
    barcodes_matrix = [[0] * 1 for i in range(len(header))]
    for i, x in enumerate(header):
        barcodes_matrix[i][0] = x
    return barcodes_matrix


def barcode_test(barcodes_matrix, input1_list):
    test_count = 0
    for x in input1_list:
        filename = os.path.basename(x)
        if dual_index == True:
            modifier = 2 
        else:
            modifier = 1
        for i, item in enumerate(barcodes_matrix):
            if item[0] == filename:
                test_count += 1
    if test_count != len(barcodes_matrix) - modifier:
        sys.exit('''
based on the configuration, the header of the barcodes file does not match
the fastq files contained in the project directory'''
        )


if __name__ == '__main__':
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('''
project directory not found
                ''')
    fastq_list, gz, pairs_list = fastq_reader(project_dir)
    input1_list, input2_list = input_sort(paired, pairs_list)
    if barcodes_file:
        barcodes_matrix = barcode_reader(project_dir, barcodes_file)
        barcode_test(barcodes_matrix, input1_list)
    if front_trim > 0:
        trim_part = partial(trimmer, front_trim, back_trim, project_dir)
        pool = Pool(threads)        
        pool.map(trim_part, fastq_list)
        pool.close()
        for i, filename in enumerate(input1_list):
            input1_list[i] = project_dir + '/trimmed_' + os.path.basename(filename)
        for i, filename in enumerate(input2_list):
            input2_list[i] = project_dir + '/trimmed_' + os.path.basename(filename)
    if barcodes_file:
        if paired == True:
            comp_part = partial(comp_piper, input1_list, input2_list, mismatch, R1_barcodes, project_dir)
        if paired == False:
            comp_part = partial(comp_piper_single, mismatch, R1_barcodes, project_dir)    
        pool = Pool(threads)
        pool.map(comp_part, input1_list)
        pool.close()
#        
#        '''
#        remove unwanted files
#        '''
#        tmp1, tmp2, = [], []
#        for x in range(len(R1_barcodes)):
#            for i, filename in enumerate(input1_list):
#                tmp1.append(project_dir + '/' + str(x + 1) + '_' + os.path.basename(filename))
#            for i, filename in enumerate(input2_list):
#                tmp2.append(project_dir + '/' + str(x + 1) + '_' + os.path.basename(filename))
#        if remove_intermediates == True and front_trim > 0:
#            for filename in input1_list:
#                os.remove(filename)
#            for filename in input2_list:
#                os.remove(filename)
#        if remove_fail == True:
#            fail1, fail2 = [], []
#            for filename in input1_list:
#                fail1.append(project_dir + '/unknown_' + os.path.basename(filename))
#            for filename in input2_list:
#                fail2.append(project_dir + '/unknown_' + os.path.basename(filename))
#            for filename in fail1:
#                os.remove(filename)
#            for filename in fail2:
#                os.remove(filename)
#            
#        input1_list, input2_list = tmp1, tmp2
#    if overhang_list:
#        hang_part = partial(overhang, project_dir, overhang_list)
#        pool = Pool(threads)
#        inputs_list = input1_list + input2_list
#        pool.map(hang_part, inputs_list)
#        pool.close()
#        
#        '''
#        remove unwanted files
#        '''
#        if remove_intermediates == True and R1_barcodes:
#            for item in input1_list:
#                os.remove(item)
#            for item in input2_list:
#                os.remove(item)
