import sys
import os
import gzip
import shutil


from multiprocessing import Pool
from functools import partial
from conf import *
from trimmer import *
from composer import *
from overhang import overhang


def initialize(project_dir, paired):
    '''
    check files in project_dir for completeness
    '''
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('project directory not found')
    fastq_list, gz, pairs_list = fastq_reader(project_dir)
    input1_list, input2_list = input_sort(paired, pairs_list)
    return input1_list, input2_list, fastq_list, pairs_list


def fastq_reader(project_dir):
    '''
    open fastq files, test first read structure, create file lists
    '''
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
    '''
    test first read structure if fastq
    '''
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
    '''
    test first read structure if fastq.gz
    '''
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
                sys.exit("sorry, gzipped functionality is not currently supported")


def is_paired(filename, line, pairs_list):
    '''
    use header info to match paired ends, if present
    '''
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
    '''
    if paired headers present, input1 and input2 lists ordered to keep pairs sequential
    '''
    input1_list, input2_list = [], []
    for values in pairs_list.values():
        if paired == True and len(values) == 2:
            input_paired(values, input1_list, input2_list)
        elif paired == True and len(values) != 2:
            sys.exit("paired forward and reverse reads don't match expected number" + "\n" +
                    "check the naming conventions of the barcodes_file")
        if paired == False:
            input_single(values, input1_list, input2_list)
    return input1_list, input2_list


def input_paired(values, input1_list, input2_list):
    '''
    form list if paired is True
    '''
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
    return input1_list, input2_list


def input_single(values, input1_list, input2_list):
    '''
    form list if paired is False, with user input if paired detection
    '''
    ignore = False
    if len(values) == 1:
        for filename in values:
            input1_list.append(filename)
    elif ignore == True:
        for filename in values:
            input1_list.append(filename)
    else:
        print("unexpected paired libraries found")
        answer = input("continue treating all files as single-end libraries?\n")
        ignore = True if answer in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()
        for filename in values:
            input1_list.append(filename)
    return input1_list, input2_list


def barcode_reader(barcodes_file):
    '''
    open user-defined barcode file and extract forward and reverse barcodes
    create matrix to pull corresponding sample ids that match barcodes
    '''
    barcodes_matrix = []
    with open(barcodes_file) as f:
        R2_barcodes = R2_barcodes_maker([], f)
        barcodes_matrix = array_maker(R2_barcodes)
        barcodes_matrix, R1_barcodes = barcodes_matrix_maker([], f, barcodes_matrix)
    return barcodes_matrix, R1_barcodes, R2_barcodes


def R2_barcodes_maker(R2_barcodes, f):
    '''
    populate list of R2 barcodes
    '''
    line = f.readline()
    for item in line.split():
        R2_barcodes.append(item)
    if len(R2_barcodes) == 1 and dual_index == True:
        sys.exit("expected barcodes file to have reverse barcodes")
    return R2_barcodes


def array_maker(R2_barcodes):
    '''
    create empty matrix with dimensions 1 by length R2 barcodes 
    '''
    barcodes_matrix = [[0] * 1 for i in range(len(R2_barcodes))]
    for i, x in enumerate(R2_barcodes):
        barcodes_matrix[i][0] = x
    return barcodes_matrix


def barcodes_matrix_maker(R1_barcodes, f, barcodes_matrix):
    '''
    create a matrix of user-defined sample ids and populate list of R1 barcodes
    '''
    for i, line in enumerate(f):
        for j, item in enumerate(line.split()):
            if j == 0:
                R1_barcodes.append(item)
            if j > 0:
                barcodes_matrix[j - 1].append(item)    
    return barcodes_matrix, R1_barcodes


def trim_muliproc(project_dir, threads, front_trim, back_trim, fastq_list):
    '''
    create user-defined number of subprocesses to trim every file in fastq_list
    '''
    project_dir_current = project_dir + '/trim'    
    os.mkdir(project_dir_current)
    trim_part = partial(trimmer, front_trim, back_trim, project_dir_current)
    pool = Pool(threads)        
    pool.map(trim_part, fastq_list)
    pool.close()
    for i, filename in enumerate(input1_list):
        input1_list[i] = project_dir_current + '/trimmed_' + os.path.basename(filename)
    for i, filename in enumerate(input2_list):
        input2_list[i] = project_dir_current + '/trimmed_' + os.path.basename(filename)
    # if barcodes_file:
        # for i in range(len(barcodes_matrix)):
            # barcodes_matrix[i][0] = 'trimmed_' + barcodes_matrix[i][0]    


def grater_multiproc():
    '''
    create user-defined number of subprocesses to demultiplex
    '''
    pass

if __name__ == '__main__':
    input1_list, input2_list, fastq_list, pairs_list = initialize(project_dir, paired)
    if barcodes_file:
        barcodes_file = project_dir + '/' + barcodes_file
        barcodes_matrix, R1_barcodes, R2_barcodes = barcode_reader(barcodes_file)
        print(R1_barcodes)
        print(R2_barcodes)
        print(barcodes_matrix)
    if front_trim > 0:
        trim_muliproc(project_dir, threads, front_trim, back_trim, fastq_list)
    # if barcodes_file:
        # os.mkdir(project_dir + '/demulti')
        # project_dir_current = project_dir + '/demulti'
        # if dual_index == True:
            # comp_part = partial(comp_pipeline_pe_di, input1_list, input2_list, mismatch, barcodes_matrix, project_dir_current)
        # TODO make a nice new function in composer that uses the barcodes as the file prefixes and doesn't duplicate barcodes
        # elif paired == True:
            # comp_part = partial(comp_pipeline_pe, input1_list, input2_list, mismatch, barcodes_matrix, project_dir_current)
        # if paired == False:
            # comp_part = partial(comp_pipeline_se, mismatch, barcodes_matrix, project_dir_current)    
        # pool = Pool(threads)
        # pool.map(comp_part, input1_list)
        # pool.close()

  
    # if overhang_list:
        # os.mkdir(project_dir + '/overhang')
        # inputs_list = []
        # try:
            # for filename in os.listdir(project_dir_current):
                # inputs_list.append(project_dir_current + '/' + filename)
        # except:
            # inputs_list = fastq_list
        # project_dir_current = project_dir + '/overhang'
        # hang_part = partial(overhang, project_dir_current, overhang_list)
        # pool = Pool(threads)
        # pool.map(hang_part, inputs_list)
        # pool.close()


    # shutil.rmtree(project_dir + '/demulti')
    # shutil.rmtree(project_dir + '/trim')
    # shutil.rmtree(project_dir + '/overhang')
