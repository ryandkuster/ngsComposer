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


if __name__ == '__main__':
    input1_list, input2_list, fastq_list, pairs_list = initialize(project_dir, paired)
    print(input1_list)
    print(input2_list)
    print(fastq_list)
    print(pairs_list)
    # if barcodes_file:
        # barcodes_matrix = barcode_reader(project_dir, barcodes_file)
        # barcode_test(barcodes_matrix, input1_list)
    # if front_trim > 0:
        # os.mkdir(project_dir + '/trim')
        # project_dir_current = project_dir + '/trim'
        # trim_part = partial(trimmer, front_trim, back_trim, project_dir_current)
        # pool = Pool(threads)        
        # pool.map(trim_part, fastq_list)
        # pool.close()
        # for i, filename in enumerate(input1_list):
            # input1_list[i] = project_dir_current + '/trimmed_' + os.path.basename(filename)
        # for i, filename in enumerate(input2_list):
            # input2_list[i] = project_dir_current + '/trimmed_' + os.path.basename(filename)
        # if barcodes_file:
            # for i in range(len(barcodes_matrix)):
                # barcodes_matrix[i][0] = 'trimmed_' + barcodes_matrix[i][0]


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
