import sys
import os
import gzip
import shutil


from multiprocessing import Pool
from functools import partial
from conf import *
from scallop import *
from anemone import *
from rotifer import *


# composer is:
# scallop - remove buffer sequences
# anemone - demultiplex based on barcodes
# krill - filter base calling with custom cutoffs
# rotifer - remove sequence artifacts from library prep


def index_reader(barcodes_index):
    '''
    open barcodes_index and create dictionary of associated barcode keyfiles
    '''
    barcodes_dict = {}
    with open(barcodes_index) as f:
        for line in f:
            for j, item in enumerate(line.split()):
                if j == 0:
                    key = item
                if j == 1:
                    value = item
                if j == 2:
                    sys.exit("barcodes_index should only contain forward reads" + "\n" +
                            "and their respective barcodes key file")
            barcodes_dict[key] = project_dir + '/' + value
    return barcodes_dict


def initialize(project_dir, paired, barcodes_dict):
    '''
    check files in project_dir for completeness
    '''
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('project directory not found')
    fastq_list, gz, pairs_list = fastq_reader(project_dir, barcodes_dict)
    input1_list, input2_list = input_sort(paired, pairs_list)
    return input1_list, input2_list, fastq_list, pairs_list


def fastq_reader(project_dir, barcodes_dict):
    '''
    bypass recognized barcode files, open fastqs, test first read structure, create file lists
    '''
    fastq_list, gz, pairs_list = [], [], {}
    for filename in os.listdir(project_dir):
        if filename == os.path.basename(barcodes_index):
            pass
        elif project_dir + '/' + filename in barcodes_dict.values():
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
                    "check the naming conventions of the barcodes_index match barcode keyfiles")
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


def trim_muliproc(project_dir, threads, front_trim, back_trim, fastq_list):
    '''
    create user-defined number of subprocesses to trim every file in fastq_list
    '''
    project_dir_current = project_dir + '/trim'    
    os.mkdir(project_dir_current)
    trim_part = partial(scallop_pipeline, front_trim, back_trim, project_dir_current)
    pool = Pool(threads)
    pool.map(trim_part, fastq_list)
    pool.close()
    for i, filename in enumerate(input1_list):
        input1_list[i] = project_dir_current + '/' + os.path.basename(filename)
    for i, filename in enumerate(input2_list):
        input2_list[i] = project_dir_current + '/' + os.path.basename(filename) 


def anemone_multiproc():
    '''
    create user-defined number of subprocesses to demultiplex
    '''
    project_dir_current = project_dir + '/demulti'
    os.mkdir(project_dir_current)
    comp_part = partial(anemone_pipeline, input1_list, input2_list, mismatch, barcodes_dict, project_dir_current)
    pool = Pool(threads)
    pool.map(comp_part, input1_list)
    pool.close()


if __name__ == '__main__':
    if barcodes_index:
        barcodes_index = project_dir + '/' + barcodes_index
        barcodes_dict = index_reader(barcodes_index)
    else:
        barcodes_dict = {}
        barcodes_index = ''
    input1_list, input2_list, fastq_list, pairs_list = initialize(project_dir, paired, barcodes_dict)
    #TODO add qc step here and let users know where to find its output
    if front_trim > 0:
        trim_muliproc(project_dir, threads, front_trim, 0, fastq_list)
    if barcodes_index:
        anemone_multiproc()

    # if overhang_list:
        # os.mkdir(project_dir + '/overhang')
        # inputs_list = []
        # try:
            # for filename in os.listdir(project_dir_current):
                # inputs_list.append(project_dir_current + '/' + filename)
        # except:
            # inputs_list = fastq_list
        # project_dir_current = project_dir + '/overhang'
        # hang_part = partial(rotifer, project_dir_current, overhang_list)
        # pool = Pool(threads)
        # pool.map(hang_part, inputs_list)
        # pool.close()


    # shutil.rmtree(project_dir + '/demulti')
    # print('\n composer is removing the dir, FYI \n')
    
    # shutil.rmtree(project_dir + '/trim')
    # shutil.rmtree(project_dir + '/overhang')
