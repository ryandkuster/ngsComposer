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
# anemone - demultiplex based on bcs
# krill - filter base calling with custom cutoffs
# rotifer - remove sequence artifacts from library prep


def index_reader(bcs_index):
    '''
    open bcs_index and create dictionary of associated bc keyfiles
    '''
    bcs_dict = {}
    with open(bcs_index) as f:
        for line in f:
            for j, item in enumerate(line.split()):
                if j == 0:
                    key = item
                if j == 1:
                    value = item
                if j == 2:
                    sys.exit("bcs_index should only contain forward reads" + "\n" +
                            "and their respective bcs key file")
            bcs_dict[key] = proj_dir + '/' + value
    return bcs_dict


def initialize(proj_dir, paired, bcs_dict):
    '''
    check files in proj_dir for completeness
    '''
    if os.path.exists(proj_dir) == True:
        proj_dir = os.path.abspath(proj_dir)
    else:
        sys.exit('project directory not found')
    fastq_ls, gz, pairs_ls = fastq_reader(proj_dir, bcs_dict)
    input1_ls, input2_ls = input_sort(paired, pairs_ls)
    return input1_ls, input2_ls, fastq_ls, pairs_ls


def fastq_reader(proj_dir, bcs_dict):
    '''
    bypass recognized bc files, open fastqs, test first read structure, create file lists
    '''
    fastq_ls, gz, pairs_ls = [], [], {}
    for filename in os.listdir(proj_dir):
        if filename == os.path.basename(bcs_index):
            pass
        elif proj_dir + '/' + filename in bcs_dict.values():
            pass
        else:
            try:
                fastq_test, pairs_ls = is_fq(proj_dir + '/' + filename, pairs_ls)
                gz.append(0)
            except UnicodeDecodeError:
                fastq_test, pairs_ls = is_gz(proj_dir + '/' + filename, pairs_ls)
                gz.append(1)
            if fastq_test is None:
                raise TypeError
            fastq_ls.append(proj_dir + '/' + filename)
    return fastq_ls, gz, pairs_ls


def is_fq(filename, pairs_ls):
    '''
    test first read structure if fastq
    '''
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            else:
                pairs_ls = is_paired(filename, line, pairs_ls)
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                return True, pairs_ls


def is_gz(filename, pairs_ls):
    '''
    test first read structure if fastq.gz
    '''
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            else:
                pairs_ls = is_paired(filename, line, pairs_ls)
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                sys.exit("sorry, gzipped functionality is not currently supported")


def is_paired(filename, line, pairs_ls):
    '''
    use header info to match paired ends, if present
    '''
    for i, x in enumerate(line):
        if x == ' ':
            space_pos = i
    header = line[:space_pos]
    if header in pairs_ls:
        pairs_ls[header].append(filename)
    else:
        pairs_ls[header] = [filename]
    return pairs_ls  


def input_sort(paired, pairs_ls):
    '''
    if paired headers present, input1 and input2 lists ordered to keep pairs sequential
    '''
    input1_ls, input2_ls = [], []
    for values in pairs_ls.values():
        if paired == True and len(values) == 2:
            input_paired(values, input1_ls, input2_ls)
        elif paired == True and len(values) != 2:
            sys.exit("paired forward and reverse reads don't match expected number" + "\n" +
                    "check the naming conventions of the bcs_index match bc keyfiles")
        if paired == False:
            input_single(values, input1_ls, input2_ls)
    return input1_ls, input2_ls


def input_paired(values, input1_ls, input2_ls):
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
                input1_ls.append(filename)
            if int(end) == 2:
                input2_ls.append(filename)
    return input1_ls, input2_ls


def input_single(values, input1_ls, input2_ls):
    '''
    form list if paired is False, with user input if paired detection
    '''
    ignore = False
    if len(values) == 1:
        for filename in values:
            input1_ls.append(filename)
    elif ignore == True:
        for filename in values:
            input1_ls.append(filename)
    else:
        print("unexpected paired libraries found")
        answer = input("continue treating all files as single-end libraries?\n")
        ignore = True if answer in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()
        for filename in values:
            input1_ls.append(filename)
    return input1_ls, input2_ls


def trim_muliproc(proj_dir, threads, front_trim, back_trim, fastq_ls):
    '''
    create user-defined number of subprocesses to trim every file in fastq_ls
    '''
    proj_dir_current = proj_dir + '/trim'    
    os.mkdir(proj_dir_current)
    trim_part = partial(scallop_pipeline, front_trim, back_trim, proj_dir_current)
    pool = Pool(threads)
    pool.map(trim_part, fastq_ls)
    pool.close()
    for i, filename in enumerate(input1_ls):
        input1_ls[i] = proj_dir_current + '/' + os.path.basename(filename)
    for i, filename in enumerate(input2_ls):
        input2_ls[i] = proj_dir_current + '/' + os.path.basename(filename) 


def anemone_multiproc():
    '''
    create user-defined number of subprocesses to demultiplex
    '''
    proj_dir_current = proj_dir + '/demulti'
    os.mkdir(proj_dir_current)
    comp_part = partial(anemone_pipeline, input1_ls, input2_ls, mismatch, bcs_dict, proj_dir_current)
    pool = Pool(threads)
    pool.map(comp_part, input1_ls)
    pool.close()
    for folder in os.listdir(proj_dir_current):
        for file in os.listdir(proj_dir_current + '/' + folder):
            if file.startswith('unknown.'):
                pass
            else:
                print(file)


if __name__ == '__main__':
    if bcs_index:
        bcs_index = proj_dir + '/' + bcs_index
        bcs_dict = index_reader(bcs_index)
    else:
        bcs_dict = {}
        bcs_index = ''
    input1_ls, input2_ls, fastq_ls, pairs_ls = initialize(proj_dir, paired, bcs_dict)
    #TODO add qc step here and let users know where to find its output
    if front_trim > 0:
        trim_muliproc(proj_dir, threads, front_trim, 0, fastq_ls)
    if bcs_index:
        anemone_multiproc()

    shutil.rmtree(proj_dir + '/trim')
    # shutil.rmtree(proj_dir + '/demulti')
    print('\n composer is removing the trim dir, FYI \n')

    # if overhang_ls:
        # os.mkdir(proj_dir + '/overhang')
        # inputs_ls = []
        # try:
            # for filename in os.listdir(proj_dir_current):
                # inputs_ls.append(proj_dir_current + '/' + filename)
        # except:
            # inputs_ls = fastq_ls
        # proj_dir_current = proj_dir + '/overhang'
        # hang_part = partial(rotifer, proj_dir_current, overhang_ls)
        # pool = Pool(threads)
        # pool.map(hang_part, inputs_ls)
        # pool.close()

    # shutil.rmtree(proj_dir + '/overhang')
