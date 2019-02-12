import sys
import os
import gzip
import shutil
from multiprocessing import Pool
from functools import partial


from conf import *


from scallop import scallop_comp
from anemone import anemone_comp
from rotifer import rotifer_comp
from krill import krill_comp
from porifera import porifera_comp
from crinoid import crinoid_comp


def initialize(proj_dir):
    '''
    check files in proj_dir for completeness
    '''
    if os.path.exists(proj_dir) is True:
        proj_dir = os.path.abspath(proj_dir)
    else:
        sys.exit('project directory not found')
    return proj_dir


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
                    sys.exit("bcs_index should only contain forward read(s)\n" +
                             "and their respective barcodes key file(s) tab\n" +
                             "in tab-delimited format")
            bcs_dict[key] = proj_dir + '/' + value
    return bcs_dict


def fastq_reader(proj_dir, bcs_dict):
    '''
    bypass recognized bc files, open fastqs, test first read structure
    create file lists and test for pairs
    '''
    fastq_ls, gz = [], []
    for filename in os.listdir(proj_dir):
        if filename == os.path.basename(bcs_index):
            pass
        elif proj_dir + '/' + filename in bcs_dict.values():
            pass
        else:
            try:
                fastq_test = is_fq(proj_dir + '/' + filename)
                gz.append(0)
            except UnicodeDecodeError:
                fastq_test = is_gz(proj_dir + '/' + filename)
                gz.append(1)
            if fastq_test is None:
                raise TypeError
            fastq_ls.append(proj_dir + '/' + filename)
    return fastq_ls, gz


def is_fq(filename):
    '''
    test first read structure if fastq
    '''
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                return True


def is_gz(filename):
    '''
    test first read structure if fastq.gz
    '''
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i == 1 and line[0] != '@':
                return
            if i == 3 and line[0] != '+':
                return
            if i == 5 and line[0] != '@':
                return
            else:
                sys.exit("sorry, .gz functionality is not currently supported")


def is_paired(fastq_ls):
    '''
    use header info to match paired ends, if present
    '''
    pairs_dict = {}
    for filename in fastq_ls:
        with open(filename) as f:
            header = f.readline()
            for i, x in enumerate(header.split(' ')):
                if i == 0:
                    header_id = x
            if header_id in pairs_dict:
                pairs_dict[header_id].append(filename)
            else:
                pairs_dict[header_id] = [filename]
    return pairs_dict


def input_sort(paired, pairs_dict):
    '''
    if headers identical, in1 and in2 lists ordered to keep pairing
    '''
    in1_ls, in2_ls = [], []
    for values in pairs_dict.values():
        if paired is True and len(values) == 2:
            input_paired(values, in1_ls, in2_ls)
        elif paired is True and len(values) != 2:
            sys.exit("paired R1 and R2 reads don't match expected number\n"
                     + "check the naming conventions of the bcs_index")
        if paired is False:
            input_single(values, in1_ls, in2_ls)
    return in1_ls, in2_ls


def input_paired(values, in1_ls, in2_ls):
    '''
    form list if paired is True
    '''
    for filename in values:
        with open(filename) as f:
            header = f.readline()
            for i, x in enumerate(header.split(' ')):
                if i == 1:
                    end = x
            if end.startswith('1'):
                in1_ls.append(filename)
            if end.startswith('2'):
                in2_ls.append(filename)
    return in1_ls, in2_ls


def input_single(values, in1_ls, in2_ls):
    '''
    form list if paired is False, with user input if paired detection
    '''
    ignore = False
    if len(values) == 1:
        for filename in values:
            in1_ls.append(filename)
    elif ignore is True:
        for filename in values:
            in1_ls.append(filename)
    else:
        print("unexpected paired libraries found")
        ans = input("continue treating files as single-end libraries?\n")
        ignore = True if ans in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()
        for filename in values:
            in1_ls.append(filename)
    return in1_ls, in2_ls


def concater(concat_dict, proj_dir_current, cat_ls):
    for root, dirs, files in os.walk(os.path.abspath(proj_dir_current)):
        for i in files:
            fullname = os.path.join(root, i)
            if i.startswith('unknown.'):
                pass
            else:
                if os.path.getsize(fullname) == 0:
                    os.remove(fullname)
#                    pass
                elif i in concat_dict:
                    concat_dict[i].append(fullname)
                else:
                    concat_dict[i] = [fullname]
    for filename in concat_dict.keys():
        with open(proj_dir_current + '/' + filename, 'w') as o1:
            cat_ls.append(o1.name)
            for i in concat_dict[filename]:
                with open(i) as obj1:
                    shutil.copyfileobj(obj1, o1)
                os.remove(i)
    return cat_ls


def pathfinder(proj_dir_current, singles_ls, fastq_ls):
    for root, dirs, files in os.walk(os.path.abspath(proj_dir_current)):
        for i in files:
            fullname = os.path.join(root, i)
            if os.path.getsize(fullname) == 0:
                os.remove(fullname)
#                pass
            elif root == str(proj_dir_current + '/single'):
                singles_ls.append(fullname)
            else:
                fastq_ls.append(fullname)
    return singles_ls, fastq_ls


def scallop_muliproc(proj_dir, procs, front_trim, back_trim, fastq_ls):
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    '''
    proj_dir_current = proj_dir + '/trimmed'
    os.mkdir(proj_dir_current)
    trim_part = partial(scallop_comp, front_trim, back_trim, proj_dir_current)
    pool = Pool(procs)
    pool.map(trim_part, fastq_ls)
    pool.close()
    for i, filename in enumerate(in1_ls):
        in1_ls[i] = proj_dir_current + '/' + os.path.basename(filename)
    for i, filename in enumerate(in2_ls):
        in2_ls[i] = proj_dir_current + '/' + os.path.basename(filename)


def anemone_multiproc(proj_dir, mismatch, bcs_dict, in1_ls, in2_ls):
    '''
    create user-defined subprocesses to demultiplex
    '''
    proj_dir_current = proj_dir + '/demultiplexed'
    os.mkdir(proj_dir_current)
    anemone_part = partial(anemone_comp, in1_ls, in2_ls, mismatch, bcs_dict,
                        proj_dir_current)
    pool = Pool(procs)
    pool.map(anemone_part, in1_ls)
    pool.close()
    concat_dict, fastq_ls, in1_ls, in2_ls = {}, [], [], []
    fastq_ls = concater(concat_dict, proj_dir_current, fastq_ls)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)
    return in1_ls, in2_ls


def rotifer_multiproc(proj_dir, in1_ls, in2_ls, bases_ls, non_genomic):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    proj_dir_current = proj_dir + '/parsed'
    os.mkdir(proj_dir_current)
    os.mkdir(proj_dir_current + '/single')
    os.mkdir(proj_dir_current + '/paired')
    rotifer_part = partial(rotifer_comp, in1_ls, in2_ls, bases_ls, non_genomic, proj_dir_current)
    pool = Pool(procs)
    pool.map(rotifer_part, in1_ls)
    pool.close()
    singles_ls, fastq_ls, in1_ls, in2_ls = [], [], [], []
    singles_ls, fastq_ls = pathfinder(proj_dir_current, singles_ls, fastq_ls)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)
    return in1_ls, in2_ls, singles_ls


def krill_multiproc(in1_ls, in2_ls, singles_ls, q_min, q_percent):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    proj_dir_current = proj_dir + '/filtered'
    os.mkdir(proj_dir_current)
    os.mkdir(proj_dir_current + '/single')
    os.mkdir(proj_dir_current + '/single/pe_lib')
    os.mkdir(proj_dir_current + '/single/se_lib')
    os.mkdir(proj_dir_current + '/paired')
    krill_part = partial(krill_comp, in1_ls, in2_ls, q_min, q_percent, proj_dir_current)
    pool = Pool(procs)
    pool.map(krill_part, in1_ls)
    pool.close()
    if singles_ls:
        krill_part = partial(krill_comp, in1_ls, in2_ls, q_min, q_percent, proj_dir_current)
        pool = Pool(procs)
        pool.map(krill_part, singles_ls)
        pool.close()
    concat_dict, singles_ls, fastq_ls, in1_ls, in2_ls = {}, [], [], [], []
    _ = concater(concat_dict, proj_dir_current + '/single', [])
    singles_ls, fastq_ls = pathfinder(proj_dir_current, singles_ls, fastq_ls)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)
    return in1_ls, in2_ls, singles_ls


if __name__ == '__main__':
    proj_dir = initialize(proj_dir)
    if bcs_index:
        bcs_index = proj_dir + '/' + bcs_index
        bcs_dict = index_reader(bcs_index)
    else:
        bcs_dict = {}
        bcs_index = ''
    fastq_ls, gz = fastq_reader(proj_dir, bcs_dict)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)

    # TODO add qc step here and let users know where to find its output

    if front_trim > 0:
        scallop_muliproc(proj_dir, procs, front_trim, 0, fastq_ls)
    if bcs_index:
        in1_ls, in2_ls = anemone_multiproc(proj_dir, mismatch, bcs_dict, in1_ls, in2_ls)
    if bases_ls:
        in1_ls, in2_ls, singles_ls = rotifer_multiproc(proj_dir, in1_ls, in2_ls, bases_ls, non_genomic)
    if q_min and q_percent:
        in1_ls, in2_ls, singles_ls = krill_multiproc(in1_ls, in2_ls, singles_ls, q_min, q_percent)




#    try:
#        shutil.rmtree(proj_dir + '/trimmed')
#    except FileNotFoundError:
#        pass

#    try:
#        shutil.rmtree(proj_dir + '/demultiplexed')
#    except FileNotFoundError:
#        pass

#    try:
#        shutil.rmtree(proj_dir + '/parsed')
#    except FileNotFoundError:
#        pass

#    try:
#        shutil.rmtree(proj_dir + '/filtered')
#    except FileNotFoundError:
#        pass

    print('\n composer is removing directories, FYI \n')









