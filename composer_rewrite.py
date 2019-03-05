import sys
import os
import gzip
import shutil
from multiprocessing import Pool
from functools import partial


from tools.scallop import scallop_comp
from tools.anemone import anemone_comp
from tools.rotifer import rotifer_comp
from tools.scallop import scallop_end
from tools.krill import krill_comp
from tools.crinoid import crinoid_comp


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
            bcs_dict[key] = proj_dir + '/' + value
    return bcs_dict


def dir_test(proj_dir, bcs_dict):
    '''
    test project directory for correct structure
    '''
    fastq_ls = []
    for filename in os.listdir(proj_dir):
        if filename == os.path.basename(bcs_index):
            pass
        elif proj_dir + '/' + filename in bcs_dict.values():
            pass
        elif filename == 'conf.py' or filename == '__pycache__':
            pass
        else:
            fastq_ls.append(proj_dir + '/' + filename)
    return fastq_ls


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
            in1_ls, in2_ls = input_paired(values, in1_ls, in2_ls)
        elif paired is True and len(values) != 2:
            sys.exit("R1 and R2 pairs not found, expect paired library\n"
                     + "check the naming conventions of the bcs_index")
        if paired is False:
            in1_ls, in2_ls = input_single(values, in1_ls, in2_ls)
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


def dir_del(rm_dirs):
    for folder in rm_dirs:
        try:
            shutil.rmtree(folder)
            dir_name = os.path.basename(folder)
            print('\n composer is removing the ' + dir_name + ' directory')
        except FileNotFoundError:
            pass


def pool_multi(pool_part, pool_ls):
    pool = Pool(procs)
    pool.map(pool_part, pool_ls)
    pool.close()


def pathfinder(proj_dir_current, singles_ls, fastq_ls):
    '''
    walk current folder and pull files as lists
    '''
    singles_ls, fastq_ls, in1_ls, in2_ls = [], [], [], []
    for root, dirs, files in os.walk(os.path.abspath(proj_dir_current)):
        for i in files:
            fullname = os.path.join(root, i)
            if os.path.getsize(fullname) == 0:
                os.remove(fullname)
            elif root == str(proj_dir_current + '/single'):
                singles_ls.append(fullname)
            else:
                fastq_ls.append(fullname)
    return singles_ls, fastq_ls


def crinoid_multiproc(proj_dir, fastq_ls):
    '''
    create user-defined subprocesses to produce base-call summary
    '''
    proj_dir_current = proj_dir + '/qc'
    os.mkdir(proj_dir_current)
    crinoid_part = partial(crinoid_comp, proj_dir_current)
    pool_multi(crinoid_part, fastq_ls)


def scallop_muliproc(proj_dir, procs, front_trim, back_trim, fastq_ls, rm_dirs):
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    '''
    proj_dir_current = proj_dir + '/trimmed'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    scallop_part = partial(scallop_comp, front_trim, back_trim, proj_dir_current)
    pool_multi(scallop_part, fastq_ls)
    singles_ls, fastq_ls = pathfinder(proj_dir_current)

    #TODO make the below commands part of pathfinder
    singles_ls, fastq_ls = pathfinder(proj_dir_current, singles_ls, fastq_ls)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)


if __name__ == '__main__':
    proj_dir = initialize(sys.argv[1])
    sys.path.append(proj_dir)
    from conf import *

######################################################
#TODO delete the following:
######################################################
    old_dirs = [proj_dir + '/qc',
                proj_dir + '/trimmed',
                proj_dir + '/demultiplexed',
                proj_dir + '/parsed',
                proj_dir + '/end_trimmed',
                proj_dir + '/filtered']
    dir_del(old_dirs)
######################################################

    if bcs_index:
        bcs_index = proj_dir + '/' + bcs_index
        bcs_dict = index_reader(bcs_index)
    else:
        bcs_dict = {}
        bcs_index = ''

    fastq_ls = dir_test(proj_dir, bcs_dict)

    for filename in fastq_ls:
        try:
            fastq_test = is_fq(filename)
        except UnicodeDecodeError:
            fastq_test = is_gz(filename)
        if fastq_test is None:
            raise TypeError

    if bcs_index and paired is True and len(fastq_ls)/2 != len(bcs_dict):
        sys.exit('incorrect number of files based on index.txt')
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(paired, pairs_dict)
    rm_dirs = []


    if initial_qc is True:
        crinoid_multiproc(proj_dir, fastq_ls)

    if front_trim > 0:
        scallop_muliproc(proj_dir, procs, front_trim, 0, fastq_ls, rm_dirs)

    print(in1_ls)
    sys.exit()

    if bcs_index:
        fastq_ls, in1_ls, in2_ls = anemone_multiproc(
                walkthrough,
                proj_dir,
                mismatch,
                bcs_dict,
                in1_ls,
                in2_ls,
                rm_dirs)
    if bases_ls:
        fastq_ls, in1_ls, in2_ls, singles_ls = rotifer_multiproc(
                walkthrough,
                proj_dir,
                in1_ls,
                in2_ls,
                bases_ls,
                non_genomic,
                rm_dirs)
    if end_trim is True:
        try:
            fastq_ls, in1_ls, in2_ls, singles_ls = scallop_end_multiproc(
                bases_ls,
                fastq_ls,
                singles_ls,
                q_min,
                rm_dirs)
        except NameError:
            fastq_ls, in1_ls, in2_ls, singles_ls = scallop_end_multiproc(
                bases_ls,
                fastq_ls,
                [],
                q_min,
                rm_dirs)
    if q_min and q_percent:
        try:
            fastq_ls, in1_ls, in2_ls, singles_ls = krill_multiproc(
                walkthrough,
                in1_ls,
                in2_ls,
                singles_ls,
                q_min,
                q_percent,
                rm_dirs)
        except NameError:
            fastq_ls, in1_ls, in2_ls, singles_ls = krill_multiproc(
                walkthrough,
                in1_ls,
                in2_ls,
                [],
                q_min,
                q_percent,
                rm_dirs)
    
