import argparse
import gzip
import os
import shutil
import subprocess
import sys
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


def conf_confirm(proj_dir):
    assert cfg.paired == True or cfg.paired == False
    assert cfg.procs >= 1
    assert os.path.exists(cfg.alt_dir) or cfg.alt_dir == False
    assert cfg.initial_qc == True or cfg.initial_qc == False
    assert cfg.walkthrough == True or cfg.walkthrough == False
    assert cfg.walkaway == True or cfg.walkaway == False
    assert isinstance(cfg.front_trim, int)
    assert os.path.exists(proj_dir + '/' + cfg.bcs_index) or\
            cfg.bcs_index == False
    assert isinstance(cfg.mismatch, int)
    assert cfg.R1_bases_ls == False or isinstance(cfg.R1_bases_ls, list)
    if isinstance(cfg.R1_bases_ls, list):
        for i in cfg.R1_bases_ls:
            for j in i:
                assert j in ['A', 'C', 'G', 'T']
    assert cfg.R2_bases_ls == False or isinstance(cfg.R2_bases_ls, list)
    if isinstance(cfg.R2_bases_ls, list):
        for i in cfg.R2_bases_ls:
            for j in i:
                assert j in ['A', 'C', 'G', 'T']
    assert cfg.non_genomic == False or isinstance(cfg.non_genomic, int)
    assert 0 <= cfg.q_min <= 40
    assert 0 <= cfg.q_percent <= 100
    if cfg.q_min and cfg.q_percent:
        pass
    else:
        sys.exit('both q_min and q_percent variables must be defined')
    assert cfg.end_trim == True or cfg.end_trim == False
    assert cfg.rm_transit == True or cfg.rm_transit == False


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
    #TODO add test for presence of all files in barcodes file
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


def bc_test(bcs_file):
    '''
    open user-defined bc file and extract forward and reverse
    bcs and test for redundancy
    '''
    #TODO handle identical barcodes...
    with open(bcs_file) as f:
        line = f.readline()
        R1_bcs, R2_bcs, R1_nest, R2_nest = {}, {}, {}, {}
        for i, item in enumerate(line.split()):
            if item in R2_bcs.keys():
                sys.exit(bcs_file + ' contains duplicated R2 barcodes')
            else:
                R2_bcs[item] = i
        for i, line in enumerate(f):
            for j, item in enumerate(line.split()):
                if j == 0:
                    if item in R1_bcs.keys():
                        sys.exit(bcs_file + ' contains duplicated R1 barcodes')
                    else:
                        R1_bcs[item] = i
    for i, (k1, v1) in enumerate(R1_bcs.items()):
        for j, (k2, v2) in enumerate(R1_bcs.items()):
            if k2.startswith(k1) and v1 != v2:
                if k1 in R1_nest.keys():
                    R1_nest[k1].append(k2)
                else:
                    R1_nest[k1] = [k1, k2]
#                sys.exit('redundancy detected with barcodes ' + k2 + ' and ' +
#                        k1 + ' in file ' + bcs_file)
    for i, (k1, v1) in enumerate(R2_bcs.items()):
        for j, (k2, v2) in enumerate(R2_bcs.items()):
            if k2.startswith(k1) and v1 != v2:
                if k1 in R2_nest.keys():
                    R2_nest[k1].append(k2)
                else:
                    R2_nest[k1] = [k1, k2]
#                sys.exit('redundancy detected with barcodes ' + k2 + ' and ' +
#                        k1 + ' in file ' + bcs_file)
    return R1_nest, R2_nest

def is_fq(filename):
    '''
    test first read structure if fastq
    '''
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 0 and line[0] != '@':
                return
            if i == 2 and line[0] != '+':
                return
            else:
                return 1


def is_gz(filename):
    '''
    test first read structure if fastq.gz
    '''
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            if i == 0 and line[0] != '@':
                return
            if i == 2 and line[0] != '+':
                return
            else:
                return 0


def gz_header(filename):
    try:
        fastq_test = is_fq(filename)
        with open(filename) as f:
            header = f.readline()
            return header
    except:
        with gzip.open(filename, 'rt') as f:
            header = f.readline()
            return header


def is_paired(fastq_ls):
    '''
    use header info to match paired ends, if present
    '''
    pairs_dict = {}
    for filename in fastq_ls:
        header = gz_header(filename)
        for i, x in enumerate(header.split(' ')):
            if i == 0:
                header_id = x
        if header_id in pairs_dict:
            pairs_dict[header_id].append(filename)
        else:
            pairs_dict[header_id] = [filename]
    return pairs_dict


def input_sort(pairs_dict):
    '''
    if headers identical, in1 and in2 lists ordered to keep pairing
    '''
    in1_ls, in2_ls = [], []
    for values in pairs_dict.values():
        if cfg.paired is True and len(values) == 2:
            in1_ls, in2_ls = input_paired(values, in1_ls, in2_ls)
        elif cfg.paired is True and len(values) != 2:
            sys.exit('R1 and R2 pairs not found, expect paired library\n'
                     + 'check the naming conventions of the bcs_index')
        if cfg.paired is False:
            in1_ls, in2_ls = input_single(values, in1_ls, in2_ls)
    return in1_ls, in2_ls


def input_paired(values, in1_ls, in2_ls):
    '''
    form list if paired is True
    '''
    for filename in values:
        header = gz_header(filename)
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
        print('unexpected paired libraries found')
        ans = input('continue treating files as single-end libraries?\n')
        ignore = True if ans in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()
        for filename in values:
            in1_ls.append(filename)
    return in1_ls, in2_ls


def dir_size(proj_dir, fastq_ls, fastq_test):
    drive_stats = os.statvfs(proj_dir)
    drive_free = drive_stats.f_frsize * drive_stats.f_bavail
    dir_used = 0
    dir_plan = 0
    dir_plan = dir_plan + 1 if cfg.front_trim else dir_plan
    dir_plan = dir_plan + 1 if cfg.bcs_index else dir_plan
    if cfg.R1_bases_ls:
        dir_plan += 1
    elif cfg.R2_bases_ls:
        dir_plan += 1
    dir_plan = dir_plan + 1 if cfg.end_trim else dir_plan
    dir_plan = dir_plan + 1 if cfg.q_min else dir_plan
    for i in fastq_ls:
        dir_used += os.path.getsize(i)
    if cfg.rm_transit:
        dir_plan = dir_used * 2 if fastq_test else dir_used * 2 * 5
    else:
        dir_plan = dir_used * dir_plan if fastq_test\
                else dir_used * dir_plan * 5
    if dir_plan >= drive_free:
        sys.exit('an estimated ' + str(dir_plan) + ' bytes are required to ' +
                'process, consider rm_transit or alt_dir variables')


def dir_del(rm_dirs):
    '''
    deletes specified list of folders
    '''
    for folder in rm_dirs:
        for root, dirs, files in os.walk(os.path.abspath(folder)):
            if 'qc' in dirs:
                dirs.remove('qc')
            for i in files:
                fullname = os.path.join(root, i)
                os.remove(fullname)


def pool_multi(pool_part, pool_ls):
    '''
    create n subprocesses of current tool
    '''
    pool = Pool(cfg.procs)
    pool.map(pool_part, pool_ls)
    pool.close()


def pathfinder(proj_dir_current):
    '''
    walk current directory and pull files as lists
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
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(pairs_dict)
    return singles_ls, fastq_ls, in1_ls, in2_ls


def concater(proj_dir_current):
    '''
    walk current directory and concatenate files with identical names
    '''
    concat_dict, cat_ls = {}, []
    for root, dirs, files in os.walk(os.path.abspath(proj_dir_current)):
        for i in files:
            fullname = os.path.join(root, i)
            if i.startswith('unknown.'):
                os.remove(fullname)
            else:
                if os.path.getsize(fullname) == 0:
                    os.remove(fullname)
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


def walkaway_opt(folder):
    ans = input('\n check qc folder(s) in ' + folder + ' - continue? (y/n)\n')
    ignore = True if ans in ('Y', 'y', 'Yes', 'yes', 'YES') else sys.exit()


def crinoid_multiproc(proj_dir, fastq_ls):
    '''
    create user-defined subprocesses to produce base-call summary
    '''
    proj_dir_current = proj_dir + '/qc'
    os.mkdir(proj_dir_current)
    crinoid_part = partial(crinoid_comp, proj_dir_current)
    pool_multi(crinoid_part, fastq_ls)


def scallop_muliproc(proj_dir, fastq_ls):
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    '''
    proj_dir_current = proj_dir + '/trimmed'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    scallop_part = partial(scallop_comp, cfg.front_trim, None,
            proj_dir_current)
    pool_multi(scallop_part, fastq_ls)
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(proj_dir_current)
    return singles_ls, fastq_ls, in1_ls, in2_ls


def anemone_multiproc(proj_dir, bcs_dict, in1_ls, in2_ls):
    '''
    create user-defined subprocesses to demultiplex
    '''
    proj_dir_current = proj_dir + '/demultiplexed'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    anemone_part = partial(anemone_comp, in1_ls, in2_ls, cfg.mismatch,
            bcs_dict, proj_dir_current)
    pool_multi(anemone_part, in1_ls)
    concater(proj_dir_current)
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(proj_dir_current)
    if cfg.walkthrough:
        crinoid_multiproc(proj_dir_current, fastq_ls)
        cfg.walkaway == False and walkaway_opt(rm_dirs[-1])
    if cfg.rm_transit is True:
        dir_del(rm_dirs[:-1])
    return singles_ls, fastq_ls, in1_ls, in2_ls


def rotifer_multiproc(proj_dir, in1_ls, in2_ls):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    proj_dir_current = proj_dir + '/parsed'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    os.mkdir(proj_dir_current + '/single')
    os.mkdir(proj_dir_current + '/paired')
    rotifer_part = partial(rotifer_comp, in1_ls, in2_ls, cfg.R1_bases_ls,
            cfg.R2_bases_ls, cfg.non_genomic, proj_dir_current)
    pool_multi(rotifer_part, in1_ls)
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(proj_dir_current)
    if cfg.walkthrough:
        crinoid_multiproc(proj_dir_current + '/single', singles_ls)
        crinoid_multiproc(proj_dir_current + '/paired', fastq_ls)
        cfg.walkaway == False and walkaway_opt(rm_dirs[-1])
    if cfg.rm_transit is True:
        dir_del(rm_dirs[:-1])
    return singles_ls, fastq_ls, in1_ls, in2_ls


def scallop_end_multiproc(fastq_ls, singles_ls):
    '''
    automated 3' end read trimming based on q_min value
    '''
    #TODO adjust for initial_qc == False and walkthrough == False...
    proj_dir_current = proj_dir + '/end_trimmed'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    if cfg.R1_bases_ls and cfg.R2_bases_ls:
        os.mkdir(proj_dir_current + '/single')
        os.mkdir(proj_dir_current + '/paired')
    if cfg.walkthrough is False and rm_dirs[-2] != proj_dir + '/qc':
        try:
            crinoid_multiproc(rm_dirs[-2] + '/single', singles_ls)
            crinoid_multiproc(rm_dirs[-2] + '/paired', fastq_ls)
        except FileNotFoundError:
            crinoid_multiproc(rm_dirs[-2], fastq_ls)
    scallop_end_part = partial(scallop_end, proj_dir_current, cfg.q_min)
    pool_multi(scallop_end_part, fastq_ls)
    if singles_ls:
        scallop_end_part = partial(scallop_end, proj_dir_current, cfg.q_min)
        pool_multi(scallop_end_part, singles_ls)
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(proj_dir_current)
    if cfg.rm_transit is True:
        dir_del(rm_dirs[:-1])
    return singles_ls, fastq_ls, in1_ls, in2_ls


def krill_multiproc(in1_ls, in2_ls, singles_ls):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    proj_dir_current = proj_dir + '/filtered'
    rm_dirs.append(proj_dir_current)
    os.mkdir(proj_dir_current)
    os.mkdir(proj_dir_current + '/single')
    os.mkdir(proj_dir_current + '/single/pe_lib')
    os.mkdir(proj_dir_current + '/single/se_lib')
    os.mkdir(proj_dir_current + '/paired')
    krill_part = partial(krill_comp, in1_ls, in2_ls, cfg.q_min, cfg.q_percent,
            proj_dir_current)
    pool_multi(krill_part, in1_ls)
    if singles_ls:
        krill_part = partial(krill_comp, in1_ls, in2_ls, cfg.q_min,
                cfg.q_percent, proj_dir_current)
        pool_multi(krill_part, singles_ls)
    concater(proj_dir_current + '/single')
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(proj_dir_current)
    shutil.rmtree(proj_dir_current + '/single/pe_lib')
    shutil.rmtree(proj_dir_current + '/single/se_lib')
    if cfg.walkthrough:
        crinoid_multiproc(proj_dir_current + '/single', singles_ls)
        crinoid_multiproc(proj_dir_current + '/paired', fastq_ls)
    if cfg.rm_transit is True:
        dir_del(rm_dirs[:-1])
    return singles_ls, fastq_ls, in1_ls, in2_ls


if __name__ == '__main__':
    singles_ls, fastq_ls, in1_ls, in2_ls, rm_dirs = [], [], [], [], []
    parser = argparse.ArgumentParser(description=
            'composer is a base-call error-filtering and read preprocessing\
            pipeline for fastq libraries - \
            see https://github.com/ryandkuster/composer for usage')
    parser.add_argument('-i', type=str,
            help='the full or relative path to the project directory')
    args = parser.parse_args()
    proj_dir = initialize(args.i)
    sys.path.append(proj_dir)
    import conf as cfg
    conf_confirm(proj_dir)

######################################################
# TODO delete the following (for ease of testing):
######################################################
    old_dirs = [proj_dir + '/qc',
                proj_dir + '/trimmed',
                proj_dir + '/demultiplexed',
                proj_dir + '/parsed',
                proj_dir + '/end_trimmed',
                proj_dir + '/filtered']
    for folder in old_dirs:
        try:
            shutil.rmtree(folder)
            dir_name = os.path.basename(folder)
            print('\n composer is removing the ' + dir_name + ' directory')
        except FileNotFoundError:
            pass
######################################################

    if cfg.bcs_index:
        bcs_index = proj_dir + '/' + cfg.bcs_index
        bcs_dict = index_reader(bcs_index)
    else:
        bcs_dict = {}
        bcs_index = ''

    fastq_ls = dir_test(proj_dir, bcs_dict)

    for i in bcs_dict.values():
        R1_nest, R2_nest = bc_test(i)

    for filename in fastq_ls:
        try:
            fastq_test = is_fq(filename)
        except UnicodeDecodeError:
            fastq_test = is_gz(filename)
        if fastq_test is None:
            sys.exit(filename + ' was not expected in project directory')

    if cfg.alt_dir:
        proj_dir = initialize(cfg.alt_dir)
        if len(os.listdir(proj_dir)) != 0:
            sys.exit('alt_dir must be an empty directory')

    dir_size(proj_dir, fastq_ls, fastq_test)

    if cfg.bcs_index and cfg.paired is True and \
            len(fastq_ls)/2 != len(bcs_dict):
        sys.exit('incorrect number of files based on index.txt')
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(pairs_dict)

    '''
    begin calling tools
    '''
    if cfg.initial_qc is True:
        crinoid_multiproc(proj_dir, fastq_ls)

    if cfg.front_trim > 0:
        singles_ls, fastq_ls, in1_ls, in2_ls = scallop_muliproc(proj_dir,
                fastq_ls)

    if cfg.bcs_index:
        singles_ls, fastq_ls, in1_ls, in2_ls = anemone_multiproc(proj_dir,
                bcs_dict, in1_ls, in2_ls)

    if cfg.R1_bases_ls and cfg.R2_bases_ls:
        singles_ls, fastq_ls, in1_ls, in2_ls = rotifer_multiproc(proj_dir,
                in1_ls, in2_ls)

    if cfg.end_trim is True:
        singles_ls, fastq_ls, in1_ls, in2_ls = scallop_end_multiproc(fastq_ls,
                singles_ls)

    if cfg.q_min and cfg.q_percent:
        singles_ls, fastq_ls, in1_ls, in2_ls = krill_multiproc(in1_ls, in2_ls,
                singles_ls)
