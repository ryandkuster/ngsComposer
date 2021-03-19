import argparse
import datetime
import gzip
import os
import shutil
import subprocess
import sys
from argparse import RawTextHelpFormatter
from functools import partial
from multiprocessing import Pool

sys.path.append(os.path.join(os.path.dirname(__file__),'tools'))

import helpers.messages as msg
from anemone import anemone_comp, bc_reader, bc_test
from crinoid import combine_matrix, crinoid_comp
from helpers.gzip_handling import gzip_test
from krill import krill_comp
from porifera import porifera_comp, reverse_comp
from rotifer import rotifer_comp
from scallop import scallop_comp
from subprocess import check_call


class Project:
    def __init__(self):
        self.start_time =  str(datetime.datetime.now()).split('.')[0]
        self.singles_ls = []
        self.fastq_ls = []
        self.fastq_dt = {}
        self.in1_ls = []
        self.in2_ls = []
        self.rm_dirs = []
        self.ignore = False
        self.bypass = False

        self.paired = False
        self.procs = 1
        self.alt_dir = False
        self.initial_qc = False
        self.all_qc = False
        self.walkaway = True
        self.front_trim = 0
        self.end_trim = False
        self.bcs_index = ''
        self.mismatch = False
        self.R1_bases_ls = []
        self.R2_bases_ls = []
        self.non_genomic = False
        self.end_score = False
        self.window = False
        self.min_len = 1
        self.adapters1 = ''
        self.adapters2 = ''
        self.adapter_match = False
        self.tiny = False
        self.q_min = False
        self.q_percent = False
        self.rm_transit = True
        self.p64 = False
        self.compress = True
        self.r_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                'tools', 'helpers', 'R_packages')


    def initialize(self, proj):
        '''
        check files in project directory for completeness
        '''
        if os.path.exists(proj) is True:
            self.proj = os.path.abspath(proj)
            sys.path.append(self.proj)
            import conf as config
            for i in config.__dict__.keys():
                if i not in c.__dict__.keys() and not i.startswith('__'):
                    print('didn\'t expect ' + i)
                    sys.exit()
            c.__dict__.update(config.__dict__)
        else:
            sys.exit(msg.initialize1)

    def index_reader(self):
        '''
        open bcs_index and create dictionary of associated bc keyfiles
        open adapters and bcs_index and test for correct format
        '''
        adapters1 = os.path.join(self.proj, 'adapters.R1.txt')
        adapters2 = os.path.join(self.proj, 'adapters.R2.txt')
        self.adapters1 = adapters1 if os.path.exists(adapters1) else ''
        self.adapters2 = adapters2 if os.path.exists(adapters2) else ''
        if self.adapters1:
            with open(self.adapters1) as f:
                adapters_ls = [line.rstrip() for line in f]
                nucleotide_test(adapters_ls)
        if self.adapters2:
            with open(self.adapters2) as f:
                adapters_ls = [line.rstrip() for line in f]
                nucleotide_test(adapters_ls)

        bc_path = os.path.join(self.proj, 'index.txt')
        self.bcs_index = bc_path if os.path.exists(bc_path) else ''
        self.bcs_dict = {}
        if not c.bcs_index:
            return
        with open(self.bcs_index) as f:
            for line in f:
                for j, item in enumerate(line.split()):
                    if j == 0:
                        key = item
                    if j == 1:
                        value = item
                self.bcs_dict[key] = os.path.join(self.proj, value)

    def conf_confirm(self):
        '''
        test user-input from configuration file
        '''
        if not isinstance(self.paired, bool):
            raise Exception(msg.conf_confirm1)

        if self.procs < 1 or not isinstance(self.procs, int):
            raise Exception(msg.conf_confirm2)

        if self.alt_dir and not os.path.exists(self.alt_dir):
            raise Exception(msg.conf_confirm3)

        if not isinstance(self.initial_qc, bool):
            raise Exception(msg.conf_confirm4)

        if self.all_qc and self.all_qc not in ['full', 'summary']:
            raise Exception(msg.conf_confirm5)

        if not isinstance(self.walkaway, bool):
            raise Exception(msg.conf_confirm6)

        if self.walkaway is False and not self.all_qc:
            raise Exception(msg.conf_confirm7)

        if self.front_trim and not isinstance(self.front_trim, int):
            raise Exception(msg.conf_confirm8)
        if self.front_trim and self.front_trim < 1:
            raise Exception(msg.conf_confirm9)
        if self.mismatch:
            if not os.path.exists(os.path.join(self.proj, 'index.txt')):
                raise Exception(msg.conf_confirm10)
            if not isinstance(self.mismatch, int) or self.mismatch < 0:
                raise Exception(msg.conf_confirm11)

        if self.R1_bases_ls:
            if not isinstance(self.R1_bases_ls, list):
                raise Exception(msg.conf_confirm12)
            nucleotide_test(self.R1_bases_ls)
        if self.R2_bases_ls:
            if not isinstance(self.R2_bases_ls, list):
                raise Exception(msg.conf_confirm13)
            nucleotide_test(self.R2_bases_ls)
        if self.non_genomic:
            if not isinstance(self.non_genomic, int) or self.non_genomic < 1:
                raise Exception(msg.conf_confirm14)

        if self.adapters1:
            if not self.adapter_match:
                raise Exception(msg.conf_confirm23)
        if self.adapter_match:
            if not isinstance(self.adapter_match, int):
                raise Exception(msg.conf_confirm24)
            if self.adapter_match < 10:
                raise Exception(msg.conf_confirm24)
            if not os.path.exists(os.path.join(self.proj, 'adapters.R2.txt')):
                raise Exception(msg.conf_confirm23)

        if self.q_min or self.q_percent:
            if not self.q_min or not self.q_percent:
                raise Exception(msg.conf_confirm18)
            if not isinstance(self.q_min, int):
                raise Exception(msg.conf_confirm19)
            if not 0 <= self.q_min <= 42:
                raise Exception(msg.conf_confirm19)
            if not isinstance(self.q_percent, int):
                raise Exception(msg.conf_confirm20)
            if not 0 <= self.q_percent <= 100:
                raise Exception(msg.conf_confirm20)

        if not isinstance(self.rm_transit, bool):
            raise Exception(msg.conf_confirm21)

        if not isinstance(self.p64, bool):
            raise Exception(msg.conf_confirm25)

    def dir_test(self):
        '''
        test project directory for correct structure
        '''
        self.fastq_ls = []
        for filename in os.listdir(self.proj):
            if filename == os.path.basename(self.bcs_index):
                pass
            elif filename == os.path.basename(self.adapters1):
                pass
            elif filename == os.path.basename(self.adapters2):
                pass
            elif os.path.join(self.proj, filename) in self.bcs_dict.values():
                pass
            elif filename == 'conf.py' or filename == '__pycache__':
                pass
            else:
                self.fastq_ls.append(os.path.join(self.proj, filename))


def r_packages():
    '''
    test R and packages and install dependencies if needed
    '''
    try:
        subprocess.check_call(['Rscript', os.path.join(os.path.dirname(
            __file__), 'tests', 'test_packages.R')] + [c.r_dir], shell=False)
    except FileNotFoundError:
        sys.exit(msg.r_packages1)


def nucleotide_test(ls):
    '''
    test list of sequences for proper formatting
    bypass is used for manual input during walkthrough mode
    '''
    test = msg.nucleotide_test1
    nuc_ls = (j for i in ls for j in i)
    for i in nuc_ls:
        if i in ['A', 'C', 'G', 'T']:
            test = True
            pass
        elif i.upper() in ['A', 'C', 'G', 'T']:
            test = msg.nucleotide_test2
            break
        else:
            test = msg.nucleotide_test3
            break
    if test is not True:
        print(test)
        if not c.bypass:
            sys.exit()
        raise Exception


def fastq_test(fastq_ls):
    '''
    test if gzipped fastq file
    '''
    fastq_dt, in1_ls, in2_ls = {}, [], []
    for filename in fastq_ls:
        compressed = gzip_test(filename)
        if compressed is None:
            sys.exit('\n\n' + filename + msg.fastq_test1)
        elif compressed is True:
            with gzip.open(filename, 'rt') as f:
                fastq_dt = fastq_structure(f, filename, fastq_dt)
        else:
            with open(filename) as f:
                fastq_dt = fastq_structure(f, filename, fastq_dt)

    for i in fastq_dt.values():
        if None not in i:
            in1_ls.append(i[0])
            in2_ls.append(i[1])
    return in1_ls, in2_ls


def fastq_structure(f, filename, fastq_dt):
    '''
    test first read for fastq structure, extract headers, find pairs
    '''
    for i, line in enumerate(f):
        if i == 0:
            try:
                header = line.rstrip().split(' ')[0]
                read_no = int(line.rstrip().split(' ')[1][0])
            except(IndexError, TypeError, ValueError) as e:
                try:
                    header = line.rstrip()[:-1]
                    read_no = int(line.rstrip().split('/')[-1])
                except(IndexError, TypeError, ValueError) as e:
                    sys.exit('expected fastq headers not found in ' +
                             os.path.basename(filename))
            if header in fastq_dt:
                fastq_dt[header][read_no -1] = filename
            else:
                fastq_dt[header] = [None, None]
                fastq_dt[header][read_no -1] = filename
            if line[0] != '@':
                sys.exit('\n\n' + filename + msg.fastq_test1)
        if i == 2 and line[0] != '+':
            sys.exit('\n\n' + filename + msg.fastq_test1)
        else:
            return fastq_dt


def demultiplex_test():
    '''
    test for correct structure of barcodes file(s)
    '''
    for k, v in c.bcs_dict.items():
        try:
            assert os.path.join(c.proj, k) in c.fastq_ls
        except AssertionError:
            sys.exit('check index.txt formatting, ' + k + ' not found')
        names_mx, R1_bcs, R2_bcs, dual_index = bc_reader(v)
        test = bc_test(R1_bcs, names_mx, True)
        if test:
            print('redundant R1 barcodes detected in ' + v)
        test = bc_test(R2_bcs, names_mx, False)
        if test:
            print('redundant R2 barcodes detected in ' + v)
        nucleotide_test(R1_bcs.values())
        if len(R2_bcs) > 1:
            nucleotide_test(R2_bcs.values())
            for k2, v2 in c.fastq_dt.items():
                if v2[0] == os.path.join(c.proj, k):
                    try:
                        assert v2[1] is not None
                    except AssertionError:
                        sys.exit('R2 for ' + k + ' expect but not found')


def adapters_test():
    '''
    assess adapters for presence of barcodes
    '''
    if c.adapters1 and c.bcs_dict:
        with open(c.adapters1) as f1:
            adapters_ls1 = [line.rstrip() for line in f1]
    else:
        return

    if c.paired is True:
        try:
            with open(c.adapters2) as f2:
                adapters_ls2 = [line.rstrip() for line in f2]
        except FileNotFoundError:
            sys.exit('please provide \'adapters.R2.txt\' containing adapter ' +
                     'sequences expected in R2 sequences.')
    else:
        return

    for k, v in c.bcs_dict.items():
        names_mx, R1_bcs, R2_bcs, dual_index = bc_reader(v)
        for i in R2_bcs.values():
            R2_test = [j for j in adapters_ls1 if i in j]
            if len(R2_test) == 0:
                print('barcode ' + i + ' not found in adapters.R1.txt')
                while True:
                    choice = input(msg.adapters_test1)
                    if choice == '1':
                        return
                    elif choice == '2':
                        sys.exit()
                    else:
                        print('\nselect an option from the list\n')

        for i in R1_bcs.values():
            R1_test = [j for j in adapters_ls2 if i in j]
            if len(R1_test) == 0:
                print('barcode ' + i + ' not found in adapters.R2.txt')
                while True:
                    choice = input(msg.adapters_test1)
                    if choice == '1':
                        return
                    elif choice == '2':
                        sys.exit()
                    else:
                        print('\nselect an option from the list\n')


def dir_size(proj, fastq_ls):
    '''
    test the specified output directory for adequate disk space
    '''
    drive_stats = os.statvfs(proj)
    drive_free = drive_stats.f_frsize * drive_stats.f_bavail
    dir_used = 0
    dir_plan = 0
    dir_plan = dir_plan + 1 if c.front_trim else dir_plan
    dir_plan = dir_plan + 1 if c.bcs_index else dir_plan
    dir_plan = dir_plan + 1 if c.R1_bases_ls or c.R2_bases_ls else dir_plan
    dir_plan = dir_plan + 1 if c.end_score else dir_plan
    dir_plan = dir_plan + 1 if c.adapters1 else dir_plan
    dir_plan = dir_plan + 1 if c.q_min else dir_plan

    for i in fastq_ls:
        if i[-2:] == 'gz':
            dir_used += os.path.getsize(i) * 5
        else:
            dir_used += os.path.getsize(i)

    dir_plan = dir_used * 2 if c.rm_transit is True else dir_used * dir_plan
    if dir_plan >= drive_free:
        while True:
            choice = input(msg.dir_size1 + str(dir_plan) + msg.dir_size2)
            if choice == '1':
                break
            elif choice == '2':
                sys.exit()
            else:
                print('\nselect an option from the list\n')


def dir_make(title):
    '''
    create new directory when pipeline tools called
    '''
    curr = os.path.join(c.proj, str(len(c.rm_dirs)) + '_' + title)
    c.rm_dirs.append(curr)
    os.mkdir(curr)
    return curr


def dir_del(rm_dirs):
    '''
    delete specified list of folders
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
    pool = Pool(c.procs)
    pool.map(pool_part, pool_ls)
    pool.close()


def pathfinder(curr):
    '''
    walk current directory and pull files as 4 distinct lists:
    - singles_ls = files stored in 'singles' directory
    - fastq_ls = files not in 'singles' directory
    - in1_ls = ordered R1 files with same index as in2_ls
    - in2_ls = ordered R2 files with same index as in1_ls
    singles_ls is a necessary evil to avoid identity issues
    '''
    fastq_ls, singles_ls = [], []
    for root, dirs, files in os.walk(os.path.abspath(curr)):
        for i in files:
            fullname = os.path.join(root, i)
            if i.startswith(('unknown.', 'redundant.')):
                pass
            elif os.path.getsize(fullname) == 0:
                os.remove(fullname)
            elif root == str(os.path.join(curr, 'single')):
                singles_ls.append(fullname)
            else:
                fastq_ls.append(fullname)
    in1_ls, in2_ls = fastq_test(fastq_ls)
    for i in fastq_ls:
        if i not in in1_ls and i not in in2_ls:
            singles_ls.append(i)
    return [singles_ls, fastq_ls, in1_ls, in2_ls]


def compress_files(curr):
    files_ls = []
    for root, dirs, files in os.walk(os.path.abspath(curr)):
        for i in files:
            fullname = os.path.join(root, i)
            if os.path.getsize(fullname) == 0:
                pass
            else:
                files_ls.append(fullname) 
    compress_part = partial(compress_multi)
    pool_multi(compress_part, files_ls)


def compress_multi(filename):
    check_call(['gzip', filename])


def concater(curr):
    '''
    walk current directory and concatenate files with identical names
    files from subdirectories are concatenated in the containing folder
    and removed one by one
    '''
    concat_dict = {}
    for root, dirs, files in os.walk(os.path.abspath(curr)):
        for i in files:
            fullname = os.path.join(root, i)
            if i.startswith(('unknown.', 'redundant.')):
                pass
            else:
                if os.path.getsize(fullname) == 0:
                    os.remove(fullname)
                elif i in concat_dict:
                    concat_dict[i].append(fullname)
                else:
                    concat_dict[i] = [fullname]

    for filename, filepaths in concat_dict.items():
        if len(filepaths) == 1:
            shutil.move(filepaths[0], os.path.join(curr, filename))
        else:
            with open(os.path.join(curr, filename), 'w') as o1:
                for i in filepaths:
                    with open(i) as obj1:
                        shutil.copyfileobj(obj1, o1)
                    os.remove(i)


def crinoid_multi(proj, ls):
    '''
    create user-defined subprocesses to produce base-call summary
    '''
    curr = os.path.join(proj, 'qc')
    os.mkdir(curr)
    all_qc = 'full' if proj == c.proj else c.all_qc
    crinoid_part = partial(crinoid_comp, curr, all_qc, 1, c.p64, c.r_dir)
    pool_multi(crinoid_part, ls)


def scallop_multi():
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    takes fastq_ls as input, as no pair order is necessary
    '''
    curr = dir_make('trimmed')
    scallop_part = partial(scallop_comp, c.in1_ls, c.in2_ls, c.front_trim, None, None, None, None, curr)
    pool_multi(scallop_part, c.fastq_ls)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    return temp_ls


def anemone_multi():
    '''
    create user-defined subprocesses to demultiplex
    takes fastq_ls as primary input and checks each file for the
    possibility of a pair
    '''
    curr = dir_make('demultiplexed')
    anemone_part = partial(anemone_comp, c.in1_ls, c.in2_ls, c.mismatch,
                           c.bcs_dict, curr, c.front_trim)
    pool_multi(anemone_part, c.fastq_ls)
    concater(curr)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, anemone_multi, temp_ls,
                              mismatch=c.mismatch)
    return temp_ls


def paired_setup(curr):
    '''
    avoids conflict of writing identically named files to the
    single directory by keeping previously single files (se_lib)
    and newly-written single files (pe_lib) separate
    '''
    os.mkdir(os.path.join(curr, 'single'))
    os.mkdir(os.path.join(curr, 'single', 'pe_lib'))
    os.mkdir(os.path.join(curr, 'single', 'se_lib'))
    os.mkdir(os.path.join(curr, 'paired'))


def paired_takedown(curr):
    '''
    combine identically named files in pe_lib and se_lib and remove
    these directories
    '''
    concater(os.path.join(curr, 'single'))
    shutil.rmtree(os.path.join(curr, 'single', 'pe_lib'))
    shutil.rmtree(os.path.join(curr, 'single', 'se_lib'))


def rotifer_multi():
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    curr = dir_make('parsed')
    paired_setup(curr)
    rotifer_part = partial(rotifer_comp, c.in1_ls, c.in2_ls, c.R1_bases_ls,
                           c.R2_bases_ls, c.non_genomic, curr)
    pool_multi(rotifer_part, c.in1_ls)
    if c.singles_ls:
        rotifer_part = partial(rotifer_comp, [], [], c.R1_bases_ls,
                               c.R2_bases_ls, c.non_genomic, curr)
        pool_multi(rotifer_part, c.singles_ls)
    elif not c.in1_ls:
        pool_multi(rotifer_part, c.fastq_ls)
    paired_takedown(curr)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, rotifer_multi, temp_ls,
                              R1_bases_ls=c.R1_bases_ls,
                              R2_bases_ls=c.R2_bases_ls,
                              non_genomic=c.non_genomic)
    return temp_ls


def scallop_end_multi():
    '''
    create user-defined subprocesses to remove low-scoring 3' ends
    '''
    curr = dir_make('end_trimmed')
    paired_setup(curr)
    scallop_part = partial(scallop_comp, c.in1_ls, c.in2_ls, None, None,
                           c.end_score, c.window, c.min_len, curr)
    pool_multi(scallop_part, c.in1_ls)
    if c.singles_ls:
        scallop_part = partial(scallop_comp, [], [], None, None, c.end_score,
                               c.window, c.min_len, curr)
        pool_multi(scallop_part, c.singles_ls)
    elif not c.in1_ls:
        pool_multi(scallop_part, c.fastq_ls)
    paired_takedown(curr)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, scallop_end_multi, temp_ls,
                              end_score=c.end_score,
                              window=c.window,
                              min_len=c.min_len)
    return temp_ls


def porifera_multi():
    '''
    create user-defined subprocesses to detect and remove adapter sequences
    '''
    curr = dir_make('adapted')
    paired_setup(curr)
    if c.R1_bases_ls or c.R2_bases_ls:
        search = len(max(c.R1_bases_ls + c.R2_bases_ls, key = len)) + c.front_trim
    else:
        search = c.front_trim
    porifera_part = partial(porifera_comp, curr, c.in1_ls, c.in2_ls,
                            c.adapters1, c.adapters2, c.bcs_dict,
                            search, c.adapter_match, c.min_len, c.tiny,
                            c.R1_bases_ls, c.R2_bases_ls)
    pool_multi(porifera_part, c.in1_ls)
    if c.singles_ls:
        porifera_part = partial(porifera_comp, curr, [], [], c.adapters1,
                                c.adapters2, c.bcs_dict, search,
                                c.adapter_match, c.min_len, c.tiny,
                                c.R1_bases_ls, c.R2_bases_ls)
        pool_multi(porifera_part, c.singles_ls)
    elif not c.in1_ls:
        pool_multi(porifera_part, c.fastq_ls)
    paired_takedown(curr)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, porifera_multi, temp_ls,
                              adapter_match=c.adapter_match,
                              min_len=c.min_len)
    return temp_ls


def krill_multi():
    '''
    create user-defined subprocesses to parse based on quality scores
    '''
    curr = dir_make('filtered')
    paired_setup(curr)
    krill_part = partial(krill_comp, c.in1_ls, c.in2_ls, c.q_min, c.q_percent,
                         c.p64, curr)
    pool_multi(krill_part, c.in1_ls)
    if c.singles_ls:
        krill_part = partial(krill_comp, [], [], c.q_min, c.q_percent, c.p64, curr)
        pool_multi(krill_part, c.singles_ls)
    elif not c.in1_ls:
        pool_multi(krill_part, c.fastq_ls)
    paired_takedown(curr)
    if c.compress:
        compress_files(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, krill_multi, temp_ls,
                              q_min=c.q_min,
                              q_percent=c.q_percent)
    return temp_ls


def walkthrough(curr, tool, temp_ls, **kwargs):
    '''
    query user for modifying or accepting current step in pipeline
    '''
    c.bypass = True # erroneous input won't halt pipeline (raises exception)

    try:
        crinoid_multi(os.path.join(curr, 'single'), temp_ls[0])
        crinoid_multi(os.path.join(curr, 'paired'), temp_ls[1])
    except FileNotFoundError:
        crinoid_multi(curr, temp_ls[1])

    if len(temp_ls[2]) >= 1:
        combine_matrix(temp_ls[2], c.r_dir, 'R1_summary.csv')
    if len(temp_ls[3]) >= 1:
        combine_matrix(temp_ls[3], c.r_dir, 'R2_summary.csv')
    if len(temp_ls[0]) >= 1:
        combine_matrix(temp_ls[0], c.r_dir, 'singles_summary.csv')

    if c.walkaway:
        return temp_ls

    while True:
        status = msg.walkthrough2 if not c.walkaway else msg.walkthrough3
        print('\nwalkthrough mode is ' + status)
        print('\nplease check the qc folders in ' + c.rm_dirs[-1])
        print('\nvariables to modify:')
        for k, v in kwargs.items():
            print('\n' + k + ' : ' + str(v))
        choice1 = input(msg.walkthrough1)
        if choice1 == '1':
            return temp_ls
        elif choice1 == '3':
            shutil.rmtree(c.rm_dirs[-1])
            c.rm_dirs.pop()
            for k, v in kwargs.items():
                v = False
                setattr(c, k, v)
            return [c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls]
        elif choice1 == '4':
            c.walkaway = True if c.walkaway is False else False
        elif choice1 == '5':
            summary_file()
            sys.exit('\nngs-composer is now exiting')
        elif choice1 == '2':
            break
        else:
            print('\nselect an option from the list\n')

    while True:
        for k, v in kwargs.items():
            print('\ncurrent value for ' + k + ' is ' + str(v))
            v_new = input('press enter to keep current value or input new ' +
                          'value for ' + k + '? > ')
            if v_new.startswith('['):
                v_new = v_new.strip('[\']').split(', ')
            elif v_new == '':
                v_new = v
            elif v_new == 'False':
                v_new = False
            elif any(i.isdigit() for i in v_new):
                v_new = int(v_new)
            setattr(c, k, v_new)
        try:
            c.conf_confirm()
            choice2 = input(msg.walkthrough4)
            if choice2 == '1':
                shutil.rmtree(c.rm_dirs[-1])
                c.rm_dirs.pop()
                print('\nrerunning step with updated variables...')
                temp_ls = tool()
                return temp_ls
            else:
                print('\nselect an option from the list\n')
        except Exception:
            print('\nnot a valid entry\n')


def tidy_up():
    '''
    walk the complete directory and remove empty subdirectories
    '''
    rm_list, empty_dir = [], False
    try:
        shutil.rmtree(os.path.join(c.proj, '__pycache__'))
    except FileNotFoundError:
        pass
    for root, dirs, files in os.walk(os.path.abspath(c.proj)):
        for i in dirs:
            if len(os.listdir(os.path.join(root, i))) == 0:
                empty_dir = True
                rm_list.append(os.path.join(root, i))
    for i in rm_list:
        shutil.rmtree(i)
    if empty_dir is True:
        tidy_up()
    else:
        return


def summary_file():
    end_time = str(datetime.datetime.now()).split('.')[0]
    log = ('ngscomposer version ' + version + '\n' +
           'see https://github.com/ryandkuster/ngscomposer/releases '\
           'for newest release info\n\n' +
           'start ' + c.start_time + '\n' +
           'end   ' + end_time + '\n\n' +
           'paired = ' + str(c.paired) + '\n' +
           'procs = ' + str(c.procs) + '\n' +
           'alt_dir = ' + str(c.alt_dir) + '\n' +
           'walkaway = ' + str(c.walkaway) + '\n' +
           'rm_transit = ' + str(c.rm_transit) + '\n' +
           'initial_qc = ' + str(c.initial_qc) + '\n' +
           'all_qc = ' + str(c.all_qc) + '\n' +
           'front_trim = ' + str(c.front_trim) + '\n' +
           'bcs_index = ' + str(c.bcs_index) + '\n' +
           'mismatch = ' + str(c.mismatch) + '\n' +
           'R1_bases_ls = ' + str(c.R1_bases_ls) + '\n' +
           'R2_bases_ls = ' + str(c.R2_bases_ls) + '\n' +
           'non_genomic = ' + str(c.non_genomic) + '\n' +
           'end_score = ' + str(c.end_score) + '\n' +
           'window = ' + str(c.window) + '\n' +
           'min_len = ' + str(c.min_len) + '\n' +
           'R1 searched for: ' + str(c.adapters1) + '\n' +
           'R2 searched for: ' + str(c.adapters2) + '\n' +
           'adapter_match = ' + str(c.adapter_match) + '\n' +
           'q_min = ' + str(c.q_min) + '\n' +
           'q_percent = ' + str(c.q_percent) + '\n' +
           'phred64 = ' + str(c.p64) + '\n' +
           'compress = ' + str(c.compress) + '\n')

    print('\n' + log + '\nthanks for using ngscomposer!')
    with open(os.path.join(c.proj, 'summary.txt'), 'w') as out:
        out.write(log)


if __name__ == '__main__':
    version = '0.4.9'

    parser = argparse.ArgumentParser(description=('#' * 50 + '\n' +
        ' ' * 15 + 'NGS-COMPOSER:\n' +
        '#' * 50 + '\n\n' +
        'version: ' + version + '\n\n' +
        'ngs-composer is an empirical, base-call filtering ' +
        'pipeline for preprocessing fastq libraries\n\n' +
        'basic usage: python3 composer.py -i <path_to_project_directory>\n\n' +
        'pipeline tools may be called individually from tools folder:\n\n' +
        '    crinoid.py  - qc stats\n' +
        '    scallop.py  - trimming\n' +
        '    anemone.py  - demultiplexing\n' +
        '    rotifer.py  - motif detection\n' +
        '    porifera.py - adapter removal\n' +
        '    krill.py    - quality filtering\n\n' +
        'see https://github.com/ryandkuster/ngscomposer for full usage notes\n\n' +
        ''), formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', type=str, required=True,
            help='the full or relative path to the project directory')
    parser.add_argument('-j', default=False, action="store_true", required=False,
            help='job mode (for advanced users): no project verification!')
    parser.add_argument('-v', default=False, action="store_true", required=False,
            help='verification mode: check project structure only (don\'t run)')
    args = parser.parse_args()

    c = Project()
    c.initialize(args.i)
    c.index_reader()
    c.conf_confirm()
    c.dir_test()
    r_packages()
    c.in1_ls, c.in2_ls = fastq_test(c.fastq_ls)
    demultiplex_test()

    if c.alt_dir:
        c.initialize(c.alt_dir)
        if len(os.listdir(c.proj)) != 0:
            sys.exit('alt_dir must be an empty directory')

    if args.j:
        print('job mode: skipping adapters test and size estimation')
        c.walkaway = True
        args.v = False
    else:
        adapters_test()
        dir_size(c.proj, c.fastq_ls)

    if args.v:
        sys.exit('verification complete, project looks good')

    if c.initial_qc:
        print(msg.crin_title)
        crinoid_multi(c.proj, c.fastq_ls)

    if c.front_trim > 0 and not c.bcs_index:
        print(msg.scal_title1)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = scallop_multi()

    if c.bcs_index:
        print(msg.nem_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = anemone_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.R1_bases_ls or c.R2_bases_ls:
        print(msg.rot_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = rotifer_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.end_score:
        print(msg.scal_title2)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = scallop_end_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.adapters1:
        print(msg.porf_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = porifera_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.q_min and c.q_percent:
        print(msg.kril_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = krill_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    tidy_up()

    summary_file()
