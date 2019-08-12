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

import tools.helpers.messages as msg
from tools.anemone import anemone_comp, bc_reader, bc_test
from tools.crinoid import combine_matrix, crinoid_comp
from tools.krill import krill_comp
from tools.porifera import porifera_comp
from tools.rotifer import rotifer_comp
from tools.scallop import scallop_comp, scallop_end


class Project:
    def __init__(self):
        self.singles_ls = []
        self.fastq_ls = []
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
        self.front_trim = False
        self.bcs_index = ''
        self.mismatch = False
        self.R1_bases_ls = False
        self.R2_bases_ls = False
        self.non_genomic = False
        self.trim_mode = False
        self.auto_trim = False
        self.q_min = False
        self.q_percent = False
        self.adapters = ''
        self.min_start = False
        self.adapter_mismatch = False
        self.rm_transit = True
        self.p64 = False

    def initialize(self, proj):
        '''
        check files in project directory for completeness
        '''
        if os.path.exists(proj) is True:
            self.proj = os.path.abspath(proj)
            sys.path.append(self.proj)
            import conf as config
            c.__dict__.update(config.__dict__)
        else:
            sys.exit(msg.initialize1)

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
        if self.auto_trim or self.trim_mode:
            if not self.auto_trim or not self.trim_mode:
                raise Exception(msg.conf_confirm15)
            if not isinstance(self.auto_trim, int) or self.auto_trim is True:
                raise Exception(msg.conf_confirm16)
            if not 0 <= self.auto_trim <= 42:
                raise Exception(msg.conf_confirm16)
            if self.trim_mode not in ['whisker', 'quartile', 'median', 'mean']:
                raise Exception(msg.conf_confirm17)
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
        if self.min_start:
            if not isinstance(self.min_start, int):
                raise Exception(msg.conf_confirm22)
            if not self.adapter_mismatch:
                raise Exception(msg.conf_confirm23)
            if not os.path.exists(os.path.join(self.proj, 'adapters.txt')):
                raise Exception(msg.conf_confirm23)
        if self.adapter_mismatch:
            if not isinstance(self.adapter_mismatch, int):
                raise Exception(msg.conf_confirm24)
            if not self.min_start:
                raise Exception(msg.conf_confirm23)
            if not os.path.exists(os.path.join(self.proj, 'adapters.txt')):
                raise Exception(msg.conf_confirm23)
        if not isinstance(self.p64, bool):
            raise Exception(msg.conf_confirm25)

    def index_reader(self):
        '''
        open bcs_index and create dictionary of associated bc keyfiles
        '''
        adapters = os.path.join(self.proj, 'adapters.txt')
        self.adapters = adapters if os.path.exists(adapters) else ''
        if self.adapters:
            with open(self.adapters) as f:
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

    def dir_test(self):
        '''
        test project directory for correct structure
        '''
        self.fastq_ls = []
        for filename in os.listdir(self.proj):
            if filename == os.path.basename(self.bcs_index):
                pass
            elif filename == os.path.basename(self.adapters):
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
            __file__), 'tests', 'test_packages.R')], shell=False)
    except FileNotFoundError:
        sys.exit(msg.r_packages1)


def nucleotide_test(ls):
    '''
    test list of sequences for proper formatting
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


def is_fq(filename):
    '''
    test first read structure if fastq
    '''
    try:
        with open(filename) as f:
            for i, line in enumerate(f):
                if i == 0 and line[0] != '@':
                    return
                if i == 2 and line[0] != '+':
                    return
                else:
                    return 1
    except IsADirectoryError:
        raise Exception(msg.is_fq1)


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
    '''
    return gz status if gzipped
    '''
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
    if headers identical, order in1 and in2 lists to keep pairing
    '''
    in1_ls, in2_ls = [], []
    for values in pairs_dict.values():
        if c.paired is True and len(values) == 2:
            in1_ls, in2_ls = input_paired(values, in1_ls, in2_ls)
        elif c.paired is True and len(values) != 2:
            sys.exit('R1 and R2 pairs not found, expect paired library\n' +
                     'check the naming conventions of the bcs_index')
        if c.paired is False:
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
    if len(values) == 1:
        for filename in values:
            in1_ls.append(filename)
    elif c.ignore is True:
        for filename in values:
            in1_ls.append(filename)
    else:
        while True:
            ans = input(msg.input_single1)
            if ans == '1':
                c.ignore = True
                break
            if ans == '2':
                c.paired = True
                break
            elif ans == '3':
                sys.exit('\nngs-composer is now exiting')
            else:
                print('\nselect an option from the list\n')
        for filename in values:
            in1_ls.append(filename)
    return in1_ls, in2_ls


def dir_size(proj, fastq_ls, fastq_test):
    '''
    test the specified output directory for adequate disk space
    '''
    drive_stats = os.statvfs(proj)
    drive_free = drive_stats.f_frsize * drive_stats.f_bavail
    dir_used = 0
    dir_plan = 0
    dir_plan = dir_plan + 1 if c.front_trim else dir_plan
    dir_plan = dir_plan + 1 if c.bcs_index else dir_plan
    if c.R1_bases_ls:
        dir_plan += 1
    elif c.R2_bases_ls:
        dir_plan += 1
    dir_plan = dir_plan + 1 if c.auto_trim else dir_plan
    dir_plan = dir_plan + 1 if c.q_min else dir_plan
    dir_plan = dir_plan + 1 if c.adapters else dir_plan
    for i in fastq_ls:
        dir_used += os.path.getsize(i)
    if c.rm_transit:
        dir_plan = dir_used * 2 if fastq_test else dir_used * 2 * 5
    else:
        dir_plan = dir_used * dir_plan if fastq_test\
                else dir_used * dir_plan * 5
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
    walk current directory and pull files as lists
    '''
    singles_ls, fastq_ls, in1_ls, in2_ls = [], [], [], []
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
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(pairs_dict)
    return [singles_ls, fastq_ls, in1_ls, in2_ls]


def concater(curr):
    '''
    walk current directory and concatenate files with identical names
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
    crinoid_part = partial(crinoid_comp, curr, all_qc, c.p64)
    pool_multi(crinoid_part, ls)


def scallop_multi():
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    '''
    curr = dir_make('trimmed')
    scallop_part = partial(scallop_comp, c.front_trim, None, curr)
    pool_multi(scallop_part, c.fastq_ls)
    temp_ls = pathfinder(curr)
    return temp_ls


def anemone_multi():
    '''
    create user-defined subprocesses to demultiplex
    '''
    curr = dir_make('demultiplexed')
    anemone_part = partial(anemone_comp, c.in1_ls, c.in2_ls, c.mismatch,
                           c.bcs_dict, curr)
    pool_multi(anemone_part, c.in1_ls)
    concater(curr)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, anemone_multi, temp_ls,
                              mismatch=c.mismatch)
    return temp_ls


def rotifer_multi():
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    curr = dir_make('parsed')
    os.mkdir(os.path.join(curr, 'single'))
    os.mkdir(os.path.join(curr, 'paired'))
    rotifer_part = partial(rotifer_comp, c.in1_ls, c.in2_ls, c.R1_bases_ls,
                           c.R2_bases_ls, c.non_genomic, curr)
    pool_multi(rotifer_part, c.in1_ls)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, rotifer_multi, temp_ls,
                              R1_bases_ls=c.R1_bases_ls,
                              R2_bases_ls=c.R2_bases_ls,
                              non_genomic=c.non_genomic)
    return temp_ls


def scallop_end_multi():
    '''
    automated 3' end read trimming based on minimum value
    '''
    curr = dir_make('auto_trimmed')

    if len(c.rm_dirs) >= 2:
        if os.path.exists(os.path.join(c.rm_dirs[-2], 'qc')):
            pass
        elif os.path.exists(os.path.join(c.rm_dirs[-2], 'single', 'qc')):
            pass
        else:
            try:
                crinoid_multi(os.path.join(c.rm_dirs[-2], 'single'),
                              c.singles_ls)
                crinoid_multi(os.path.join(c.rm_dirs[-2], 'paired'),
                              c.fastq_ls)
            except FileNotFoundError:
                crinoid_multi(c.rm_dirs[-2], c.fastq_ls)
    else:
        if os.path.exists(os.path.join(c.proj, 'qc')):
            pass
        else:
            crinoid_multi(c.proj, c.fastq_ls)

    if os.path.exists(os.path.join(c.rm_dirs[-2], 'single')):
        os.mkdir(os.path.join(curr, 'single'))
        os.mkdir(os.path.join(curr, 'paired'))

    scallop_end_part = partial(scallop_end, curr, c.auto_trim, c.trim_mode)
    pool_multi(scallop_end_part, c.fastq_ls)

    if c.singles_ls:
        scallop_end_part = partial(scallop_end, curr, c.auto_trim, c.trim_mode)
        pool_multi(scallop_end_part, c.singles_ls)
    temp_ls = pathfinder(curr)

    if c.all_qc:
        temp_ls = walkthrough(curr, scallop_end_multi, temp_ls,
                              auto_trim=c.auto_trim,
                              trim_mode=c.trim_mode)
    return temp_ls


def krill_multi():
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    curr = dir_make('filtered')
    os.mkdir(os.path.join(curr, 'single'))
    os.mkdir(os.path.join(curr, 'single', 'pe_lib'))
    os.mkdir(os.path.join(curr, 'single', 'se_lib'))
    os.mkdir(os.path.join(curr, 'paired'))
    krill_part = partial(krill_comp, c.in1_ls, c.in2_ls, c.q_min, c.q_percent,
                         c.p64, curr)
    pool_multi(krill_part, c.in1_ls)
    if c.singles_ls:
        krill_part = partial(krill_comp, c.in1_ls, c.in2_ls, c.q_min,
                             c.q_percent, c.p64, curr)
        pool_multi(krill_part, c.singles_ls)
    concater(os.path.join(curr, 'single'))
    temp_ls = pathfinder(curr)
    shutil.rmtree(os.path.join(curr, 'single', 'pe_lib'))
    shutil.rmtree(os.path.join(curr, 'single', 'se_lib'))
    if c.all_qc:
        temp_ls = walkthrough(curr, krill_multi, temp_ls,
                              q_min=c.q_min,
                              q_percent=c.q_percent)
    return temp_ls


def porifera_multi():
    '''
    create user-defined subprocesses to detect and remove adapter sequences
    '''
    curr = dir_make('adapted')

    if len(c.rm_dirs) >= 2:
        if os.path.exists(os.path.join(c.rm_dirs[-2], 'single')):
            os.mkdir(os.path.join(curr, 'single'))
            os.mkdir(os.path.join(curr, 'paired'))

    porifera_part = partial(porifera_comp, curr, c.adapters, c.min_start,
                            c.adapter_mismatch)
    pool_multi(porifera_part, c.fastq_ls)
    if c.singles_ls:
        pool_multi(porifera_part, c.singles_ls)
    temp_ls = pathfinder(curr)
    return temp_ls


def walkthrough(curr, tool, temp_ls, **kwargs):
    '''
    query user for modifying or accepting current step in pipeline
    '''
    c.bypass = True
    try:
        crinoid_multi(os.path.join(curr, 'single'), temp_ls[0])
        crinoid_multi(os.path.join(curr, 'paired'), temp_ls[1])
    except FileNotFoundError:
        crinoid_multi(curr, temp_ls[1])

    if len(temp_ls[2]) >= 1:
        combine_matrix(temp_ls[2], 'R1_summary.txt')
    if len(temp_ls[3]) >= 1:
        combine_matrix(temp_ls[3], 'R2_summary.txt')
    if len(temp_ls[0]) >= 1:
        combine_matrix(temp_ls[0], 'singles_summary.txt')

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('#' * 50 + '\n' +
        ' ' * 15 + 'NGS-COMPOSER:\n' +
        '#' * 50 + '\n' +
        'ngs-composer is an empirical, base-call filtering ' +
        'pipeline for preprocessing fastq libraries\n\n' +
        'basic usage: python3 composer.py -i <path_to_project_directory>\n\n' +
        'pipeline tools may be called individually from tools folder:\n\n' +
        '    crinoid.py  - qc stats\n' +
        '    scallop.py  - trimming\n' +
        '    anemone.py  - demultiplexing\n' +
        '    rotifer.py  - motif detection\n' +
        '    krill.py    - quality filtering\n' +
        '    porifera.py - adapter removal\n\n' +
        'see https://github.com/ryandkuster/ngs-composer for full usage notes\n\n' +
        ''), formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', type=str, help='the full or relative path to ' +
        'the project directory')
    args = parser.parse_args()
    c = Project()
    c.initialize(args.i)
    c.conf_confirm()
    c.index_reader()
    c.dir_test()
    r_packages()

    for k, v in c.bcs_dict.items():
        try:
            assert os.path.join(c.proj, k) in c.fastq_ls
        except AssertionError:
            sys.exit('check index.txt formatting')
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

    for filename in c.fastq_ls:
        try:
            fastq_test = is_fq(filename)
        except UnicodeDecodeError:
            fastq_test = is_gz(filename)
        if fastq_test is None:
            sys.exit(filename + ' was not expected in project directory')

    if c.alt_dir:
        c.initialize(c.alt_dir)
        if len(os.listdir(c.proj)) != 0:
            sys.exit('alt_dir must be an empty directory')

    dir_size(c.proj, c.fastq_ls, fastq_test)

    if c.bcs_index and c.paired is True and \
            len(c.fastq_ls)/2 != len(c.bcs_dict):
        sys.exit('incorrect number of files based on index.txt')
    c.pairs_dict = is_paired(c.fastq_ls)
    c.in1_ls, c.in2_ls = input_sort(c.pairs_dict)

    if c.initial_qc:
        print(msg.crin_title)
        crinoid_multi(c.proj, c.fastq_ls)

    if c.front_trim > 0:
        print(msg.scal_title1)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = scallop_multi()

    if c.bcs_index:
        print(msg.nem_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = anemone_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.R1_bases_ls or c.R2_bases_ls or c.non_genomic:
        print(msg.rot_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = rotifer_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.auto_trim:
        print(msg.scal_title2)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = scallop_end_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.q_min and c.q_percent:
        print(msg.kril_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = krill_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])

    if c.adapters:
        print(msg.porf_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = porifera_multi()
        if c.rm_transit is True:
            dir_del(c.rm_dirs[:-1])



    log = (str(datetime.datetime.now()).split('.')[0] + '\n\n' +
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
           'auto_trim = ' + str(c.auto_trim) + '\n' +
           'trim_mode = ' + str(c.trim_mode) + '\n' +
           'q_min = ' + str(c.q_min) + '\n' +
           'q_percent = ' + str(c.q_percent) + '\n' +
           'adapters = ' + str(c.adapters) + '\n' +
           'min_start = ' + str(c.min_start) + '\n' +
           'adapter_mismatch = ' + str(c.adapter_mismatch) + '\n' +
           'phred64 = ' + str(c.p64) + '\n')

    print(log)
    with open(os.path.join(c.proj, 'summary.txt'), 'w') as out:
        out.write(log)
