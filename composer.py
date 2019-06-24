import argparse
import gzip
import os
import shutil
import subprocess
import sys
from functools import partial
from multiprocessing import Pool

from tools.anemone import anemone_comp, bc_reader, bc_test
from tools.crinoid import crinoid_comp
from tools.krill import krill_comp
from tools.rotifer import rotifer_comp
from tools.scallop import scallop_comp, scallop_end


class Project:
    def __init__(self):
        self.singles_ls = []
        self.fastq_ls = []
        self.in1_ls = []
        self.in2_ls = []
        self.rm_dirs = []

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
            sys.exit('project directory not found')

    def conf_confirm(self):
        '''
        test user-input from configuration file
        '''
        assert self.paired is True or self.paired is False
        assert self.procs >= 1
        assert os.path.exists(self.alt_dir) or self.alt_dir is False
        assert self.initial_qc is True or self.initial_qc is False
        assert self.walkthrough is True or self.walkthrough is False
        assert self.walkaway is True or self.walkaway is False
        assert isinstance(self.front_trim, int)
        if self.bcs_index:
            assert os.path.exists(os.path.join(self.proj, self.bcs_index))
        assert isinstance(self.mismatch, int)
        assert self.R1_bases_ls is False or isinstance(self.R1_bases_ls, list)
        if isinstance(self.R1_bases_ls, list):
            for i in self.R1_bases_ls:
                for j in i:
                    assert j in ['A', 'C', 'G', 'T']
        assert self.R2_bases_ls is False or isinstance(self.R2_bases_ls, list)
        if isinstance(self.R2_bases_ls, list):
            for i in self.R2_bases_ls:
                for j in i:
                    assert j in ['A', 'C', 'G', 'T']
        assert self.non_genomic is False or isinstance(self.non_genomic, int)
        if self.end_trim:
            assert 0 <= self.end_trim <= 40 and self.end_trim is not True
        if self.q_min or self.q_percent:
            assert 0 <= self.q_min <= 40 and self.q_min is not True
            assert 0 <= self.q_percent <= 100 and self.q_percent is not True
            if self.q_min and self.q_percent:
                pass
            else:
                sys.exit('both q_min and q_percent variables must be defined')
        assert self.rm_transit is True or self.rm_transit is False

    def index_reader(self):
        '''
        open bcs_index and create dictionary of associated bc keyfiles
        '''
        c.bcs_index = os.path.join(c.proj, c.bcs_index) if c.bcs_index else ''
        self.bcs_dict = {}
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
            elif os.path.join(self.proj, filename) in self.bcs_dict.values():
                pass
            elif filename == 'conf.py' or filename == '__pycache__':
                pass
            else:
                self.fastq_ls.append(os.path.join(self.proj, filename))


def r_packages():
    try:
        subprocess.check_call(['Rscript', os.path.dirname(
            os.path.abspath(sys.argv[0])) +
            '/tests/test_packages.R'], shell=False)
    except FileNotFoundError:
        sys.exit('please install latest version of R')


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
    dir_plan = dir_plan + 1 if c.end_trim else dir_plan
    dir_plan = dir_plan + 1 if c.q_min else dir_plan
    for i in fastq_ls:
        dir_used += os.path.getsize(i)
    if c.rm_transit:
        dir_plan = dir_used * 2 if fastq_test else dir_used * 2 * 5
    else:
        dir_plan = dir_used * dir_plan if fastq_test\
                else dir_used * dir_plan * 5
    if dir_plan >= drive_free:
        sys.exit('an estimated ' + str(dir_plan) + ' bytes are required to ' +
                 'process, consider rm_transit or alt_dir variables')


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
            elif root == str(curr + '/single'):
                singles_ls.append(fullname)
            else:
                fastq_ls.append(fullname)
    pairs_dict = is_paired(fastq_ls)
    in1_ls, in2_ls = input_sort(pairs_dict)
    return singles_ls, fastq_ls, in1_ls, in2_ls


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


def crinoid_multi():
    '''
    create user-defined subprocesses to produce base-call summary
    '''
    curr = os.path.join(c.proj, 'qc')
    os.mkdir(curr)
    crinoid_part = partial(crinoid_comp, curr)
    pool_multi(crinoid_part, c.fastq_ls)


def scallop_multi():
    '''
    create user-defined subprocesses to trim every file in fastq_ls
    '''
    curr = os.path.join(c.proj, 'trimmed')
    c.rm_dirs.append(curr)
    os.mkdir(curr)
    scallop_part = partial(scallop_comp, c.front_trim, None, curr)
    pool_multi(scallop_part, c.fastq_ls)
    singles_ls, fastq_ls, in1_ls, in2_ls = pathfinder(curr)
    return singles_ls, fastq_ls, in1_ls, in2_ls


def anemone_multi():
    '''
    create user-defined subprocesses to demultiplex
    '''
    curr = os.path.join(c.proj, 'demultiplexed')
    rm_dirs.append(curr)
    os.mkdir(curr)
    anemone_part = partial(anemone_comp, c.in1_ls, c.in2_ls, c.mismatch,
                           c.bcs_dict, curr)
    pool_multi(anemone_part, c.in1_ls)
    concater(curr)
    t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = pathfinder(curr)
    if c.walkthrough:
        crinoid_multi(curr, t_fastq_ls)
        update = input(msg.anemone_qc) if c.walkaway is False else 'y'
        if update not in (msg.confirm):
            prompt = input(msg.anemone_up)
            if prompt in (msg.confirm):
                c.mismatch = int(input(msg.anemone_in))
                shutil.rmtree(rm_dirs[-1])
                _, _, _, _ = anemone_multi(proj, bcs_dict, in1_ls, in2_ls)
                rm_dirs.pop()
            else:
                exit = input('\nexit ngsComposer? (y/n)\n')
                if exit in (msg.confirm):
                    sys.exit('\nngsComposer is now exiting')
    return t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls


def rotifer_multi(proj, fastq_ls, in1_ls, in2_ls):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    curr = proj + '/parsed'
    rm_dirs.append(curr)
    os.mkdir(curr)
    os.mkdir(curr + '/single')
    os.mkdir(curr + '/paired')
    rotifer_part = partial(rotifer_comp, in1_ls, in2_ls, c.R1_bases_ls,
                           c.R2_bases_ls, c.non_genomic, curr)
    pool_multi(rotifer_part, in1_ls)
    t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = pathfinder(curr)
    if c.walkthrough:
        crinoid_multi(curr + '/single', t_singles_ls)
        crinoid_multi(curr + '/paired', t_fastq_ls)
        update = input(msg.rotifer_qc) if c.walkaway is False else 'y'
        if update not in (msg.confirm):
            prompt = input(msg.rotifer_up)
            if prompt in (msg.confirm):
                shutil.rmtree(rm_dirs[-1])
                rm_dirs.pop()
                R1 = input(msg.rotifer_in1)
                c.R1_bases_ls = False if R1 == '' else R1.split()
                R2 = input(msg.rotifer_in2)
                c.R2_bases_ls = False if R2 == '' else R2.split()
                if c.R1_bases_ls is False and c.R2_bases_ls is False:
                    return singles_ls, fastq_ls, in1_ls, in2_ls
                t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = rotifer_multi(proj, fastq_ls, in1_ls, in2_ls)
            else:
                exit = input('\nexit ngsComposer? (y/n)\n')
                if exit in (msg.confirm):
                    sys.exit('\nngsComposer is now exiting')
    return t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls


def scallop_end_multi(proj, fastq_ls, singles_ls):
    '''
    automated 3' end read trimming based on minimum value
    '''
    curr = proj + '/end_trimmed'
    rm_dirs.append(curr)
    os.mkdir(curr)
    if c.R1_bases_ls or c.R2_bases_ls:
        os.mkdir(curr + '/single')
        os.mkdir(curr + '/paired')
    if os.path.exists(rm_dirs[-2] + '/qc'):
        pass
    elif os.path.exists(rm_dirs[-2] + '/single/qc'):
        pass
    else:
        try:
            crinoid_multi(rm_dirs[-2] + '/single', singles_ls)
            crinoid_multi(rm_dirs[-2] + '/paired', fastq_ls)
        except FileNotFoundError:
            crinoid_multi(rm_dirs[-2], fastq_ls)
    scallop_end_part = partial(scallop_end, curr, c.end_trim)
    pool_multi(scallop_end_part, fastq_ls)
    if singles_ls:
        scallop_end_part = partial(scallop_end, curr, c.end_trim)
        pool_multi(scallop_end_part, singles_ls)
    t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = pathfinder(curr)
    if c.walkthrough:
        if c.R1_bases_ls or c.R2_bases_ls:
            crinoid_multi(curr + '/single', t_singles_ls)
            crinoid_multi(curr + '/paired', t_fastq_ls)
        else:
            crinoid_multi(curr, t_fastq_ls)
        update = input(msg.scallop_qc) if c.walkaway is False else 'y'
        if update not in (msg.confirm):
            prompt = input(msg.scallop_up)
            if prompt in (msg.confirm):
                shutil.rmtree(rm_dirs[-1])
                rm_dirs.pop()
                c.end_trim = input(msg.scallop_in)
                c.end_trim = False if c.end_trim == '' else int(c.end_trim)
                if c.end_trim is False:
                    return singles_ls, fastq_ls, in1_ls, in2_ls
                t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = scallop_end_multi(proj, fastq_ls, singles_ls)
            else:
                exit = input('\nexit ngsComposer? (y/n)\n')
                if exit in (msg.confirm):
                    sys.exit('\nngsComposer is now exiting')
    return t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls


def krill_multi(in1_ls, in2_ls, singles_ls):
    '''
    create user-defined subprocesses to parse based on expected sequences
    '''
    curr = proj + '/filtered'
    rm_dirs.append(curr)
    os.mkdir(curr)
    os.mkdir(curr + '/single')
    os.mkdir(curr + '/single/pe_lib')
    os.mkdir(curr + '/single/se_lib')
    os.mkdir(curr + '/paired')
    krill_part = partial(krill_comp, in1_ls, in2_ls, c.q_min, c.q_percent,
                         curr)
    pool_multi(krill_part, in1_ls)
    if singles_ls:
        krill_part = partial(krill_comp, in1_ls, in2_ls, c.q_min,
                             c.q_percent, curr)
        pool_multi(krill_part, singles_ls)
    concater(curr + '/single')
    t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = pathfinder(curr)
    shutil.rmtree(curr + '/single/pe_lib')
    shutil.rmtree(curr + '/single/se_lib')
    if c.walkthrough:
        crinoid_multi(curr + '/single', t_singles_ls)
        crinoid_multi(curr + '/paired', t_fastq_ls)
        update = input(msg.krill_qc) if c.walkaway is False else 'y'
        if update not in (msg.confirm):
            prompt = input(msg.krill_up)
            if prompt in (msg.confirm):
                shutil.rmtree(rm_dirs[-1])
                rm_dirs.pop()
                c.q_min = input(msg.krill_in1)
                c.q_percent = input(msg.krill_in2)
                c.q_min = False if c.q_min == '' else int(c.q_min)
                c.q_percent = False if c.q_percent == '' else int(c.q_percent)
                if c.q_min is False or c.q_percent is False:
                    return singles_ls, fastq_ls, in1_ls, in2_ls
                t_singles_ls, t_fastq_ls, t_in1_ls, t_in2_ls = krill_multi(in1_ls, in2_ls, singles_ls)
            else:
                exit = input('\nexit ngsComposer? (y/n)\n')
                if exit in (msg.confirm):
                    sys.exit('\nngsComposer is now exiting')
    return singles_ls, fastq_ls, in1_ls, in2_ls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='composer is a base-call\
        error-filtering and read preprocessing pipeline for fastq libraries - \
        see https://github.com/ryandkuster/composer for usage')
    parser.add_argument('-i', type=str, help='the full or relative path to \
        the project directory')
    args = parser.parse_args()
    c = Project()
    c.initialize(args.i)
    c.conf_confirm()

######################################################
# TODO delete the following (for ease of testing):
######################################################
    old_dirs = [c.proj + '/qc',
                c.proj + '/trimmed',
                c.proj + '/demultiplexed',
                c.proj + '/parsed',
                c.proj + '/end_trimmed',
                c.proj + '/filtered']
    for folder in old_dirs:
        try:
            shutil.rmtree(folder)
            dir_name = os.path.basename(folder)
            print('\ncomposer is removing the ' + dir_name + ' directory')
        except FileNotFoundError:
            pass
######################################################

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

    for filename in c.fastq_ls:
        try:
            fastq_test = is_fq(filename)
        except UnicodeDecodeError:
            fastq_test = is_gz(filename)
        if fastq_test is None:
            sys.exit(filename + ' was not expected in project directory')

    if c.alt_dir:
        c.proj = initialize(c.alt_dir)
        if len(os.listdir(c.proj)) != 0:
            sys.exit('alt_dir must be an empty directory')

    dir_size(c.proj, c.fastq_ls, fastq_test)

    if c.bcs_index and c.paired is True and \
            len(fastq_ls)/2 != len(bcs_dict):
        sys.exit('incorrect number of files based on index.txt')
    pairs_dict = is_paired(c.fastq_ls)
    in1_ls, in2_ls = input_sort(pairs_dict)

    '''
    begin calling tools
    '''
    import tools.helpers.messages as msg
    if c.initial_qc is True:
        print(msg.crin_title)
        crinoid_multi()

    if c.front_trim > 0:
        print(msg.scal_title1)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = scallop_multi()

    if c.bcs_index:
        print(msg.nem_title)
        c.singles_ls, c.fastq_ls, c.in1_ls, c.in2_ls = anemone_multi()
        if c.rm_transit is True:
            dir_del(rm_dirs[:-1])

    if c.R1_bases_ls or c.R2_bases_ls:
        print(msg.rot_title)
        singles_ls, fastq_ls, in1_ls, in2_ls = rotifer_multi(proj, fastq_ls,
                in1_ls, in2_ls)
        if c.rm_transit is True:
            dir_del(rm_dirs[:-1])

    if c.end_trim:
        print(msg.scal_title2)
        singles_ls, fastq_ls, in1_ls, in2_ls = scallop_end_multi(proj,
                fastq_ls, singles_ls)
        if c.rm_transit is True:
            dir_del(rm_dirs[:-1])

    if c.q_min and c.q_percent:
        print(msg.kril_title)
        singles_ls, fastq_ls, in1_ls, in2_ls = krill_multi(in1_ls, in2_ls,
                singles_ls)
        if c.rm_transit is True:
            dir_del(rm_dirs[:-1])

    print('\n',
          'paired =', c.paired, '\n',
          'procs =', c.procs, '\n',
          'alt_dir =', c.alt_dir, '\n',
          'initial_qc =', c.initial_qc, '\n',
          'walkthrough =', c.walkthrough, '\n',
          'walkaway =', c.walkaway, '\n',
          'front_trim =', c.front_trim, '\n',
          'bcs_index =', c.bcs_index, '\n',
          'mismatch =', c.mismatch, '\n',
          'R1_bases_ls =', c.R1_bases_ls, '\n',
          'R2_bases_ls =', c.R2_bases_ls, '\n',
          'non_genomic =', c.non_genomic, '\n',
          'q_min =', c.q_min, '\n',
          'q_percent =', c.q_percent, '\n',
          'end_trim =', c.end_trim, '\n',
          'rm_transit =', c.rm_transit,  '\n')





