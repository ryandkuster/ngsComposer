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
        self.bypass = False

        self.paired = False
        self.procs = 1
        self.alt_dir = False
        self.initial_qc = False
        self.all_qc = False
        self.walkaway = True
        self.front_trim = False
        self.end_trim = False
        self.bcs_index = False
        self.mismatch = False
        self.R1_bases_ls = False
        self.R2_bases_ls = False
        self.non_genomic = False
        self.q_min = False
        self.q_percent = False
        self.auto_trim = False
        self.rm_transit = False

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
        if self.alt_dir:
            assert os.path.exists(self.alt_dir)
        if self.initial_qc:
            assert self.initial_qc is True
        if self.all_qc:
            assert self.all_qc is True
        if self.walkaway:
            assert self.walkaway is True
        else:
            self.all_qc = True
        if self.front_trim:
            assert isinstance(self.front_trim, int)
            assert self.front_trim > 0
        if self.end_trim:
            assert isinstance(self.end_trim, int)
            assert self.end_trim > 0
        if self.end_trim and self.auto_trim:
            print(msg.conf_end)
            raise AssertionError
        if self.bcs_index or self.mismatch:
            assert os.path.exists(os.path.join(self.proj, self.bcs_index))
            assert isinstance(self.mismatch, int)
        if self.R1_bases_ls:
            assert isinstance(self.R1_bases_ls, list)
            nucleotide_test(self.R1_bases_ls)
        if self.R2_bases_ls:
            assert isinstance(self.R2_bases_ls, list)
            nucleotide_test(self.R2_bases_ls)
        if self.non_genomic:
            assert isinstance(self.non_genomic, int)
        if self.auto_trim:
            assert isinstance(self.auto_trim, int)
            assert 0 <= self.auto_trim <= 40 and self.auto_trim is not True
        if self.q_min or self.q_percent:
            assert isinstance(self.q_min, int)
            assert isinstance(self.q_percent, int)
            assert 0 <= self.q_min <= 40 and self.q_min is not True
            assert 0 <= self.q_percent <= 100 and self.q_percent is not True
            if self.q_min and self.q_percent:
                pass
            else:
                sys.exit(msg.q_vars)
        assert self.rm_transit is True or self.rm_transit is False

    def index_reader(self):
        '''
        open bcs_index and create dictionary of associated bc keyfiles
        '''
        c.bcs_index = os.path.join(c.proj, c.bcs_index) if c.bcs_index else ''
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
            elif os.path.join(self.proj, filename) in self.bcs_dict.values():
                pass
            elif filename == 'conf.py' or filename == '__pycache__':
                pass
            else:
                self.fastq_ls.append(os.path.join(self.proj, filename))


def nucleotide_test(ls):
    nuc_ls = (j for i in ls for j in i)
    for i in nuc_ls:
        if i in ['A', 'C', 'G', 'T']:
            test = True
            pass
        elif i.upper() in ['A', 'C', 'G', 'T']:
            test = msg.nucs1
            break
        else:
            test = msg.nucs2
            break
    if test is not True:
        print(test)
        if not c.bypass:
            sys.exit()
        raise AssertionError


def r_packages():
    try:
        subprocess.check_call(['Rscript', os.path.join(os.path.dirname(
            __file__), 'tests', 'test_packages.R')], shell=False)
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
    dir_plan = dir_plan + 1 if c.auto_trim else dir_plan
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


def dir_make(title):
    '''
    create new directory when pipeline tools called
    '''
    # TODO add following at release:
    # curr = os.path.join(c.proj, str(len(c.rm_dirs)) + '_' + title)
    curr = os.path.join(c.proj, title)
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


def trim_assist():
    # TODO make function for the complex decisions with end-trimming:
    # collapse qc data and plot
    # change autotrim between q1, mean, median, lw
    pass


def combine_matrix():
    with open(args.a1) as f1, open(args.a2) as f2:
        a = [line.rstrip().split(',') for line in f1]
        b = [line.rstrip().split(',') for line in f2]
        if b:
            pass
        else:
            b = [[0 for j in range(len(a[0]))] for i in range(args.l)]
        for i, col in enumerate(b):
            for j, score in enumerate(col):
                try:
                    b[i][j] = int(b[i][j]) + int(a[i][j])
                except IndexError:
                    break
    with open(args.o, 'w') as o1:
        for row in b:
            o1.write(",".join(str(x) for x in row))
            o1.write("\n")


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

    if c.walkaway:
        return temp_ls

    while True:
        status = msg.walk2 if not c.walkaway else msg.walk3
        print('\nwalkthrough mode is ' + status)
        print('\nplease check the qc folders in ' + c.rm_dirs[-1])
        choice1 = input(msg.walk1)
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

    shutil.rmtree(c.rm_dirs[-1])
    c.rm_dirs.pop()

    while True:
        for k, v in kwargs.items():
            print('\ncurrent value for ' + k + ' is ' + str(v))
            v = input('new value for ' + k + '? > ')
            if v.startswith('['):
                v = v.strip('[\']').split(', ')
            elif v == 'False':
                v = False
            elif any(i.isdigit() for i in v):
                v = int(v)
            setattr(c, k, v)
        try:
            c.conf_confirm()
            break
        except AssertionError:
            print('\nnot a valid entry\n')

    print('\nrerunning step with updated variables...')
    temp_ls = tool()
    return temp_ls


def crinoid_multi(proj, ls):
    '''
    create user-defined subprocesses to produce base-call summary
    '''
    curr = os.path.join(proj, 'qc')
    os.mkdir(curr)
    crinoid_part = partial(crinoid_comp, curr)
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
    if c.R1_bases_ls or c.R2_bases_ls:
        os.mkdir(os.path.join(curr, 'single'))
        os.mkdir(os.path.join(curr, 'paired'))

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

    scallop_end_part = partial(scallop_end, curr, c.auto_trim)
    pool_multi(scallop_end_part, c.fastq_ls)
    if c.singles_ls:
        scallop_end_part = partial(scallop_end, curr, c.auto_trim)
        pool_multi(scallop_end_part, c.singles_ls)
    temp_ls = pathfinder(curr)
    if c.all_qc:
        temp_ls = walkthrough(curr, scallop_end_multi, temp_ls,
                              auto_trim=c.auto_trim)
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
                         curr)
    pool_multi(krill_part, c.in1_ls)
    if c.singles_ls:
        krill_part = partial(krill_comp, c.in1_ls, c.in2_ls, c.q_min,
                             c.q_percent, curr)
        pool_multi(krill_part, c.singles_ls)
    concater(os.path.join(curr, 'single'))
    temp_ls = pathfinder(curr)
    shutil.rmtree(os.path.join(curr, 'single', 'pe_lib'))
    shutil.rmtree(os.path.join(curr, 'single', 'se_lib'))
    if c.all_qc:
        temp_ls = walkthrough(curr, krill_multi, temp_ls, q_min=c.q_min,
                              q_percent=c.q_percent)
    return temp_ls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='composer is a base-call\
        error-filtering and read preprocessing pipeline for fastq libraries - \
        see https://github.com/ryandkuster/composer for usage')
    parser.add_argument('-i', type=str, help='the full or relative path to \
        the project directory')
    args = parser.parse_args()
    c = Project()
    c.initialize(args.i)
    import tools.helpers.messages as msg
    c.conf_confirm()

######################################################
# TODO delete the following (for ease of testing):
######################################################
    old_dirs = [os.path.join(c.proj, 'qc'),
                os.path.join(c.proj, 'trimmed'),
                os.path.join(c.proj, 'demultiplexed'),
                os.path.join(c.proj, 'parsed'),
                os.path.join(c.proj, 'auto_trimmed'),
                os.path.join(c.proj, 'filtered')]
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

    if c.front_trim > 0 or c.end_trim > 0:
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

    print('\n',
          'paired =', c.paired, '\n',
          'procs =', c.procs, '\n',
          'alt_dir =', c.alt_dir, '\n',
          'initial_qc =', c.initial_qc, '\n',
          'all_qc =', c.all_qc, '\n',
          'walkaway =', c.walkaway, '\n',
          'front_trim =', c.front_trim, '\n',
          'end_trim =', c.end_trim, '\n',
          'bcs_index =', c.bcs_index, '\n',
          'mismatch =', c.mismatch, '\n',
          'R1_bases_ls =', c.R1_bases_ls, '\n',
          'R2_bases_ls =', c.R2_bases_ls, '\n',
          'non_genomic =', c.non_genomic, '\n',
          'q_min =', c.q_min, '\n',
          'q_percent =', c.q_percent, '\n',
          'auto_trim =', c.auto_trim, '\n',
          'rm_transit =', c.rm_transit,  '\n')
