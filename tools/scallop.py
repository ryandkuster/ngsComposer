import argparse
import gzip
import os
import sys

from helpers.gzip_handling import gzip_test


def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    in1 = args.r1
    in2 = args.r2 if args.r2 else None
    front_trim = args.f
    end_trim = None if args.b == 0 else args.b
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    pe_1 = os.path.join(proj, 'trimmed_pe.' + os.path.basename(in1))
    se_1 = os.path.join(proj, 'trimmed_se.' + os.path.basename(in1))
    if args.e or args.w or args.l:
        if not args.e or not args.w or not args.l:
            sys.exit('-e, -w, and -l must be defined to end trim')
        end_score = args.e
        window = args.w
        min_l = args.l
        if front_trim:
            min_l = front_trim + min_l
        if in2:
            pe_2 = os.path.join(proj, 'trimmed_pe.' + os.path.basename(in2))
            se_2 = os.path.join(proj, 'trimmed_se.' + os.path.basename(in2))
            scallop_end_open(in1, in2, pe_1, pe_2, se_1, se_2, front_trim, end_score, window, min_l)
        else:
            scallop_single_end_open(in1, se_1, front_trim, end_score, window, min_l) 
    else:
        scallop_open(in1, front_trim, end_trim, se_1)


def scallop_comp(in1_ls, in2_ls, front_trim, end_trim, end_score, window,
                 min_l, curr, in1):
    '''
    composer entry point to scallop
    '''
    if end_score:
        try:
            in2 = in2_ls[in1_ls.index(in1)]
            pe_1 = os.path.join(curr, 'paired', os.path.basename(in1))
            se_1 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in1))
            pe_2 = os.path.join(curr, 'paired', os.path.basename(in2))
            se_2 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in2))
        except (IndexError, ValueError) as e:
            se_1 = os.path.join(curr, 'single', 'se_lib', os.path.basename(in1))
            in2 = None
        if front_trim:
            min_l = front_trim + min_l
        if in2:
            scallop_end_open(in1, in2, pe_1, pe_2, se_1, se_2, front_trim,
                            end_score, window, min_l)
        else:
            scallop_single_end_open(in1, se_1, front_trim, end_score, window, min_l) 
    else:
        se_1 = os.path.join(curr, os.path.basename(in1))
        end_trim = None if end_trim == 0 else end_trim
        scallop_open(in1, front_trim, end_trim, se_1)


def scallop_open(in1, front_trim, end_trim, out1):
    '''
    test if gzipped, then open single-end files
    '''
    compressed = gzip_test(in1)
    if compressed:
        out1, _ = os.path.splitext(out1)
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            scallop(front_trim, end_trim, f, o)
    else:
        with open(in1) as f, open(out1, 'w') as o:
            scallop(front_trim, end_trim, f, o)


def scallop(front_trim, end_trim, f, o):
    i = 0
    for line in f:
        i += 1
        if i % 2 == 0:
            o.write(line.rstrip()[front_trim:end_trim] + "\n")
        else:
            o.write(line)


def scallop_end_open(in1, in2, pe_1, pe_2, se_1, se_2, front_trim, end_score, window, min_l):
    '''
    test if gzipped, then open paired-end files
    '''
    compressed = gzip_test(in1)
    if compressed:
        pe_1, _ = os.path.splitext(pe_1)
        pe_2, _ = os.path.splitext(pe_2)
        se_1, _ = os.path.splitext(se_1)
        se_2, _ = os.path.splitext(se_2)
        with gzip.open(in1, 'rt') as f1, gzip.open(in2, 'rt') as f2,\
                open(pe_1, "w") as pe_o1,\
                open(pe_2, "w") as pe_o2,\
                open(se_1, "w") as se_o1,\
                open(se_2, "w") as se_o2:
            scallop_end_line(f1, f2, pe_o1, pe_o2, se_o1, se_o2, front_trim, end_score, window, min_l)
    else:
        with open(in1) as f1, open(in2) as f2,\
                open(pe_1, "w") as pe_o1,\
                open(pe_2, "w") as pe_o2,\
                open(se_1, "w") as se_o1,\
                open(se_2, "w") as se_o2:
            scallop_end_line(f1, f2, pe_o1, pe_o2, se_o1, se_o2, front_trim, end_score, window, min_l)


def scallop_end_line(f1, f2, pe_o1, pe_o2, se_o1, se_o2, front_trim, end_score, window, min_l):
    '''
    automated end-trimming for paired-end reads
    '''
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    val = dict(zip(scores[:43], range(0, 43)))
    good_ls = [k for k, v in val.items() if int(v) >= int(end_score)]
    i = 0
    entry1, entry2 = [], []
    for line1, line2 in zip(f1, f2):
        entry1.append(line1.rstrip())
        entry2.append(line2.rstrip())
        i += 1
        if i == 4:
            end_trim1 = viewfinder(line1.rstrip(), good_ls, window)
            end_trim2 = viewfinder(line2.rstrip(), good_ls, window)
            if end_trim1 >= min_l and end_trim2 >= min_l:
                scallop_writer(entry1, front_trim, end_trim1, pe_o1)
                scallop_writer(entry2, front_trim, end_trim2, pe_o2)
            elif end_trim1 >= min_l:
                scallop_writer(entry1, front_trim, end_trim1, se_o1)
            elif end_trim2 >= min_l:
                scallop_writer(entry2, front_trim, end_trim2, se_o2)
            entry1, entry2 = [], []
            i = 0


def scallop_single_end_open(in1, out1, front_trim, end_score, window, min_l):
    '''
    test if gzipped, then open single-end files
    '''
    compressed = gzip_test(in1)
    if compressed:
        out1, _ = os.path.splitext(out1)
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            scallop_single_end_line(f, o, front_trim, end_score, window, min_l)
    else:
        with open(in1) as f, open(out1, 'w') as o:
            scallop_single_end_line(f, o, front_trim, end_score, window, min_l)


def scallop_single_end_line(f, o, front_trim, end_score, window, min_l):
    '''
    automated end-trimming for single-end reads
    '''
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    val = dict(zip(scores[:43], range(0, 43)))
    good_ls = [k for k, v in val.items() if int(v) >= int(end_score)]
    i = 0
    entry = []
    for line in f:
        entry.append(line.rstrip())
        i += 1
        if i == 4:
            end_trim = viewfinder(line.rstrip(), good_ls, window)
            if end_trim >= min_l:
                scallop_writer(entry, front_trim, end_trim, o)
            entry = []
            i = 0


def viewfinder(line, good_ls, window):
    '''
    search defined window for passing scores
    '''
    for i in range(len(line) - window, -1, -1):
        for pos, j in enumerate(line[i:i+window]):
            if j not in good_ls:
                break
            elif pos == window - 1:
                return i+window
    return 0


def scallop_writer(entry, front_trim, end_trim, o):
    '''
    write passing entries to file
    '''
    entry[1] = entry[1][front_trim:end_trim]
    entry[3] = entry[3][front_trim:end_trim]
    for element in entry:
        o.write('%s\n' % element)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='trim ends of fastq reads')
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-r2', type=str, metavar='', 
            help='the full or relative path to R2 fastq file (optional)')
    parser.add_argument('-f', type=int, metavar='',
            help='number of bases to remove from beginning of read (integer)')
    parser.add_argument('-b', type=int, metavar='',
            help='final position to keep within a read (integer)')
    parser.add_argument('-e', type=int, metavar='',
            help='end-trim where entire window >= score (integer)')
    parser.add_argument('-w', type=int, metavar='',
            help='use with \'e\', size of window consisting of >= e (integer)')
    parser.add_argument('-l', type=int, metavar='',
            help='use with \'e\' & \'w\', minimum read length to keep (integer)')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    scallop_main()
