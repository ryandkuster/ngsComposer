import argparse
import gzip
import os
import sys


def rotifer_main():
    '''
    standalone, command line entry point to rotifer using stdin
    '''
    in1 = args.r1
    if args.o is None:
        proj_dir = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj_dir = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    pe_1 = proj_dir + '/pe.' + os.path.basename(in1)
    se_1 = proj_dir + '/se.' + os.path.basename(in1)
    R1_bases_ls = args.m1 if args.m1 else None
    R2_bases_ls = args.m2 if args.m2 else None
    if R1_bases_ls is None and R2_bases_ls is None:
        sys.exit('at least one motifs list must be defined')
    try:
        in2 = args.r2
        pe_2 = proj_dir + '/pe.' + os.path.basename(in2)
        se_2 = proj_dir + '/se.' + os.path.basename(in2)
        rotifer_open(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)
    except TypeError:
        rotifer_single_open(R1_bases_ls, in1, se_1)


def rotifer_comp(in1_ls, in2_ls, R1_bases_ls, R2_bases_ls, non_genomic,
            proj_dir_current, in1):
    '''
    composer entry point to rotifer
    '''
    pe_1 = proj_dir_current + '/paired/' + os.path.basename(in1)
    se_1 = proj_dir_current + '/single/' + os.path.basename(in1)
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_2 = proj_dir_current + '/paired/' + os.path.basename(in2)
        se_2 = proj_dir_current + '/single/' + os.path.basename(in2)
        rotifer_open(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)
    except (IndexError, ValueError) as e:
        rotifer_single_open(R1_bases_ls, in1, se_1)


def rotifer_open(R1_bases_ls, R2_bases_ls, in1, in2, pe_1, pe_2, se_1, se_2):
    '''
    parse single and paired-end reads for recognized motifs
    '''
    try:
        with gzip.open(in1, 'rt') as f1, gzip.open(in2, 'rt') as f2,\
                open(pe_1, "w") as pe_o1,\
                open(pe_2, "w") as pe_o2,\
                open(se_1, "w") as se_o1,\
                open(se_2, "w") as se_o2:
            rotifer(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2, se_o1, se_o2)
    except OSError:
        with open(in1) as f1, open(in2) as f2,\
                open(pe_1, "w") as pe_o1,\
                open(pe_2, "w") as pe_o2,\
                open(se_1, "w") as se_o1,\
                open(se_2, "w") as se_o2:
            rotifer(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2, se_o1, se_o2)


def rotifer(R1_bases_ls, R2_bases_ls, f1, f2, pe_o1, pe_o2, se_o1, se_o2):
    y, entry1, entry2 = 0, "", ""
    for line1, line2 in zip(f1, f2):
        y += 1
        line1 = line1.rstrip()
        line2 = line2.rstrip()
        entry1 = entry1 + line1 + "\n"
        entry2 = entry2 + line2 + "\n"
        if y == 2:
            rotifer1 = rotifer_test(line1, R1_bases_ls) if R1_bases_ls else False
            rotifer2 = rotifer_test(line2, R2_bases_ls) if R2_bases_ls else False
        if y == 4:
            if rotifer1 is False and rotifer2 is False:
                pe_o1.write(entry1)
                pe_o2.write(entry2)
            elif rotifer1 is False and rotifer2 is True:
                se_o1.write(entry1)
            elif rotifer1 is True and rotifer2 is False:
                se_o2.write(entry2)
            elif rotifer1 is True and rotifer2 is True:
                pass
            y, entry1, entry2 = 0, "", ""


def rotifer_single_open(R1_bases_ls, in1, se_1):
    '''
    parse single-end reads for recognized motifs
    '''
    try:
        with gzip.open(in1, 'rt') as f1, open(se_1, "w") as se_o1:
            rotifer_single(R1_bases_ls, f1, se_o1)
    except OSError:
        with open(in1) as f1, open(se_1, "w") as se_o1:
            rotifer_single(R1_bases_ls, f1, se_o1)


def rotifer_single():
    y, entry1 = 0, ""
    for line1 in f1:
        y += 1
        line1 = line1.rstrip()
        entry1 = entry1 + line1 + "\n"
        if y == 2:
            rotifer1 = rotifer_test(line1, R1_bases_ls)
        if y == 4:
            if rotifer1 is False:
                se_o1.write(entry1)
            elif rotifer1 is True:
                pass
            y, entry1 = 0, ""


def rotifer_test(line, bases_ls):
    '''
    test for motif at beginning of read sequence
    '''
    rotifer = True
    for seq in bases_ls:
        if line.startswith(seq):
            rotifer = False
    return rotifer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='filter fastq reads for motif')
    parser.add_argument('-r1', type=str,
            help='the full or relative path to R1 fastq file')
    parser.add_argument('-r2', type=str,
            help='the full or relative path to R2 fastq file (optional)')
    parser.add_argument('-m1', type=str, nargs='+',
            help='space separated list of motifs expected in R1')
    parser.add_argument('-m2', type=str, nargs='+',
            help='space separated list of motifs expected in R2')
    parser.add_argument('-o', type=str,
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    rotifer_main()
