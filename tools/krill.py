import argparse
import gzip
import os
import sys

from helpers.gzip_handling import gzip_test


def krill_main():
    '''
    standalone, command line entry point to krill using stdin
    '''
    q_min = args.q
    q_percent = args.p
    in1 = args.r1
    in2 = args.r2 if args.r2 else None
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    p64 = args.s if args.s else False
    pe_1 = os.path.join(proj, 'pe.' + os.path.basename(in1))
    se_1 = os.path.join(proj, 'se.' + os.path.basename(in1))
    if in2:
        pe_2 = os.path.join(proj, 'pe.' + os.path.basename(in2))
        se_2 = os.path.join(proj, 'se.' + os.path.basename(in2))
        krill_open(q_min, q_percent, p64, in1, in2, pe_1, pe_2, se_1, se_2)
    else:
        krill_single_open(q_min, q_percent, p64, in1, se_1)


def krill_comp(in1_ls, in2_ls, q_min, q_percent, p64, curr, in1):
    '''
    composer entry point to krill
    '''
    pe_1 = os.path.join(curr, 'paired', os.path.basename(in1))
    se_1 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in1))
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_2 = os.path.join(curr, 'paired', os.path.basename(in2))
        se_2 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in2))
        krill_open(q_min, q_percent, p64, in1, in2, pe_1, pe_2, se_1, se_2)
    except (IndexError, ValueError) as e:
        se_1 = os.path.join(curr, 'single', 'se_lib', os.path.basename(in1))
        krill_single_open(q_min, q_percent, p64, in1, se_1)


def krill_open(q_min, q_percent, p64, in1, in2, pe_1, pe_2, se_1, se_2):
    '''
    quality filter test for single and paired-end reads
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
            krill(q_min, q_percent, p64, f1, f2, pe_o1, pe_o2, se_o1, se_o2)
    else:
        with open(in1) as f1, open(in2) as f2,\
                open(pe_1, "w") as pe_o1,\
                open(pe_2, "w") as pe_o2,\
                open(se_1, "w") as se_o1,\
                open(se_2, "w") as se_o2:
            krill(q_min, q_percent, p64, f1, f2, pe_o1, pe_o2, se_o1, se_o2)


def krill(q_min, q_percent, p64, f1, f2, pe_o1, pe_o2, se_o1, se_o2):
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    if p64:
        val = dict(zip(scores[31:95], range(0, 43)))
    else:
        val = dict(zip(scores[:43], range(0, 43)))
    y, entry1, entry2 = 0, "", ""
    for line1, line2 in zip(f1, f2):
        y += 1
        line1 = line1.rstrip()
        line2 = line2.rstrip()
        entry1 = entry1 + line1 + "\n"
        entry2 = entry2 + line2 + "\n"
        if y == 4:
            krill1 = krill_test(line1, q_min, q_percent, val)
            krill2 = krill_test(line2, q_min, q_percent, val)
            if krill1 is False and krill2 is False:
                pe_o1.write(entry1)
                pe_o2.write(entry2)
            elif krill1 is False and krill2 is True:
                se_o1.write(entry1)
            elif krill1 is True and krill2 is False:
                se_o2.write(entry2)
            elif krill1 is True and krill2 is True:
                pass
            y, entry1, entry2 = 0, "", ""


def krill_single_open(q_min, q_percent, p64, in1, se_1):
    '''
    quality filter test for single-end reads
    '''
    compressed = gzip_test(in1)
    if compressed:
        se_1, _ = os.path.splitext(se_1)
        with gzip.open(in1, 'rt') as f1, open(se_1, "w") as se_o1:
            krill_single(q_min, q_percent, p64, f1, se_o1)
    else:
        with open(in1) as f1, open(se_1, "w") as se_o1:
            krill_single(q_min, q_percent, p64, f1, se_o1)



def krill_single(q_min, q_percent, p64, f1, se_o1):
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    if p64:
        val = dict(zip(scores[31:95], range(0, 43)))
    else:
        val = dict(zip(scores[:43], range(0, 43)))
    y, entry1 = 0, ""
    for line1 in f1:
        y += 1
        line1 = line1.rstrip()
        entry1 = entry1 + line1 + "\n"
        if y == 4:
            krill1 = krill_test(line1, q_min, q_percent, val)
            if krill1 is False:
                se_o1.write(entry1)
            elif krill1 is True:
                pass
            y, entry1 = 0, ""


def krill_test(line, q_min, q_percent, val):
    x, fail_count, krill = 0, 0, False
    fail_base = ((100-q_percent)/100)*(len(line))
    for x in range((len(line) - 1), -1, -1):
        if int(val[line[x]]) < q_min:
            fail_count += 1
        if fail_count > fail_base:
            krill = True
            break
    return krill


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='filter fastq reads for quality')
    parser.add_argument('-r1', type=str, required=True,
            help='the full or relative path to R1 fastq file')
    parser.add_argument('-r2', type=str,
            help='the full or relative path to R2 fastq file (optional)')
    parser.add_argument('-q', type=int, required=True,
            help='the minimum qscore required at a position (integer)')
    parser.add_argument('-p', type=int, required=True,
            help='the percent frequency minimum qscore must occur per read (integer)')
    parser.add_argument('-o', type=str,
            help='the full path to output directory (optional)')
    parser.add_argument('-s', type=str,
            help='type True for phred 64 samples (optional, default phred 33)')
    args = parser.parse_args()
    krill_main()
