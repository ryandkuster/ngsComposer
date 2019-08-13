import argparse
import gzip
import os
import subprocess
import sys


def crinoid_main():
    '''
    standalone, command line entry point to crinoid using stdin
    '''
    in1 = args.r1
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    p64 = False if args.s is None else True
    out1 = os.path.join(proj, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(proj, 'qscores.' + os.path.basename(in1))
    if args.a:
        visualizer(out1, out2)
        sys.exit()
    subprocess.check_call(['Rscript',
            os.path.abspath(sys.argv[0])[:-17] +
            '/tests/test_packages.R'], shell=False)
    crinoid_open(in1, out1, out2, p64)
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/helpers/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_comp(curr, all_qc, p64, in1):
    '''
    composer entry point to crinoid
    '''
    out1 = os.path.join(curr, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(curr, 'qscores.' + os.path.basename(in1))
    crinoid_open(in1, out1, out2, p64)
    if all_qc == 'full':
        subprocess.check_call(['Rscript',
            os.path.dirname(os.path.abspath(sys.argv[0])) +
            '/tools/helpers/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_open(in1, out1, out2, p64):
    '''
    produce raw counts of nucleotide and qscore occurrences
    '''
    try:
        with gzip.open(in1, 'rt') as f:
            crinoid(f, out1, out2, p64)
    except OSError:
        with open(in1) as f:
            crinoid(f, out1, out2, p64)


def crinoid(f, out1, out2, p64):
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    if p64:
        score_dt = dict(zip(scores[31:95], range(0, 43)))
    else:
        score_dt = dict(zip(scores[:43], range(0, 43)))
    base_dt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    score_mx = [[0 for j in range(43)]]
    base_mx = [[0 for j in range(5)]]
    y = 0
    for line in f:
        y += 1
        if y == 2:
            for i, base in enumerate(line.rstrip()):
                try:
                    base_mx[i][base_dt[base]] += 1
                except IndexError:
                    base_mx = bespoke_matrix(base_mx, 5)
                    base_mx[i][base_dt[base]] += 1
        if y == 4:
            for i, base in enumerate(line.rstrip()):
                try:
                    score_mx[i][score_dt[base]] += 1
                except IndexError:
                    score_mx = bespoke_matrix(score_mx, 43)
                    score_mx[i][score_dt[base]] += 1
            y = 0
    with open(out1, "w") as o1, open(out2, "w") as o2:
        matrix_print(base_mx, o1)
        matrix_print(score_mx, o2)


def bespoke_matrix(mx, max_len):
    mx.append([0 for j in range(max_len)])
    return mx


def matrix_print(mx, outfile):
    '''
    write matrix of raw counts to csv file
    '''
    for row in mx:
        outfile.write(",".join(str(x) for x in row))
        outfile.write("\n")


def visualizer(out1, out2):
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(__file__)) +
           '/helpers/qc_plots.R'] + [out1, out2], shell=False)


def combine_matrix(in_ls, out):
    '''
    combine values from previously generated qc matrices
    '''
    # TODO update for alternate scoring systems
    out_q = os.path.join(os.path.dirname(in_ls[0]), 'qc', 'qscores.' + out)
    out_n = os.path.join(os.path.dirname(in_ls[0]), 'qc', 'nucleotides.' + out)
    q2 = [[0 for j in range(41)] for i in range(1000)]
    n2 = [[0 for j in range(5)] for i in range(1000)]

    for i in in_ls:
        q1 = os.path.join(os.path.dirname(i), 'qc', 'qscores.' +
                          os.path.basename(i))
        q2 = matrix_mash(in_ls, q1, q2)
        n1 = os.path.join(os.path.dirname(i), 'qc', 'nucleotides.' +
                          os.path.basename(i))
        n2 = matrix_mash(in_ls, n1, n2)

    for i in reversed(q2):
        if sum(i) == 0:
            q2.pop()

    for i in reversed(n2):
        if sum(i) == 0:
            n2.pop()

    with open(out_q, 'w') as o1, open(out_n, 'w') as o2:
        for row in q2:
            o1.write(",".join(str(x) for x in row))
            o1.write("\n")

        for row in n2:
            o2.write(",".join(str(x) for x in row))
            o2.write("\n")

    visualizer(out_n, out_q)


def matrix_mash(in_ls, a1, a2):
    with open(a1) as f1:
        a1 = [line.rstrip().split(',') for line in f1]
        for i, col in enumerate(a2):
            for j, score in enumerate(col):
                try:
                    a2[i][j] = int(a2[i][j]) + int(a1[i][j])
                except IndexError:
                    break
    return a2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='qc summary statistics')
    parser.add_argument('-r1', type=str,
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-o', type=str,
            help='the full path to output directory (optional)')
    parser.add_argument('-a', type=str,
            help='create visualizations on existing qscore and nucleotide \
                raw data (optional, use "True", requires base name for -r1)')
    parser.add_argument('-s', type=str,
            help='type True for phred 64 samples (optional, default phred 33)')

    args = parser.parse_args()
    crinoid_main()
