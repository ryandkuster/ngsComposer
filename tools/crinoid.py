import argparse
import gzip
import os
import subprocess
import sys


def crinoid_main():
    '''
    standalone, command line entry point to crinoid using stdin
    '''
    if args.a:
        visualizer(args.a)
        sys.exit()
    in1 = args.r1
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    out1 = os.path.join(proj, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(proj, 'qscores.' + os.path.basename(in1))
    subprocess.check_call(['Rscript',
            os.path.abspath(sys.argv[0])[:-17] +
            '/tests/test_packages.R'] , shell=False)
    crinoid_open(in1, out1, out2)
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/helpers/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_comp(curr, in1):
    '''
    composer entry point to crinoid
    '''
    out1 = os.path.join(curr, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(curr, 'qscores.' + os.path.basename(in1))
    crinoid_open(in1, out1, out2)
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/tools/helpers/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_open(in1, out1, out2):
    '''
    produce raw counts of nucleotide and qscore occurrences
    '''
    try:
        with gzip.open(in1, 'rt') as f:
            crinoid(f, out1, out2)
    except OSError:
        with open(in1) as f:
            crinoid(f, out1, out2)


def crinoid(f, out1, out2):
    #TODO avoid using '500' as default max_len
    #TODO add option to evaluate raw file alone (should R fail)

    max_len = 500
    score_ref_dt, score_dt, base_ref_dt, base_dt = dict_maker(max_len)

    y = 0
    for line in f:
        y += 1
        if y == 2:
            for i, base in enumerate(line.rstrip()):
                base_dt[i][base] += 1
        if y == 4:
            for i, base in enumerate(line.rstrip()):
                score_dt[i][base] += 1
            y = 0

    with open(out1, "w") as o1, open(out2, "w") as o2:
        base_mx = [[0] * len(base_dt[0]) for i in range(max_len)]
        score_mx = [[0] * len(score_dt[0]) for i in range(max_len)]
        base_mx = output_prep(base_dt, base_mx, base_ref_dt, 0)
        score_mx = output_prep(score_dt, score_mx, score_ref_dt, 0)
        base_mx = matrix_succinct(base_mx)
        score_mx = matrix_succinct(score_mx)
        matrix_print(base_mx, o1)
        matrix_print(score_mx, o2)


def dict_maker(max_len):
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    score_ref_dt = dict(zip(scores[:41], scores[-int(len(scores)/2):
                        -int(len(scores)/2)+41]))
    score_dt = {i: {j: 0 for j in score_ref_dt.keys()}
                for i in range(max_len)}
    base_ref_dt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    base_dt = {i: {j: 0 for j in base_ref_dt.keys()}
               for i in range(max_len)}
    return score_ref_dt, score_dt, base_ref_dt, base_dt


def output_prep(dt1, mx, dt2, col):
    '''
    use recursion to extract counts from nested dictionaries
    '''
    for k1, v1 in dt1.items():
        if isinstance(v1, dict):
            output_prep(v1, mx, dt2, k1)
        else:
            for k2, v2 in dt2.items():
                if k1 == k2:
                    mx[int(col)][int(v2)] = v1
    return mx


def matrix_succinct(mx):
    '''
    identify and rows in matrix containing zeroes
    '''
    for i, row in enumerate(mx):
        if sum(row) == 0:
            del mx[i:]
            break
    return mx


def matrix_print(mx, outfile):
    '''
    write matrix of raw counts to csv file
    '''
    for row in mx:
        outfile.write(",".join(str(x) for x in row))
        outfile.write("\n")


def visualizer(a):
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/helpers/qc_plots_score_only.R'] + [a, a], shell=False)


def combine_matrix(a1, a2, l):
    with open(a1) as f1, open(a2) as f2:
        a = [line.rstrip().split(',') for line in f1]
        b = [line.rstrip().split(',') for line in f2]
        if b:
            pass
        else:
            b = [[0 for j in range(len(a[0]))] for i in range(l)]
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='qc summary statistics')
    parser.add_argument('-r1', type=str,
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-o', type=str,
            help='the full path to output directory (optional)')
    parser.add_argument('-a', type=str,
            help='qscores file for visualization (optional)')
    args = parser.parse_args()
    crinoid_main()
