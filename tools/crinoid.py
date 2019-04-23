import sys
import os
import subprocess
import argparse


def crinoid_main():
    in1 = args.r1
    proj_dir = os.path.dirname(os.path.abspath(in1))
    out1 = proj_dir + '/nucleotides.' + os.path.basename(in1)
    out2 = proj_dir + '/qscores.' + os.path.basename(in1)
    crinoid(in1, out1, out2)
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_comp(proj_dir, in1):
    '''
    composer entry point to crinoid
    '''
    out1 = proj_dir + '/nucleotides.' + os.path.basename(in1)
    out2 = proj_dir + '/qscores.' + os.path.basename(in1)
    crinoid(in1, out1, out2)
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(sys.argv[0])) +
           '/tools/qc_plots.R'] + [out1, out2], shell=False)


def crinoid(in1, out1, out2):
    '''
    produce raw counts of nucleotide and qscore occurrences
    '''
    #TODO avoid using '500' as default max_len
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/scores.txt').read().split()
    score_ref_dt = dict(zip(scores[:41], scores[-int(len(scores)/2):
                        -int(len(scores)/2)+41]))
    max_len = 500
    with open(in1) as f, open(out1, "w") as o1, open(out2, "w") as o2:
        y = 0
        base_ref_dt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
        base_dt = {i: {j: 0 for j in base_ref_dt.keys()}
                   for i in range(max_len)}
        score_dt = {i: {j: 0 for j in score_ref_dt.keys()}
                    for i in range(max_len)}
        for line in f:
            y += 1
            if y == 2:
                for i, base in enumerate(line.rstrip()):
                    base_dt[i][base] += 1
            if y == 4:
                for i, base in enumerate(line.rstrip()):
                    score_dt[i][base] += 1
                y = 0
        base_mx = [[0] * len(base_dt[0]) for i in range(max_len)]
        score_mx = [[0] * len(score_dt[0]) for i in range(max_len)]
        base_mx = output_prep(base_dt, base_mx, base_ref_dt, 0)
        score_mx = output_prep(score_dt, score_mx, score_ref_dt, 0)
        base_mx = matrix_succinct(base_mx)
        score_mx = matrix_succinct(score_mx)
        matrix_print(base_mx, o1)
        matrix_print(score_mx, o2)


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='qc summary statistics')
    parser.add_argument('-r1', type=str,
            help='the full or relative path to R1 or R2 fastq file')
    args = parser.parse_args()
    crinoid_main()
