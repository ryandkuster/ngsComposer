import argparse
import gzip
import math
import multiprocessing
import os
import subprocess
import sys
from functools import partial
from itertools import islice, zip_longest

from helpers.gzip_handling import gzip_test


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
    procs = args.t if args.t else 1
    p64 = False if args.s is None else True
    out1 = os.path.join(proj, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(proj, 'qscores.' + os.path.basename(in1))
    if args.a:
        visualizer(out1, out2)
        sys.exit()
    subprocess.check_call(['Rscript',
            os.path.abspath(os.path.join(os.path.dirname(__file__),
            '..', 'tests', 'test_packages.R'))], shell=False)
    crinoid_open(in1, out1, out2, procs, p64)
    subprocess.check_call(['Rscript',
            os.path.abspath(os.path.join(os.path.dirname(__file__),
            'helpers', 'qc_plots.R'))] + [out1, out2], shell=False)


def crinoid_comp(curr, all_qc, procs, p64, in1):
    '''
    composer entry point to crinoid
    '''
    out1 = os.path.join(curr, 'nucleotides.' + os.path.basename(in1))
    out2 = os.path.join(curr, 'qscores.' + os.path.basename(in1))
    crinoid_open(in1, out1, out2, procs, p64)
    if all_qc == 'full':
        subprocess.check_call(['Rscript',
            os.path.dirname(os.path.abspath(sys.argv[0])) +
            '/tools/helpers/qc_plots.R'] + [out1, out2], shell=False)


def crinoid_open(in1, out1, out2, procs, p64):
    '''
    open as gzipped file object if gzipped
    '''
    compressed = gzip_test(in1)
    if compressed:
        with gzip.open(in1, 'rt') as f:
            k_score = kmer_test(f)
        with gzip.open(in1, 'rt') as f:
            crinoid(f, out1, out2, procs, p64, k_score)
    else:
        with open(in1) as f:
            k_score = kmer_test(f)
        with open(in1) as f:
            crinoid(f, out1, out2, procs, p64, k_score)


def kmer_test(f):
    '''
    read first __ lines and evaluate score composition
    '''
    char_count = []
    for i, line in enumerate(f):
        if (i + 1) % 4 == 0:
            for j in line.rstrip():
                if j not in char_count:
                    char_count.append(j)
        if i == 4000:
            break
    k_score = 3 if len(char_count) > 4 else 9
    return k_score


def crinoid(f, out1, out2, procs, p64, k_score):
    '''
    produce raw counts of nucleotide and qscore occurrences
    '''
    k_seq = 6

    if procs > 1:
        results = parallel_crinoid(f, k_seq, k_score, procs)
    else:
        results = motif_counter(k_seq, k_score, f)

    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    if p64:
        score_dt = dict(zip(scores[31:95], range(0, 43)))
    else:
        score_dt = dict(zip(scores[:43], range(0, 43)))
    base_dt = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    score_mx = [[0 for j in range(43)]]
    base_mx = [[0 for j in range(5)]]
    base_mx = dicto_iter(results[0], k_seq, base_mx, base_dt)
    score_mx = dicto_iter(results[1], k_score, score_mx, score_dt)
    base_mx = matrix_succinct(base_mx)
    score_mx = matrix_succinct(score_mx)

    with open(out1, "w") as o1, open(out2, "w") as o2:
        matrix_print(base_mx, o1)
        matrix_print(score_mx, o2)


class Test:
    pass


def parallel_crinoid(f, k_seq, k_score, procs):
    kmer_seq_dt = {}
    kmer_score_dt = {}
    total = []
    subset = round(1000000/procs) * 4
    for i in range(procs):
        text_slice(f, subset)
        total.append(Test.next_lines)
    pool = multiprocessing.Pool(procs)
    party = partial(motif_counter, k_seq, k_score)
    sub_vals = pool.map(party, total)
    pool.close()

    for i in sub_vals:
        sub_seqs = i[0]
        sub_scores = i[1]
        kmer_seq_dt = motif_total(sub_seqs, kmer_seq_dt)
        kmer_score_dt = motif_total(sub_scores, kmer_score_dt)

    while Test.next_lines != []:
        total = []
        for i in range(procs):
            text_slice(f, subset)
            total.append(Test.next_lines)
        pool = multiprocessing.Pool(procs)
        party = partial(motif_counter, k_seq, k_score)
        sub_vals = pool.map(party, total)
        pool.close()
        for i in sub_vals:
            sub_seqs = i[0]
            sub_scores = i[1]
            kmer_seq_dt = motif_total(sub_seqs, kmer_seq_dt)
            kmer_score_dt = motif_total(sub_scores, kmer_score_dt)
    return [kmer_seq_dt, kmer_score_dt]


def text_slice(file_opened, n):
    '''
    grab next n lines from input fastq
    '''
    Test.next_lines = [x for x in islice(file_opened, n)]


def motif_counter(k_seq, k_score, total):
    '''
    subprocess of crinoid to count seq, score motifs
    '''
    kmer_seq_dt, kmer_score_dt = {}, {}
    y = 0
    for line in total:
        y += 1
        if y == 2:
            okay_mer(line.rstrip(), k_seq, kmer_seq_dt)
        if y == 4:
            okay_mer(line.rstrip(), k_score, kmer_score_dt)
            y = 0
    return [kmer_seq_dt, kmer_score_dt]


def motif_total(i, dt):
    '''
    add newly-discovered motif counts to grand total dictionaries
    '''
    for k, v in i.items():
        if k in dt:
            dt[k] = [a + b for a, b in zip_longest(v, dt[k], fillvalue=0)]
        else:
            dt[k] = v
    return dt


def okay_mer(line, k, motif_dt):
    '''
    subset line by k and count motifs in dictionaries
    '''
    segment = math.ceil(len(line.rstrip())/k)
    for i in range(segment):
        s = i * k
        f = s + k
        motif = line[s:f]
        if motif in motif_dt:
            try:
                motif_dt[motif][i] += 1
            except IndexError:
                diff = i - (len(motif_dt[motif])) + 1
                for l in range(diff):
                    motif_dt[motif].append(0)
                motif_dt[motif][i] += 1
        else:
            motif_dt[motif] = [0 for seg in range(segment)]
            motif_dt[motif][i] += 1


def dicto_iter(motif_dt, k, mx, ref_dt):
    '''
    transform the dictionaries into usable matrices for R plots
    '''
    for motif, count in motif_dt.items():
        for i, window in enumerate(count):
            for j, char in enumerate(motif):
                pos = (i * k) + j

                try:
                    mx[pos][ref_dt[char]] += window
                except IndexError:
                    mx = bespoke_matrix(mx, pos)
                    mx[pos][ref_dt[char]] += window
    return mx


def bespoke_matrix(mx, pos):
    '''
    add row to matrix to accomodate variance in read lengths
    '''
    while len(mx) < pos + 1:
        mx.append([0 for j in range(len(mx[0]))])
    return mx


def matrix_print(mx, outfile):
    '''
    write matrix of raw counts to csv file
    '''
    for row in mx:
        outfile.write(",".join(str(x) for x in row))
        outfile.write("\n")


def matrix_succinct(mx):
    '''
    identify and remove rows containing zeroes
    '''
    for i, row in enumerate(mx):
        if sum(row) == 0:
            del mx[i:]
            break
    return mx


def visualizer(out1, out2):
    '''
    produce images using existing, raw nucleotide and qscore files
    '''
    subprocess.check_call(['Rscript',
           os.path.dirname(os.path.abspath(__file__)) +
           '/helpers/qc_plots.R'] + [out1, out2], shell=False)


def combine_matrix(in_ls, out):
    '''
    combine values from previously generated qc matrices
    '''
    out_q = os.path.join(os.path.dirname(in_ls[0]), 'qc', 'qscores.' + out)
    out_n = os.path.join(os.path.dirname(in_ls[0]), 'qc', 'nucleotides.' + out)
    q2 = [[0 for j in range(45)] for i in range(1000)]
    n2 = [[0 for j in range(5)] for i in range(1000)]

    for i in in_ls:
        q1 = os.path.join(os.path.dirname(i), 'qc', 'qscores.' +
                          os.path.basename(i))
        q2 = matrix_mash(q1, q2)
        n1 = os.path.join(os.path.dirname(i), 'qc', 'nucleotides.' +
                          os.path.basename(i))
        n2 = matrix_mash(n1, n2)

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


def matrix_mash(a1, a2):
    '''
    modify existing qscore files into matrices, combine with a2
    '''
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
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    parser.add_argument('-a', type=str, metavar='',
            help='create visualizations on existing qscore and nucleotide \
                raw data (optional, use "True", requires base name for -r)')
    parser.add_argument('-s', type=str, metavar='',
            help='type True for phred 64 samples (optional, default phred 33)')
    parser.add_argument('-t', type=int, metavar='',
            help='number of subprocesses')
    args = parser.parse_args()
    crinoid_main()
