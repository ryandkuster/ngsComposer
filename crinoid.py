import sys
import os


def crinoid_main():
    in1 = sys.argv[1] # input sequences
    proj_dir = os.path.dirname(os.path.abspath(in1))
    out1 = proj_dir + '/nucleotides_' + os.path.basename(in1) + '.txt'
    out2 = proj_dir + '/qscores_' + os.path.basename(in1) + '.txt'
    crinoid(in1, out1, out2)


def crinoid_comp(proj_dir, in1):
    out1 = proj_dir + '/nucleotides_' + os.path.basename(in1) + '.txt'
    out2 = proj_dir + '/qscores_' + os.path.basename(in1) + '.txt'
    crinoid(in1, out1, out2)


def crinoid(in1, out1, out2):
    scores = open(os.path.dirname(os.path.abspath(__file__)) + '/scores.txt').read().split()
    score_ref_dt = dict(zip(scores[:41],scores[-int(len(scores)/2):-int(len(scores)/2)+41]))
    max_len = max_len_count(in1)
    with open(in1) as f, open(out1, "w") as o1, open(out2, "w") as o2:
        y = 0
        base_ref_dt = {'A': 0, 'C': 1, 'G': 2, 'T' : 3, 'N': 4}
        base_dt = {i: {j: 0 for j in base_ref_dt.keys()} for i in range(max_len)}
        score_dt = {i: {j: 0 for j in score_ref_dt.keys()} for i in range(max_len)}
        for line in f:
            y += 1
            if y == 2:
                for i, base in enumerate(line.rstrip()):
                    base_dt[i][base] +=1
            if y == 4:
                for i, base in enumerate(line.rstrip()):
                    score_dt[i][base] +=1
                y = 0
        base_mx = [[0] * len(base_dt[0]) for i in range(max_len)]
        score_mx = [[0] * len(score_dt[0]) for i in range(max_len)]
        base_mx = output_prep(base_dt, base_mx, base_ref_dt, 0)
        score_mx = output_prep(score_dt, score_mx, score_ref_dt, 0)
        matrix_print(base_mx, o1)
        matrix_print(score_mx, o2)


def max_len_count(in1):
    with open(in1) as f:
        max_len, y = 0, 0
        for line in f:
            y += 1
            if y == 2:
                len_test = len(line.rstrip())
                if len_test > max_len:
                    max_len = len_test
            if y == 4:
                y = 0
    return max_len


def output_prep(dt1, mx, dt2, col):
    for k1, v1 in dt1.items():
        if isinstance(v1, dict):
            output_prep(v1, mx, dt2, k1)
        else:
            for k2, v2 in dt2.items():
                if k1 == k2:
                    mx[int(col)][int(v2)] = v1
    return mx


def matrix_print(mx, outfile):
    for row in mx:
        outfile.write(",".join(str(x) for x in row))
        outfile.write("\n")


if __name__ == "__main__":
    crinoid_main()
