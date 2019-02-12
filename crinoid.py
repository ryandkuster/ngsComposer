import sys
import os


def crinoid_main():
    in1 = sys.argv[1] # input sequences
    proj_dir = os.path.dirname(os.path.abspath(in1))
    crinoid_init(in1, proj_dir)


def crinoid_comp():
    pass
    crinoid_init(in1, proj_dir)


def crinoid_init(in1, proj_dir):
    out1 = proj_dir + '/nucleotides_' + os.path.basename(in1) + '.txt'
    out2 = proj_dir + '/qscores_' + os.path.basename(in1) + '.txt'
    scores = open(os.path.dirname(os.path.abspath(__file__)) + '/scores.txt').read().split()
    score_dict = dict(zip(scores[:41],scores[-int(len(scores)/2):-int(len(scores)/2)+41]))
    max_len = max_len_count(in1)
    crinoid(in1, out1, out2, max_len, score_dict)


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


def crinoid(in1, out1, out2, max_len, score_dict):
    with open(in1) as f, open(out1, "w") as o1, open(out2, "w") as o2:
        print('processing reads')
        y = 0
        base_dict = {i: {'A': 0, 'C': 0, 'G': 0, 'T' : 0, 'N': 0} for i in range(max_len)}
        qual_dict = {i: {j: 0 for j in score_dict.keys()} for i in range(max_len)}
        for line in f:
            y += 1
            if y == 2:
                for i, base in enumerate(line.rstrip()):
                    base_dict[i][base] +=1
            if y == 4:
                for i, base in enumerate(line.rstrip()):
                    qual_dict[i][base] +=1
                y = 0
        base_mx = [[0] * len(base_dict[0]) for i in range(max_len)]
        score_mx = [[] for i in range(max_len)]
        score_mx = output_prep(qual_dict, score_mx, score_dict)
        print(score_mx)


def output_prep(dict1, mx, dict2):
    for k1, v1 in dict1.items():
        if isinstance(v1, dict):
            output_prep(v1, mx, dict2)
        else:
            for k2, v2 in dict2.items():
                if k1 == k2:
                    print(k1)
                    print(k2)
                    mx[int(v2)].append(v1)
    return mx


if __name__ == "__main__":
    crinoid_main()
