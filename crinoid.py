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
    val = dict(zip(scores[:41],scores[-int(len(scores)/2):-int(len(scores)/2)+41]))
    max_len = max_len_count(in1)
    crinoid(in1, out1, out2, max_len, val)


def max_len_count(in1):
    with open(in1) as f:
        max_len, y = 0, 0
        for line in f:
            y += 1
            if y == 2:
                len_test = len(line.rstrip())
                if len_test > max_len:
                    max_len = len_test
    return max_len


def crinoid(in1, out1, out2, max_len, val):
    with open(in1) as f, open(out1, "w") as o1, open(out2, "w") as o2:
        y = 0
        base_dict = {i: {'A': 0, 'C': 0, 'G': 0, 'T' : 0, 'N': 0} for i in range(max_len)}
        score_dict = {i: {j: 0 for j in val.keys()} for i in range(max_len)}
        for line in f:
            y += 1
            if y == 2:
                for i, base in enumerate(line.rstrip()):
                    base_dict[i][base] +=1
            if y == 4:
                for i, base in enumerate(line.rstrip()):
                    score_dict[i][base] +=1
                y = 0
        print(base_dict)
        print(score_dict)


if __name__ == "__main__":
    crinoid_main()
