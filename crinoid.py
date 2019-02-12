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
    out1 = proj_dir + '/whatis_out1.txt'
    out2 = proj_dir + '/whatis_out2.txt'
    scores = open(os.path.dirname(os.path.abspath(__file__)) + '/scores.txt').read().split()
    val = dict(zip(scores[:41],scores[-int(len(scores)/2):-int(len(scores)/2)+41]))
    val_alpha = {"A": 0 , "C" : 1 , "G" : 2 , "T" : 3 , "N" : 4}
    chunk =  40000 # this determines the size of in1 scores to hold in a matrix, use 20000 for now
    max_len = max_len_count(in1)
    score_mx = [[None] * chunk for i in range(max_len)]
    base_mx = [[None] * chunk for i in range(max_len)]
    crinoid(in1, out1, out2, chunk, max_len, score_mx, base_mx)


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


def crinoid(in1, out1, out2, chunk, max_len, score_mx, base_mx):
    with open(in1) as f, open(out1, "w") as o1, open(out1, "w") as o2:
        y, i, base_dict, score_dict = 0, 0, {}, {}
        for line in f:
            y += 1
            if y == 2:
                for j in range(max_len):
                    base_mx[j][i] = line[j]
            if y == 4:
                for j in range(max_len):
                    score_mx[j][i] = line[j]
                y = 0
                i += 1
            if i == chunk:
                i = 0
                base_dict = vertical_mine(max_len, base_mx, base_dict)
                score_dict = vertical_mine(max_len, score_mx, score_dict)
                score_mx = [[None] * chunk for i in range(max_len)]
                base_mx = [[None] * chunk for i in range(max_len)]
        base_dict = vertical_mine(max_len, base_mx, base_dict)
        score_dict = vertical_mine(max_len, score_mx, score_dict)
#        for row in value_matrix:
#            o1.write(",".join(str(x) for x in row))
#            o1.write("\n")
#        for row in value_matrix_alpha:
#            o2.write(",".join(str(x) for x in row))
#            o2.write("\n")
        print(base_dict)
        print(score_dict)


def vertical_mine(max_len, mx, dictionary):
    for col in range(max_len):
        dictionary = counter(mx[col], dictionary)
    return dictionary
#        for k1, v1 in score_dict.items():
#            for k2, v2 in val.items():
#                if k1 == k2:
#                    value_matrix[row][int(v2)] += v1


def counter(mx_col, dictionary):
    for char in mx_col:
        if char in dictionary.keys():
            dictionary[char] += 1
        else:
            dictionary[char] = 1
    return dictionary


if __name__ == "__main__":
    crinoid_main()
