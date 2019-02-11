import linecache, sys, os


def crinoid_main():
    fastq = sys.argv[1] # input sequences
    chunk = int(sys.argv[2]) # this determines the size of fastq scores to hold in a matrix, use 20000 for now
    output1 = sys.argv[3] # name of output stats file for quality score counts
    output2 = sys.argv[4] # name of output stats file for nucleotide counts


with open(fastq) as s:
    line1 = s.readline()
    line2 = s.readline()
seq_len = len(line2) - 1

scores = open(os.path.dirname(os.path.abspath(__file__)) + '/scores.txt').read().split()
val = dict(zip(scores[:41],scores[-int(len(scores)/2):-int(len(scores)/2)+41]))
val_alpha = {"A": 0 , "C" : 1 , "G" : 2 , "T" : 3 , "N" : 4}
y, i = 0, 0
count_matrix = [[None] * chunk for i in range(seq_len)]
value_matrix = [[0] * len(val) for i in range(seq_len)]
count_matrix_alpha = [[None] * chunk for i in range(seq_len)]
value_matrix_alpha = [[0] * len(val_alpha) for i in range(seq_len)]

def vertical_mine():
    for row in range(seq_len):
        count_dict = counter(count_matrix[row])
        for k1, v1 in count_dict.items():
            for k2, v2 in val.items():
                if k1 == k2:
                    value_matrix[row][int(v2)] += v1
        count_dict_alpha = counter(count_matrix_alpha[row])
        for k1, v1 in count_dict_alpha.items():
            for k2, v2 in val_alpha.items():
                if k1 == k2:
                    value_matrix_alpha[row][int(v2)] += v1

def counter(row):
    count_dict = {}
    for character in row:
        if character in count_dict:
            count_dict[character] += 1
        else:
            count_dict[character] = 1
    return count_dict

with open(fastq) as f, open(output1, "w") as o1, open(output1, "w") as o2:
    for line in f:
        y += 1
        if y == 2:
            for j in range(seq_len):
                count_matrix_alpha[j][i] = line[j]
        if y == 4:
            for j in range(seq_len):
                count_matrix[j][i] = line[j]
            y = 0
            i += 1
        if i == chunk:
            i = 0
            vertical_mine()
            count_matrix = [[None] * chunk for i in range(seq_len)]
            count_matrix_alpha = [[None] * chunk for i in range(seq_len)]
    vertical_mine()
    for row in value_matrix:
        o1.write(",".join(str(x) for x in row))
        o1.write("\n")
    for row in value_matrix_alpha:
        o2.write(",".join(str(x) for x in row))
        o2.write("\n")

if __name__ == "__main__":
    crinoid_main()



