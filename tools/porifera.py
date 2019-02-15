import sys
import os


def main():
    proj_dir = sys.argv[1]
    in_list = []
    print(proj_dir)
    for filename in os.listdir(proj_dir):
        if filename.startswith('qscores_'):
            in_list.append(proj_dir + filename)
    
    for i in in_list:
        with open(i) as f:
            for j, line in enumerate(f):
                print(line)
        print(j + 1)
    max_len = 5
    header = ['pos', 'mean', 'lower', 'q1', 'med', 'q3', 'upper', 'A', 'C', 'T', 'G', 'N']
    output_matrix = [[None] * len(header) for i in range(max_len)]
    output_matrix[0] = header
#    porifera(score_mx)


def porifera(score_mx):
    median_index = (sum(score_mx[0])+1)/2
    print(median_index)
    if median_index % 1 == 0:
        print('^^^ this is an integer')
        q1_index = median_index - (median_index/2)
        q3_index = median_index + (median_index/2)
    else:
        q1_index = median_index - ((median_index - 0.5)/2)
        q3_index = median_index + ((median_index - 0.5)/2)
    print(q1_index)
    print(q3_index)


def stats(target_index, col):
    for j in range(seq_len):
        index = 0
        x = 0
        x_adjust = 0
        for i in value_matrix[j]:
            index = index + i
            if (x_adjust == 0 and index == target_index - 0.5):
                x_adjust = x
            if index >= target_index:
                if x_adjust == 0:
                    target = x
                    break
                else:
                    target = (x + x_adjust)/2
                    break
            x += 1
        output_matrix[j+1][col] = target


if __name__ == "__main__":
    main()
