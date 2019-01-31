import sys
import os


def rotifer_main():
    '''
    standalone, command line entry point to rotifer using stdin
    '''
    in1 = sys.argv[1]
    proj_dir = os.path.dirname(os.path.abspath(in1))
    overhang_ls = ['TCC', 'TCT']
    rotifer(proj_dir, overhang_ls, in1)


def rotifer_comp():
    '''
    composer entry point to rotifer
    '''
    pass


def rotifer(proj_dir, overhang_ls, in1):
    with open(in1) as f:
        for i, line in enumerate(f):
            if i == 1:
                seq_len = len(line.rstrip())
    in1_base = os.path.basename(in1)
    with open(in1) as f, \
            open(proj_dir + '/hit_' + str(seq_len) +
                 '_' + in1_base, 'w') as o1, \
            open(proj_dir + '/miss_' + str(seq_len) +
                 '_' + in1_base, 'w') as o2:
        y, entry = 0, ""
        for line in f:
            y += 1
            if y == 2:
                test = 0
                for seq in overhang_ls:
                    if line.startswith(seq):
                        test = 1
            entry = entry + line
            if y == 4 and test == 1:
                o1.write(entry)
                y, test, entry = 0, 0, ""
            if y == 4 and test == 0:
                o2.write(entry)
                y, test, entry = 0, 0, ""


if __name__ == "__main__":
    rotifer_main()
