import sys
import os


def rotifer_main():
    '''
    standalone, command line entry point to rotifer using stdin
    '''
    input1 = sys.argv[1] # name of input file (these MUST be identical in length)
    project_dir = os.path.dirname(os.path.abspath(input1))
    overhang_list = ['TCC','TCT']
    rotifer(project_dir, overhang_list, input1)


def rotifer_pipeline():
    '''
    composer entry point to rotifer
    '''
    pass


def rotifer(project_dir, overhang_list, input1):
    with open(input1) as s:
        line1 = s.readline()
        line2 = s.readline()
    seq_len = len(line2) - 1
    input1_base = os.path.basename(input1)
    with open(input1) as f, open(project_dir + '/hit_' + str(seq_len) + '_' + input1_base, 'w') as o1, \
    open(project_dir + '/miss_' + str(seq_len) + '_' + input1_base, 'w') as o2:
        y, entry = 0, ""
        for line in f:
            y += 1
            if y == 2:
                test = 0
                for seq in overhang_list:
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