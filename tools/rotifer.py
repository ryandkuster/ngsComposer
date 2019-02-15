import sys
import os


def rotifer_main():
    '''
    standalone, command line entry point to rotifer using stdin
    '''
    in1 = sys.argv[1]
    proj_dir = os.path.dirname(os.path.abspath(in1))
    in1 = sys.argv[2]
    try:
        in2 = sys.argv[3]
    except IndexError:
        in2 = None
    bases_ls = ['CCC', 'CCT']
    pe_1 = 'pe.' + os.path.basename(in1)
    pe_2 = 'pe.' + os.path.basename(in2)
    se_1 = 'se.' + os.path.basename(in1)
    se_2 = 'se.' + os.path.basename(in2)
    rotifer(proj_dir, bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)


def rotifer_comp(in1_ls, in2_ls, bases_ls, non_genomic, proj_dir_current, in1):
    '''
    composer entry point to rotifer
    '''
    pe_1 = proj_dir_current + '/paired/' + os.path.basename(in1)
    se_1 = proj_dir_current + '/single/' + os.path.basename(in1)
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_2 = proj_dir_current + '/paired/' + os.path.basename(in2)
        se_2 = proj_dir_current + '/single/' + os.path.basename(in2)
        rotifer(bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)
    except ValueError:
        rotifer(bases_ls, in1, pe_1, se_1)


def rotifer(bases_ls, in1, in2, pe_1, pe_2, se_1, se_2):
    '''
    parse single and paired-end reads for recognized motifs
    '''
    with open(in1) as f1, open(in2) as f2,\
            open(pe_1, "w") as pe_o1,\
            open(pe_2, "w") as pe_o2,\
            open(se_1, "w") as se_o1,\
            open(se_2, "w") as se_o2:
        y, entry1, entry2 = 0, "", ""
        for line1, line2 in zip(f1, f2):
            y += 1
            line1 = line1.rstrip()
            line2 = line2.rstrip()
            entry1 = entry1 + line1 + "\n"
            entry2 = entry2 + line2 + "\n"
            if y == 2:
                rotifer1 = rotifer_test(line1, bases_ls)
                rotifer2 = rotifer_test(line2, bases_ls)
            if y == 4:
                if rotifer1 is False and rotifer2 is False:
                    pe_o1.write(entry1)
                    pe_o2.write(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is False and rotifer2 is True:
                    se_o1.write(entry1)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is True and rotifer2 is False:
                    se_o2.write(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is True and rotifer2 is True:
                    y, entry1, entry2 = 0, "", ""


def rotifer_single(bases_ls, in1, pe_1, se_1):
    '''
    parse single-end reads for recognized motifs
    '''
    with open(in1) as f1, open(pe_1, "w") as pe_o1, open(se_1, "w") as se_o1:
        y, entry1 = 0, ""
        for line1 in f1:
            y += 1
            line1 = line1.rstrip()
            entry1 = entry1 + line1 + "\n"
            if y == 2:
                rotifer1 = rotifer_test(line1, bases_ls)
            if y == 4:
                if rotifer1 is False:
                    pe_o1.write(entry1)
                    y, entry1 = 0, ""
                elif rotifer1 is True:
                    se_o1.write(entry1)
                    y, entry1 = 0, ""


def rotifer_test(line, bases_ls):
    rotifer = True
    for seq in bases_ls:
        if line.startswith(seq):
            rotifer = False
    return rotifer


if __name__ == "__main__":
    rotifer_main()
