import sys
import os


def rotifer_main():
    '''
    standalone, command line entry point to rotifer using stdin
    '''
    non_genomic = sys.argv[1]
    in1 = sys.argv[2]
    try:
        in2 = sys.argv[3]
    except IndexError:
        in2 = None
    proj_dir = os.path.dirname(os.path.abspath(in1))
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
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_2 = proj_dir_current + '/pe.' + os.path.basename(in2)
        se_2 = proj_dir_current + '/se.' + os.path.basename(in2)
    except IndexError:
        in2 = False
    pe_1 = proj_dir_current + '/pe.' + os.path.basename(in1)
    se_1 = proj_dir_current + '/se.' + os.path.basename(in1)
    if in2 is False:
        pass
        #TODO create single function
    else:
        rotifer(proj_dir_current, bases_ls, in1, in2, pe_1, pe_2, se_1, se_2)


def rotifer(proj_dir, bases_ls, in1, in2, pe_1, pe_2, se_1, se_2):

    #TODO implement single version; create in-line description

    with open(in1) as f1, open(in2) as f2,\
            open(pe_1, "w") as pe_o1,\
            open(pe_2, "w") as pe_o2,\
            open(se_1, "w") as se_o1,\
            open(se_2, "w") as se_o2:
        i, y, entry1, entry2 = 0, 0, "", ""
        chunk = 25000000
        pe_1_ls, pe_2_ls, se_1_ls, se_2_ls = [], [], [], []
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
                    pe_1_ls.append(entry1)
                    pe_2_ls.append(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is False and rotifer2 is True:
                    se_1_ls.append(entry1)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is True and rotifer2 is False:
                    se_2_ls.append(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif rotifer1 is True and rotifer2 is True:
                    y, entry1, entry2 = 0, "", ""
            if i == chunk:
                rotifer_unload(pe_1_ls, pe_o1)
                rotifer_unload(pe_2_ls, pe_o2)
                rotifer_unload(se_1_ls, se_o1)
                rotifer_unload(se_2_ls, se_o2)
                i = 0
        rotifer_unload(pe_1_ls, pe_o1)
        rotifer_unload(pe_2_ls, pe_o2)
        rotifer_unload(se_1_ls, se_o1)
        rotifer_unload(se_2_ls, se_o2)


def rotifer_test(line, bases_ls):
    rotifer = True
    for seq in bases_ls:
        if line.startswith(seq):
            rotifer = False
    return rotifer


def rotifer_unload(ls, f):
    f.write(str(''.join(ls)))


if __name__ == "__main__":
    rotifer_main()
