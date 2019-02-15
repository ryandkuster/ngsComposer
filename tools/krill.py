import sys
import os


def krill_main():
    '''
    standalone, command line entry point to krill using stdin
    '''
    q_min = int(sys.argv[1]) # the minimum phred33 value to accept
    q_percent = int(sys.argv[2]) # the percentage of the read that must have the above value (allows variable length reads)
    in1 = sys.argv[3]
    try:
        in2 = sys.argv[4]
    except IndexError:
        in2 = None
    pe_1 = 'pe.' + os.path.basename(in1)
    pe_2 = 'pe.' + os.path.basename(in2)
    se_1 = 'se.' + os.path.basename(in1)
    se_2 = 'se.' + os.path.basename(in2)
    krill(q_min, q_percent, in1, in2, pe_1, pe_2, se_1, se_2)


def krill_comp(in1_ls, in2_ls, q_min, q_percent, proj_dir_current, in1):
    '''
    composer entry point to krill
    '''
    pe_1 = proj_dir_current + '/paired/' + os.path.basename(in1)
    se_1 = proj_dir_current + '/single/pe_lib/' + os.path.basename(in1)
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_2 = proj_dir_current + '/paired/' + os.path.basename(in2)
        se_2 = proj_dir_current + '/single/pe_lib/' + os.path.basename(in2)
        krill(q_min, q_percent, in1, in2, pe_1, pe_2, se_1, se_2)
    except ValueError:
        se_1 = proj_dir_current + '/single/se_lib/' + os.path.basename(in1)
        krill_single(q_min, q_percent, in1, se_1)


def krill(q_min, q_percent, in1, in2, pe_1, pe_2, se_1, se_2):
    '''
    quality filter test for single and paired-end reads
    '''
    with open(in1) as f1, open(in2) as f2,\
            open(pe_1, "w") as pe_o1,\
            open(pe_2, "w") as pe_o2,\
            open(se_1, "w") as se_o1,\
            open(se_2, "w") as se_o2:
        scores = open(os.path.dirname(os.path.abspath(__file__)) +
                      '/scores.txt').read().split()
        val = dict(zip(scores[:int(len(scores)/2)],scores[-int(len(scores)/2):]))
        y, entry1, entry2 = 0, "", ""
        for line1, line2 in zip(f1, f2):
            y += 1
            line1 = line1.rstrip()
            line2 = line2.rstrip()
            entry1 = entry1 + line1 + "\n"
            entry2 = entry2 + line2 + "\n"
            if y == 4:
                krill1 = krill_test(line1, q_min, q_percent, val)
                krill2 = krill_test(line2, q_min, q_percent, val)
                if krill1 is False and krill2 is False:
                    pe_o1.write(entry1)
                    pe_o2.write(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif krill1 is False and krill2 is True:
                    se_o1.write(entry1)
                    y, entry1, entry2 = 0, "", ""
                elif krill1 is True and krill2 is False:
                    se_o2.write(entry2)
                    y, entry1, entry2 = 0, "", ""
                elif krill1 is True and krill2 is True:
                    y, entry1, entry2 = 0, "", ""


def krill_single(q_min, q_percent, in1, se_1):
    '''
    quality filter test for single-end reads
    '''
    with open(in1) as f1, open(se_1, "w") as se_o1:
        scores = open(os.path.dirname(os.path.abspath(__file__)) +
                      '/scores.txt').read().split()
        val = dict(zip(scores[:int(len(scores)/2)],scores[-int(len(scores)/2):]))
        y, entry1 = 0, ""
        for line1 in f1:
            y += 1
            line1 = line1.rstrip()
            entry1 = entry1 + line1 + "\n"
            if y == 4:
                krill1 = krill_test(line1, q_min, q_percent, val)
                if krill1 is False:
                    se_o1.write(entry1)
                    y, entry1 = 0, ""
                elif krill1 is True:
                    y, entry1 = 0, ""


def krill_test(line, q_min, q_percent, val):
    x, fail_count, krill = 0, 0, False
    fail_base = ((100-q_percent)/100)*(len(line))
    for x in range((len(line) - 1),-1,-1):
        if int(val[line[x]]) < q_min:
            fail_count += 1
        if fail_count > fail_base:
            krill = True
            break
    return krill


if __name__ == "__main__":
    krill_main()
