import sys
import os


def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    in1 = sys.argv[1]
    front_trim = int(sys.argv[2])
    back_trim = int(sys.argv[3])
    out1 = sys.argv[4]
    proj_dir = os.path.abspath(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


def scallop_comp(front_trim, back_trim, proj_dir, in1):
    '''
    composer entry point to scallop
    '''
    out1 = proj_dir + '/' + os.path.basename(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


def scallop(in1, front_trim, back_trim, proj_dir, out1):
    '''
    trim defined base numbers from the front or from the end of reads
    '''
    if back_trim == 0:
        back_trim = None
    with open(in1) as f, open(out1, 'w') as o:
        i = 0
        for line in f:
            i += 1
            if i % 2 == 0:
                o.write(line.rstrip()[front_trim:back_trim] + "\n")
            else:
                o.write(line)


def scallop_end(proj_dir, q_min, in1):
    '''
    composer entry point for trimming based on qc stats
    '''
    in1_path, in1_file = os.path.split(in1)
    _, in1_pe_se = os.path.split(in1_path)
    proj_dir += '/' + in1_pe_se
    in1_scores = in1_path + '/qc/qscores.' + in1_file
    in1_mx = []
    with open(in1_scores) as f:
        for i in f:
            tmp = []
            i = i.rstrip()
            for j in i.split(','):
                tmp.append(int(j))
            in1_mx.append(tmp)
    max_len = len(in1_mx)
    len_ls = [sum(i) for i in in1_mx]
    try:
        if uniform:
            max_len = scallop_uniform(len_ls, max_len)
    except NameError:
        pass
    for index, i in enumerate(len_ls):
        len_ls[index] = (i + 1)/2
        if len_ls[index] % 1 == 0:
            len_ls[index] = (len_ls[index]/2)
        else:
            len_ls[index] = (len_ls[index] - 0.5)/2
    back_trim = scallop_auto(len_ls, in1_mx, max_len, q_min)
    scallop_comp(0, back_trim, proj_dir, in1)


def scallop_auto(len_ls, in1_mx, max_len, q_min):
    '''
    trim defined slice of reads, producing constant read length
    '''
    back_trim = None
    for j in range((max_len - 1), -1, -1):
        index, x, x_adjust = 0, 0, 0
        for i in in1_mx[j]:
            index += i
            if x_adjust == 0 and index == len_ls[j] - 0.5:
                x_adjust = x
            if index >= len_ls[j]:
                if x_adjust == 0:
                    target = x
                    break
                else:
                    target = (x + x_adjust)/2
                    break
            x += 1
        if target >= q_min:
            back_trim = j + 1
            break
    return back_trim


def scallop_uniform(len_ls, max_len):
    uniform = max_len
    for index, i in enumerate(len_ls):
        if i != max(len_ls):
            uniform = index
            break
    return uniform


if __name__ == '__main__':
    scallop_main()
