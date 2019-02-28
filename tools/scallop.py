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
    else:
        back_trim = - back_trim
    with open(in1) as f, open(out1, 'w') as o:
        for i, line in enumerate(f):
            i += 1
            if i % 2 == 0:
                line = line.rstrip()[front_trim:back_trim] + "\n"
                o.write(line)
            else:
                o.write(line)


def scallop_end(): #TODO modify this to be second composer entry point
    '''
    composer entry point for trimming based on qc stats
    '''
    in1 = sys.argv[1]
    qmin = int(sys.argv[2])
    in1_path, in1_file = os.path.split(in1)
    in1_scores = in1_path + '/qc/qscores.' + in1_file
    print(in1_scores)
    in1_mx = []
    with open(in1_scores) as f:
        for i in f:
            tmp = []
            i = i.rstrip()
            for j in i.split(','):
                tmp.append(int(j))
            in1_mx.append(tmp)
    max_len = len(in1_mx)
    uniform = max_len
    len_ls = [sum(i) for i in in1_mx]
    for index, i in enumerate(len_ls):
        if i != max(len_ls):
            uniform = index
            break
    print(str(uniform) + ' is the last base out of ' + str(max_len) + ' bases to keep for uniform length')
    for index, i in enumerate(len_ls):
        len_ls[index] = (i + 1)/2
        if len_ls[index] % 1 == 0:
            len_ls[index] = (len_ls[index]/2)
        else:
            len_ls[index] = (len_ls[index] - 0.5)/2
    scallop_uniform(len_ls, in1_mx, uniform, qmin)


def scallop_uniform(len_ls, in1_mx, uniform, qmin):
    '''
    trim defined slice of reads, producing constant read length
    '''
    #TODO run essentially identical function as 'scallop' without negative indexing
    for j in range((uniform - 1), -1, -1):
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
        if target >= qmin:
            print(str(j + 1) + ' is score ' + str(target) + ' and is acceptable')
        else:
            print(str(j + 1) + ' is score ' + str(target))


if __name__ == '__main__':
    scallop_main()
