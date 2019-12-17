import argparse
import gzip
import os
import sys

from helpers.gzip_handling import gzip_test


def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    in1 = args.r1
    front_trim = args.f
    end_trim = None if args.b == 0 else args.b
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    out1 = os.path.join(proj, 'trimmed.' + os.path.basename(in1))
    if args.e or args.w:
        if not args.e and args.w:
            sys.exit('-e and -w must be defined to end trim')
        end_score = args.e
        window = args.w
        scallop_end_open(in1, out1, front_trim, end_score, window) 
    else:
        scallop_open(in1, front_trim, end_trim, out1)


def scallop_comp(front_trim, end_trim, curr, in1):
    '''
    composer entry point to scallop
    '''
    out1 = os.path.join(curr, os.path.basename(in1))
    end_trim = None if end_trim == 0 else end_trim
    scallop_open(in1, front_trim, end_trim, out1)


def scallop_open(in1, front_trim, end_trim, out1):
    '''
    trim defined base numbers from the front or from the end of reads
    '''
    compressed = gzip_test(in1)
    if compressed:
        out1, _ = os.path.splitext(out1)
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            scallop(front_trim, end_trim, f, o)
    else:
        with open(in1) as f, open(out1, 'w') as o:
            scallop(front_trim, end_trim, f, o)


def scallop(front_trim, end_trim, f, o):
    i = 0
    for line in f:
        i += 1
        if i % 2 == 0:
            o.write(line.rstrip()[front_trim:end_trim] + "\n")
        else:
            o.write(line)


def scallop_end_open(in1, out1, front_trim, end_score, window):
    '''
    trim defined base numbers from the front or from the end of reads
    '''
    compressed = gzip_test(in1)
    if compressed:
        out1, _ = os.path.splitext(out1)
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            scallop_end_line(in1, front_trim, end_score, window, f, o)
    else:
        with open(in1) as f, open(out1, 'w') as o:
            scallop_end_line(in1, front_trim, end_score, window, f, o)


def scallop_end_line(in1, front_trim, end_score, window, f, o):
    scores = open(os.path.dirname(os.path.abspath(__file__)) +
                  '/helpers/scores.txt').read().split()
    val = dict(zip(scores[:43], range(0, 43)))
    good_ls = [k for k, v in val.items() if int(v) >= int(end_score)]
    i = 0
    entry = []
    for line in f:
        entry.append(line.rstrip())
        i += 1
        if i == 4:
            #TODO trim line 2 after line 4 is tested
            #implement paired ends sytem
            #add minimum length to keep
            end_trim = viewfinder(line.rstrip(), good_ls, window)
            entry[1] = entry[1][front_trim:end_trim]
            entry[3] = entry[3][front_trim:end_trim]
            for element in entry:
                o.write('%s\n' % element)
            entry = []
            i = 0


def viewfinder(line, good_ls, window):
    for i in range(len(line) - window, -1, -1):
        for j in line[i:i+window]:
            if j not in good_ls:
                break
            else:
                return i+window-1
    return 0

def scallop_end(curr, auto_trim, trim_mode, in1):
    '''
    composer entry point for trimming based on qc stats
    '''
    in1_path, in1_file = os.path.split(in1)
    _, in1_pe_se = os.path.split(in1_path)
    if in1_pe_se == 'paired' or in1_pe_se == 'single':
        curr += '/' + in1_pe_se
    in1_scores = os.path.join(in1_path, 'qc', 'qscores.' + in1_file)
    with open(in1_scores) as f:
        in1_mx = [[int(j) for j in i.rstrip().split(',')] for i in f]
    max_len = len(in1_mx)
    len_ls = [sum(i) for i in in1_mx]
    end_trim = None
    for base, i in enumerate(reversed(in1_mx)):
        total = 0
        for score, j in enumerate(i):
            total += score * j
        mean = total/sum(i)
        med = (sum(i) + 1)/2
        delta = med/2 if med % 1 == 0 else (med - 0.5)/2
        q1 = scallop_stats(med - delta, i)
        q3 = scallop_stats(med + delta, i)
        med = scallop_stats(med, i)
        lw = q1 - (1.5 * (q3 - q1))
        trim_dict = {'whisker': lw,
                     'quartile': q1,
                     'median': med,
                     'mean': mean}
        if trim_dict[trim_mode] >= auto_trim:
            end_trim = max_len - base
            break
    scallop_comp(0, end_trim, curr, in1)


def scallop_stats(target_index, pos):
    '''
    traverse score matrix per position in read for desired index
    '''
    index, x_adjust, target = 0, 0, 0
    for x, count in enumerate(pos):
        index += count
        if (x_adjust == 0 and index == target_index - 0.5):
            x_adjust = x
        if index >= target_index:
            if x_adjust == 0:
                target = x
                break
            else:
                target = (x + x_adjust)/2
                break
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='trim ends of fastq reads')
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-f', type=int, metavar='',
            help='number of bases to remove from beginning of read (integer)')
    parser.add_argument('-b', type=int, metavar='',
            help='final position to keep within a read (integer)')
    parser.add_argument('-e', type=int, metavar='',
            help='end-trim where entire window >= score (integer)')
    parser.add_argument('-w', type=int, metavar='',
            help='use with \'e\', size of window consisting of >= e (integer)')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    scallop_main()
