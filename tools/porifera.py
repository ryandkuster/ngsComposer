import argparse
import gzip
import os
import sys


def porifera_main():
    '''
    standalone, command line entry point to porifera using stdin
    '''
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    in1 = args.r1
    dapt_path = args.a
    min_start = args.n
    mismatch = args.m
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    out1 = os.path.join(proj, 'adapted.' + os.path.basename(in1))
    porifera_open(in1, dapt_path, min_start, mismatch, out1)


def porifera_comp(curr, dapt_path, min_start, mismatch, in1):
    '''
    composer entry point to porifera
    '''
    if os.path.basename(os.path.dirname(in1)) == 'paired':
        curr = os.path.join(curr, 'paired')
    elif os.path.basename(os.path.dirname(in1)) == 'single':
        curr = os.path.join(curr, 'single')

    out1 = os.path.join(curr, os.path.basename(in1))
    porifera_open(in1, dapt_path, min_start, mismatch, out1)


def porifera_open(in1, dapt_path, min_start, mismatch, out1):
    '''
    open files for adapter detection
    '''
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    with open(dapt_path) as f:
        adapter = [line.rstrip() for line in f]
    for i, seq in enumerate(adapter):
        new = ''
        for j, base in enumerate(reversed(seq)):
            new += revc[base]
        adapter[i] = new
    try:
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            porifera(f, adapter, min_start, mismatch, o)
    except OSError:
        with open(in1) as f, open(out1, 'w') as o:
            porifera(f, adapter, min_start, mismatch, o)


def porifera(f, adapter, min_start, mismatch, o):
    y, entry = 0, ""
    for line in f:
        y += 1
        line = line.rstrip()
        if y == 2:
            z = adapt_align(line, adapter, min_start, mismatch)
        if y == 2 or y == 4:
            line = line[:z]
        entry = entry + line + "\n"
        if y == 4:
            o.write(entry)
            y, entry = 0, ""


def adapt_align(line, adapter, min_start, mismatch):
    pos_scores, z = [], len(line)
    for seq in adapter:
        for x in range(len(line) - min_start):
            hamm = 0
            for y in range(len(line[-(x + min_start + 1):])):
                try:
                    if line[-(x + min_start + 1):][y] != seq[y]:
                        hamm += 1
                except IndexError:
                    break
                if hamm > mismatch:
                    break
#            pos_scores.append(hamm) # TODO consider multimatch scenario
            if hamm <= mismatch:
                z = len(line) - min_start - x - 1
                return z
#    print(pos_scores)
    return z


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='trim known adapter sequences')
    parser.add_argument('-r1', type=str,
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-a', type=str,
            help='the full or relative path to adapter sequences file')
    parser.add_argument('-n', type=int,
            help='number of bases from end of read to begin adapter alignment')
    parser.add_argument('-m', type=int,
            help='mismatch value for adapter hamming distance (integer)')
    parser.add_argument('-o', type=str,
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    porifera_main()
