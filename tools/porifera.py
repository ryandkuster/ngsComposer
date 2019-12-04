import argparse
import gzip
import math
import os
import sys
from re import finditer
from helpers.gzip_handling import gzip_test


def porifera_main():
    '''
    standalone, command line entry point to porifera using stdin
    '''
    in1 = args.r1
    in2 = args.r2 if args.r2 else None
    a1 = args.a1
    a2 = args.a2 if args.a2 else None
    adapter1 = reverse_comp(a1)
    k = args.k if args.k else 8
    r = args.r if args.r else 1
    rounds = r * (len(max(adapter1, key=len))//k)
    match = args.m if args.m else 12
    subseqs1 = simple_seeker_non_contig(adapter1, k)
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    out1 = os.path.join(proj, 'adapted.' + os.path.basename(in1))
    if in2:
        out2 = os.path.join(proj, 'adapted.' + os.path.basename(in2))
        if a2:
            adapter2 = reverse_comp(a2)
            subseqs2 = simple_seeker_non_contig(adapter2, k)
        else:
            subseqs2 = subseqs1
        porifera_open(in1, in2, subseqs1, subseqs2, out1, out2, k, rounds, match)
    else:
        porifera_single_open(in1, subseqs1, out1, k, rounds, match)


def porifera_comp(curr, dapt_path, match, in1):
    '''
    composer entry point to porifera
    '''
    if os.path.basename(os.path.dirname(in1)) == 'paired':
        curr = os.path.join(curr, 'paired')
    elif os.path.basename(os.path.dirname(in1)) == 'single':
        curr = os.path.join(curr, 'single')
    #TODO add se vs. pe open
    k = 8
    rounds = 1
    adapter = reverse_comp(dapt_path)
    subseqs = simple_seeker_non_contig(adapter, k)
    out1 = os.path.join(curr, os.path.basename(in1))
    porifera_single_open(in1, subseqs, out1, k, rounds, match)


def reverse_comp(dapt_path):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    with open(dapt_path) as f:
        adapter = [line.rstrip() for line in f]
    for i, seq in enumerate(adapter):
        new = ''
        for j, base in enumerate(reversed(seq)):
            new += revc[base]
        adapter[i] = new
    return adapter


def simple_seeker_non_contig(adapter, k):
    '''
    create subseqs dictionary of adapters size k
    each key is the distance to the beginning of the sequence
    keys are ordered such that if k = 3, all index 0 would be
    followed by all index 3 etc., then back to index 1, then 4... 
    '''
    subseqs = {}
    subs_mx = [[] for i in range(len(max(adapter, key=len)))]
    for seq in adapter:
        for i, base in enumerate(seq):
            if len(seq[i:i+k]) == k:
                subs_mx[i].append(seq[i:i+k])
    for i in range(k):
        for j, motif_ls in enumerate(subs_mx):
            if ((i-j)%k) == 0:
                subseqs[j] = motif_ls
    subseqs = no_empty_lists(subseqs)
    return subseqs


def no_empty_lists(subseqs):
    '''
    clean up subseqs dictionary of empty values
    '''
    rm_seqs = []
    for key, val in subseqs.items():
        if len(val) == 0:
            rm_seqs.append(key)
        subseqs[key] = list(set(val))
    for i in rm_seqs:
        del subseqs[i]
    return subseqs


def porifera_open(in1, in2, subseqs1, subseqs2, out1, out2, k, rounds, match):
    '''
    open paired end files for adapter detection
    '''
    compressed = gzip_test(in1)
    if compressed:
        with gzip.open(in1, 'rt') as f1, gzip.open(in2, 'rt') as f2,\
                open(out1, 'w') as o1, open(out2, 'w') as o2:
            porifera(f1, f2, subseqs1, subseqs2, o1, o2, k, rounds, match)
    else:
        with open(in1, 'rt') as f1, open(in2, 'rt') as f2,\
                open(out1, 'w') as o1, open(out2, 'w') as o2:
            porifera(f1, f2, subseqs1, subseqs2, o1, o2, k, rounds, match)


def porifera(f1, f2, subseqs1, subseqs2, o1, o2, k, rounds, match):
    y, entry1, entry2 = 0, "", ""
    for line1, line2 in zip(f1, f2):
        y += 1
        line1 = line1.rstrip()
        line2 = line2.rstrip()
        if y == 2:
            pass
            z1 = compromiser(line1, subseqs1, k, rounds, match)
            z2 = compromiser(line2, subseqs2, k, rounds, match)
        if y == 2 or y == 4:
            line1 = line1[:z1]
            line2 = line2[:z2]
        entry1 = entry1 + line1 + "\n"
        entry2 = entry2 + line2 + "\n"
        if y == 4:
            o1.write(entry1)
            o2.write(entry2)
            y, entry1, entry2 = 0, "", ""


def porifera_single_open(in1, subseqs, out1, k, rounds, match):
    '''
    open files for adapter detection
    '''
    compressed = gzip_test(in1)
    if compressed:
        with gzip.open(in1, 'rt') as f, open(out1, 'w') as o:
            porifera_single(f, subseqs, o, k, rounds, match)
    else:
        with open(in1) as f, open(out1, 'w') as o:
            porifera_single(f, subseqs, o, k, rounds, match)


def porifera_single(f, subseqs, o, k, rounds, match):
    y, entry = 0, ""
    for line in f:
        y += 1
        line = line.rstrip()
        if y == 2:
            z = compromiser(line, subseqs, k, rounds, match)
        if y == 2 or y == 4:
            line = line[:z]
        entry = entry + line + "\n"
        if y == 4:
            o.write(entry)
            y, entry = 0, ""


def compromiser(ref, subseqs, k, rounds, match):
    '''
    search reference (line) for kmers, store possible adapter start
    positions in focus list
    '''
    focus, finder = {}, 0
    for counter, (pos, motif_ls) in enumerate(subseqs.items()):
        for motif in motif_ls:
            if motif in ref:
                for j in finditer(motif, ref):
                    point = j.start() - pos
                    if point in focus:
                        focus[point].append(pos)
                    else:
                        focus[point] = [pos]
                finder = len(max(focus.values(), key=len))
                if finder > 4: #TODO allow for different value
                    test_z, z = match_analyzer(focus, k)
                    if test_z >= match:
                        return max(0, z)
        if finder < 2 and counter > rounds: #TODO allow for different value
            break
    return len(ref)


def match_analyzer(focus, k):
    '''
    check focus list to see if estimated adapter position
    is found a desired ratio of times
    '''
    z = max(focus, key=lambda x:len(focus[x]))
    positions = focus[z]
    pos_list = []
    for i in positions:
        pos_list.extend([j for j in range(i, (i+k))])
    pos_list = set(pos_list)
    return len(pos_list), z


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='trim known adapter sequences')
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-r2', type=str, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-a1', type=str, required=True, metavar='',
            help='the full or relative path to R1 adapter sequences file')
    parser.add_argument('-a2', type=str,  metavar='',
            help='the full or relative path to R2 adapter sequences file')
    parser.add_argument('-k', type=int, metavar='',
            help='k-mer size (for testing)')
    parser.add_argument('-r', type=int, metavar='',
            help='k-mer detection (for testing)')
    parser.add_argument('-m', type=int, metavar='',
            help='number of matching bases to detect adapter')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    porifera_main()
