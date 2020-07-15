import argparse
import gzip
import os
import sys
from re import finditer
from helpers.gzip_handling import gzip_test
from anemone import anemone_comp, bc_reader


def porifera_main():
    '''
    standalone, command line entry point to porifera using stdin
    '''
    in1 = args.r1
    in2 = args.r2 if args.r2 else None
    adapters1 = args.a1
    adapters2 = args.a2 if args.a2 else None
    with open(adapters1) as f:
        adapters_ls1 = [line.rstrip() for line in f]
    adapters_ls1 = reverse_comp(adapters_ls1)
    k = args.k if args.k else 8
    r = args.r if args.r else 1
    min_l = args.l if args.l else 0
    rounds = r * (len(max(adapters_ls1, key=len))//k)
    match = args.m if args.m else 12
    tiny_ls1 = set([i[:match] for i in adapters_ls1]) if args.t else None
    tiny = args.t if args.t else None
    subseqs1 = simple_seeker_non_contig(adapters_ls1, k)
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    if in2:
        se_1 = os.path.join(proj, 'se.adapted.' + os.path.basename(in1))
        pe_1 = os.path.join(proj, 'pe.adapted.' + os.path.basename(in1))
        se_2 = os.path.join(proj, 'se.adapted.' + os.path.basename(in2))
        pe_2 = os.path.join(proj, 'pe.adapted.' + os.path.basename(in2))
        if adapters2:
            with open(adapters2) as f:
                adapters_ls2 = [line.rstrip() for line in f]
            adapters_ls2 = reverse_comp(adapters_ls2)
            tiny_ls2 = set([i[:match] for i in adapters_ls2]) if args.t else None
            subseqs2 = simple_seeker_non_contig(adapters_ls2, k)
        else:
            subseqs2 = subseqs1
        porifera_open(in1, in2, subseqs1, subseqs2, se_1, pe_1, se_2, pe_2, k,
                      rounds, match, min_l, tiny_ls1, tiny_ls2, tiny)
    else:
        se_1 = os.path.join(proj, 'adapted.' + os.path.basename(in1))
        porifera_single_open(in1, subseqs1, se_1, k, rounds, match, min_l, tiny_ls1, tiny)


def porifera_comp(curr, in1_ls, in2_ls, adapters1, adapters2, bcs_dict,
                  search, match, min_l, tiny, in1):
    '''
    composer entry point to porifera
    '''
    k = 8
    r = 1

    with open(adapters1) as f:
        adapters_ls1 = [line.rstrip() for line in f]

    if bcs_dict:
        r1_barcodes, r2_barcodes = custom_adapters(bcs_dict, in1)
        subset_ls1 = [i for i in adapters_ls1 for j in r2_barcodes if j in i[-(search + len(j)):]]
        adapters_ls1 = subset_ls1[:] if len(subset_ls1) > 0 else adapters_ls1

    adapt1 = reverse_comp(adapters_ls1)
    tiny_ls1 = [i[:match + k] for i in adapt1] if tiny else []
    subseqs1 = simple_seeker_non_contig(adapt1, k)

    if adapters2:
        with open(adapters2) as f:
            adapters_ls2 = [line.rstrip() for line in f]
        if bcs_dict:
            subset_ls2 = [i for i in adapters_ls2 for j in r1_barcodes if j in i[-(search + len(j)):]]
            adapters_ls2 = subset_ls2[:] if len(subset_ls2) > 0 else adapters_ls2
        adapt2 = reverse_comp(adapters_ls2)
        subseqs2 = simple_seeker_non_contig(adapt2, k)
        tiny_ls2 = [i[:match + k] for i in adapt2] if tiny else []
        if in2_ls == []:
            adapters_ls1.extend(adapters_ls2)
            tiny_ls1.extend(tiny_ls2)
    else:
        subseqs2 = subseqs1
        tiny_ls2 = tiny_ls1


    rounds = r * (len(max(adapters_ls1, key=len))//k)
    try:
        in2 = in2_ls[in1_ls.index(in1)]
        pe_1 = os.path.join(curr, 'paired', os.path.basename(in1))
        se_1 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in1))
        pe_2 = os.path.join(curr, 'paired', os.path.basename(in2))
        se_2 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(in2))
        porifera_open(in1, in2, subseqs1, subseqs2, se_1, pe_1, se_2, pe_2, k,
                      rounds, match, min_l, tiny_ls1, tiny_ls2, tiny)
    except (IndexError, ValueError) as e:
        se_1 = os.path.join(curr, 'single', 'se_lib', os.path.basename(in1))
        porifera_single_open(in1, subseqs1, se_1, k, rounds, match, min_l,
                             tiny_ls1, tiny)


def custom_adapters(bcs_dict, in1):
    '''
    reconstruct adapters by combining universal adapter with barcodes
    and RE sites if included in pipeline
    '''
    name = os.path.basename(in1)[:-9]
    for i in bcs_dict.values():
        r1_barcodes = []
        r2_barcodes = []
        names_mx, r1_bcs, r2_bcs, dual_index = bc_reader(i)
        for bc1, j in enumerate(names_mx):
            for bc2, k in enumerate(j):
                if k == name:
                    r1_barcodes.append(r1_bcs[bc1])
                    r2_barcodes.append(r2_bcs[bc2])
    return set(r1_barcodes), set(r2_barcodes)


def reverse_comp(adapters_ls):
    revc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for i, seq in enumerate(adapters_ls):
        new = ''
        for j, base in enumerate(reversed(seq)):
            new += revc[base]
        adapters_ls[i] = new
    return adapters_ls


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


def porifera_open(in1, in2, subseqs1, subseqs2, se_1, pe_1, se_2, pe_2, k,
                  rounds, match, min_l, tiny_ls1, tiny_ls2, tiny):
    '''
    open paired end files for adapter detection
    '''
    compressed = gzip_test(in1)
    if compressed:
        pe_1, _ = os.path.splitext(pe_1)
        pe_2, _ = os.path.splitext(pe_2)
        se_1, _ = os.path.splitext(se_1)
        se_2, _ = os.path.splitext(se_2)
        with gzip.open(in1, 'rt') as f1, gzip.open(in2, 'rt') as f2,\
                open(pe_1, 'w') as pe_o1,\
                open(pe_2, 'w') as pe_o2,\
                open(se_1, 'w') as se_o1,\
                open(se_2, 'w') as se_o2:
            porifera(f1, f2, subseqs1, subseqs2, pe_o1, pe_o2, se_o1, se_o2, k,
                     rounds, match, min_l, tiny_ls1, tiny_ls2, tiny)
    else:
        with open(in1, 'rt') as f1, open(in2, 'rt') as f2,\
                open(pe_1, 'w') as pe_o1,\
                open(pe_2, 'w') as pe_o2,\
                open(se_1, 'w') as se_o1,\
                open(se_2, 'w') as se_o2:
            porifera(f1, f2, subseqs1, subseqs2, pe_o1, pe_o2, se_o1, se_o2, k,
                     rounds, match, min_l, tiny_ls1, tiny_ls2, tiny)


def porifera(f1, f2, subseqs1, subseqs2, pe_o1, pe_o2, se_o1, se_o2, k,
             rounds, match, min_l, tiny_ls1, tiny_ls2, tiny):
    y, entry1, entry2 = 0, "", ""

    for line1, line2 in zip(f1, f2):
        y += 1
        line1 = line1.rstrip()
        line2 = line2.rstrip()
        if y == 2:
            z1 = compromiser(line1, subseqs1, k, rounds, match)
            if tiny_ls1 and z1 == len(line1):
                z1 = tiny_handler(line1, match, tiny_ls1, tiny)
            z2 = compromiser(line2, subseqs2, k, rounds, match)
            if tiny_ls2 and z2 == len(line2):
                z2 = tiny_handler(line2, match, tiny_ls2, tiny)
        if y == 2 or y == 4:
            line1 = line1[:z1]
            line2 = line2[:z2]
        entry1 = entry1 + line1 + "\n"
        entry2 = entry2 + line2 + "\n"
        if y == 4:
            if z1 >= min_l and z2 >= min_l:
                pe_o1.write(entry1)
                pe_o2.write(entry2)
            elif z1 >= min_l:
                se_o1.write(entry1)
            elif z2 >= min_l:
                se_o2.write(entry2)
            else:
                pass
            y, entry1, entry2 = 0, "", ""


def porifera_single_open(in1, subseqs, se_1, k, rounds, match, min_l, tiny_ls, tiny):
    '''
    open files for adapter detection
    '''
    compressed = gzip_test(in1)
    if compressed:
        se_1, _ = os.path.splitext(se_1)
        with gzip.open(in1, 'rt') as f, open(se_1, 'w') as o:
            porifera_single(f, subseqs, o, k, rounds, match, min_l, tiny_ls, tiny)
    else:
        with open(in1) as f, open(se_1, 'w') as o:
            porifera_single(f, subseqs, o, k, rounds, match, min_l, tiny_ls, tiny)


def porifera_single(f, subseqs, o, k, rounds, match, min_l, tiny_ls, tiny):
    y, entry = 0, ""
    for line in f:
        y += 1
        line = line.rstrip()
        if y == 2:
            z = compromiser(line, subseqs, k, rounds, match)
            if tiny_ls and z == len(line):
                z = tiny_handler(line, match, tiny_ls, tiny)
        if y == 2 or y == 4:
            line = line[:z]
        entry = entry + line + "\n"
        if y == 4:
            if z >= min_l:
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


def tiny_handler(line, match, tiny_ls, tiny):
    '''
    perform smith-waterman alignment on reduced substring of adapters
    '''
    for tiny_motif in tiny_ls:
        z = smith_waterman(tiny_motif, line[-match:], tiny)
        if z is not None:
            return len(line) - z
    return len(line)


def smith_waterman(query, ref, tiny):
    '''
    use smith-waterman dynamic programming matrix to search for substring of
    adapter with minimal overlap in the 3' end of read
    (evaluates scores in the final row of matrix, assuming substring will be
    exhausted)
    '''
    ref_len = len(ref)
    matrix = [[0 for j in range(len(query) + 1)] for i in range(len(ref) + 1)]

    for i, ref_base in enumerate(ref):
        for j, query_base in enumerate(query):
            # use max value
            match, mismatch, gap_ref, gap_query = 0, 0, 0, 0
            if ref_base == query_base:
                match = 3 + matrix[i][j]
            else:
                mismatch = -3 + matrix[i][j]
                gap_ref = -2 + matrix[i + 1][j]
                gap_query = -2 + matrix[i][j + 1]
            matrix[i + 1][j + 1] = max(match, mismatch, gap_ref, gap_query)

    exp_scores = [(i/idx)+idx/ref_len if idx >= tiny else 0 for idx, i in enumerate(matrix[-1])]

    if max(exp_scores) > 2.5:
        x_tr = exp_scores.index(max(exp_scores))
        y_tr = ref_len-2 if x_tr == 1 else traceback(matrix, x_tr, ref_len)
        return ref_len - y_tr - 1
    return None

def traceback(matrix, x, y):
    while True:
        now = matrix[y][x]
        prev = [matrix[y - 1][x - 1], matrix[y][x - 1], matrix[y - 1][x]]
        if now + 3 == prev[0]:
            x, y = x - 1, y - 1
        elif now + 2 == prev[1]:
            x, y = x - 1, y
        elif now + 2 == prev[2]:
            x, y = x, y - 1
        else:
            x, y = x - 1, y - 1
        if now == 0:
            return y


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
    parser.add_argument('-l', type=int, metavar='',
            help='minimum read length to keep (integer)')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    parser.add_argument('-t', type=int, metavar='',
            help='\"tiny" mode, determine smallest length to call an adapter\
                 (optional; takes considerably longer)')
    args = parser.parse_args()
    porifera_main()
