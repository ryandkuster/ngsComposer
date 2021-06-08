import argparse
import gzip
import os
import sys
from re import finditer
from helpers.gzip_handling import gzip_test
from anemone import anemone_comp, bc_reader


def get_args():
    parser = argparse.ArgumentParser(description='trim known adapter sequences')
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-r2', type=str, metavar='',
            help='the full or relative path to R1 or R2 fastq file')
    parser.add_argument('-a1', type=str, required=True, metavar='',
            help='the full or relative path to R1 adapter sequences file')
    parser.add_argument('-a2', type=str,  metavar='',
            help='the full or relative path to R2 adapter sequences file')
    parser.add_argument('-m1', type=str, nargs='+', metavar='',
            help='space separated list of motifs expected in R1 (no brackets)')
    parser.add_argument('-m2', type=str, nargs='+', metavar='',
            help='space separated list of motifs expected in R2 (no brackets)')
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
    return parser.parse_args()


class Pipeline_args:

    def __init__(self, match, min_l, tiny, R1_bases_ls, R2_bases_ls):
        self.k = 8
        self.r = 1
        self.match = match
        self.min_l = min_l
        self.tiny = tiny
        if R1_bases_ls:
            #self.R2_bases_ls = R1_bases_ls #TODO
            self.R2_bases_ls = None #TODO add adapter_validate
        else:
            self.R2_bases_ls = None
        if R2_bases_ls:
            #self.R1_bases_ls = R2_bases_ls #TODO
            self.R1_bases_ls = None #TODO add adapter_validate
        else:
            self.R1_bases_ls = None


def porifera_main():
    '''
    standalone, command line entry point to porifera using stdin
    '''
    args = get_args()
    args.in1 = args.r1
    args.in2 = args.r2 if args.r2 else None
    adapters1 = args.a1
    adapters2 = args.a2 if args.a2 else None

    with open(adapters1) as f:
        adapters_ls1 = [line.rstrip() for line in f]

    adapters_ls1 = reverse_comp(adapters_ls1)
    args.k = args.k if args.k else 8
    r = args.r if args.r else 1
    args.min_l = args.l if args.l else 0
    args.rounds = r * (len(max(adapters_ls1, key=len))//args.k)
    args.match = args.m if args.m else 12
    args.tiny_ls1 = set([i[:args.match + args.k] for i in adapters_ls1]) if args.t else None
    args.tiny = args.t if args.t else None
    args.subseqs1 = simple_seeker_non_contig(adapters_ls1, args.k)
    args.R1_bases_ls = args.m1 if args.m1 else None
    args.R2_bases_ls = args.m2 if args.m2 else None

    if args.o is None:
        proj = os.path.dirname(os.path.abspath(args.in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))

    if args.in2:
        args.se_1 = os.path.join(proj, 'se.adapted.' + os.path.basename(args.in1))
        args.pe_1 = os.path.join(proj, 'pe.adapted.' + os.path.basename(args.in1))
        args.se_2 = os.path.join(proj, 'se.adapted.' + os.path.basename(args.in2))
        args.pe_2 = os.path.join(proj, 'pe.adapted.' + os.path.basename(args.in2))
        if adapters2:
            with open(adapters2) as f:
                adapters_ls2 = [line.rstrip() for line in f]
            adapters_ls2 = reverse_comp(adapters_ls2)
            args.tiny_ls2 = set([i[:args.match + args.k] for i in adapters_ls2]) if args.t else None
            args.subseqs2 = simple_seeker_non_contig(adapters_ls2, args.k)
        else:
            args.subseqs2 = args.subseqs1
        porifera_open(args)
    else:
        args.se_1 = os.path.join(proj, 'adapted.' + os.path.basename(args.in1))
        porifera_single_open(args)


def porifera_comp(curr, in1_ls, in2_ls, adapters1, adapters2, bcs_dict,
                  search, match, min_l, tiny, R1_bases_ls, R2_bases_ls, in1):
    '''
    composer entry point to porifera
    '''
    args = Pipeline_args(match, min_l, tiny, R1_bases_ls, R2_bases_ls)
    args.in1 = in1
    with open(adapters1) as f:
        adapters_ls1 = [line.rstrip() for line in f]

    if bcs_dict:
        r1_barcodes, r2_barcodes = custom_adapters(bcs_dict, args.in1)
        subset_ls1 = [i for i in adapters_ls1 for j in r2_barcodes if j in i[-(search + len(j)):]]
        adapters_ls1 = subset_ls1[:] if len(subset_ls1) > 0 else adapters_ls1

    adapt1 = reverse_comp(adapters_ls1)
    args.tiny_ls1 = [i[:args.match + args.k] for i in adapt1] if args.tiny else []
    args.subseqs1 = simple_seeker_non_contig(adapt1, args.k)

    if adapters2:
        with open(adapters2) as f:
            adapters_ls2 = [line.rstrip() for line in f]
        if bcs_dict:
            subset_ls2 = [i for i in adapters_ls2 for j in r1_barcodes if j in i[-(search + len(j)):]]
            adapters_ls2 = subset_ls2[:] if len(subset_ls2) > 0 else adapters_ls2
        adapt2 = reverse_comp(adapters_ls2)
        args.subseqs2 = simple_seeker_non_contig(adapt2, args.k)
        args.tiny_ls2 = [i[:args.match + args.k] for i in adapt2] if args.tiny else []
        if in2_ls == []:
            adapters_ls1.extend(adapters_ls2)
            args.tiny_ls1.extend(args.tiny_ls2)
    else:
        args.subseqs2 = args.subseqs1
        args.tiny_ls2 = args.tiny_ls1

    args.rounds = args.r * (len(max(adapters_ls1, key=len))//args.k)
    try:
        args.in2 = in2_ls[in1_ls.index(args.in1)]
        args.pe_1 = os.path.join(curr, 'paired', os.path.basename(args.in1))
        args.se_1 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(args.in1))
        args.pe_2 = os.path.join(curr, 'paired', os.path.basename(args.in2))
        args.se_2 = os.path.join(curr, 'single', 'pe_lib', os.path.basename(args.in2))
        porifera_open(args)
    except (IndexError, ValueError) as e:
        args.se_1 = os.path.join(curr, 'single', 'se_lib', os.path.basename(args.in1))
        porifera_single_open(args)


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


def porifera_open(args):
    '''
    open paired end files for adapter detection
    '''
    compressed = gzip_test(args.in1)
    if compressed:
        pe_1, _ = os.path.splitext(args.pe_1)
        pe_2, _ = os.path.splitext(args.pe_2)
        se_1, _ = os.path.splitext(args.se_1)
        se_2, _ = os.path.splitext(args.se_2)
        with gzip.open(args.in1, 'rt') as f1, gzip.open(args.in2, 'rt') as f2,\
                open(pe_1, 'w') as pe_o1,\
                open(pe_2, 'w') as pe_o2,\
                open(se_1, 'w') as se_o1,\
                open(se_2, 'w') as se_o2:
            porifera(args, f1, f2, pe_o1, pe_o2, se_o1, se_o2)
    else:
        with open(args.in1, 'rt') as f1, open(args.in2, 'rt') as f2,\
                open(args.pe_1, 'w') as pe_o1,\
                open(args.pe_2, 'w') as pe_o2,\
                open(args.se_1, 'w') as se_o1,\
                open(args.se_2, 'w') as se_o2:
            porifera(args, f1, f2, pe_o1, pe_o2, se_o1, se_o2)


def porifera(args, f1, f2, pe_o1, pe_o2, se_o1, se_o2):
    y, entry1, entry2 = 0, "", ""

    for line1, line2 in zip(f1, f2):
        y += 1
        line1 = line1.rstrip()
        line2 = line2.rstrip()

        if y == 2:
            if args.R1_bases_ls:
                z1 = seed_search(line1, args.subseqs1, args.R1_bases_ls, args.k, args.match, args.rounds)
            else:
                z1 = adapter_search(line1, args.subseqs1, args.k, args.match, args.rounds)
            if args.tiny_ls1 and z1 == len(line1):
                z1 = tiny_handler(line1, args.tiny_ls1, args.tiny)

            if args.R2_bases_ls:
                z2 = seed_search(line2, args.subseqs2, args.R2_bases_ls, args.k, args.match, args.rounds)
            else:
                z2 = adapter_search(line2, args.subseqs2, args.k, args.match, args.rounds)
            if args.tiny_ls2 and z2 == len(line2):
                z2 = tiny_handler(line2, args.tiny_ls2, args.tiny)

        if y % 2 == 0:
            line1 = line1[:z1]
            line2 = line2[:z2]

        entry1 = entry1 + line1 + "\n"
        entry2 = entry2 + line2 + "\n"

        if y == 4:
            if z1 >= args.min_l and z2 >= args.min_l:
                pe_o1.write(entry1)
                pe_o2.write(entry2)
            elif z1 >= args.min_l:
                se_o1.write(entry1)
            elif z2 >= args.min_l:
                se_o2.write(entry2)
            else:
                pass
            y, entry1, entry2 = 0, "", ""


def porifera_single_open(args):
    '''
    open files for adapter detection
    '''
    compressed = gzip_test(args.in1)
    if compressed:
        se_1, _ = os.path.splitext(args.se_1)
        with gzip.open(args.in1, 'rt') as f, open(se_1, 'w') as o:
            porifera_single(args, f, o)
    else:
        with open(args.in1) as f, open(args.se_1, 'w') as o:
            porifera_single(args, f, o)


def porifera_single(args, f, o):
    y, entry = 0, ""
    for line in f:
        y += 1
        line = line.rstrip()
        if y == 2:
            if args.R1_bases_ls:
                z = seed_search(line, args.subseqs1, args.R1_bases_ls, args.k, args.match, args.rounds)
            else:
                z = adapter_search(line, args.subseqs1, args.k, args.match, args.rounds)
            if args.tiny_ls1 and z == len(line):
                z = tiny_handler(line, args.tiny_ls1, args.tiny)

        if y == 2 or y == 4:
            line = line[:z]
        entry = entry + line + "\n"

        if y == 4:
            if z >= args.min_l:
                o.write(entry)
            y, entry = 0, ""


def seed_search(ref, subseqs, re_ls, k, match, rounds):
    '''
    search reference (line) for kmers, store possible adapter start
    positions in focus dictionary
    '''
    focus, finder = {len(ref): []}, 0
    seed_ls = []
    for re_motif in re_ls:
        seed_ls.extend([i.start()+len(re_motif) for i in finditer(re_motif, ref)])

    for counter, (pos, motif_ls) in enumerate(subseqs.items()):
        for motif in motif_ls:
            if motif in ref:
                for i in finditer(motif, ref):
                    point = i.start() - pos
                    if point in focus:
                        focus[point] = focus[point] + [j for j in range(pos, (pos+k)) if j not in focus[point]]
                    elif point in seed_ls:
                        focus[point] = [j for j in range(pos, (pos+k))]

        finder = len(max(focus.values(), key=len))
        if finder >= match:
            return max(focus, key=lambda x:len(focus[x]))
        if finder < k and counter > rounds:
            break
    return len(ref)


def adapter_search(ref, subseqs, k, match, rounds):
    '''
    search reference (line) for kmers, store possible adapter start
    positions in focus list
    '''
    focus, finder = {len(ref): []}, 0
    for counter, (pos, motif_ls) in enumerate(subseqs.items()):
        for motif in motif_ls:
            if motif in ref:
                for j in finditer(motif, ref):
                    point = j.start() - pos
                    if point in focus:
                        focus[point].append(pos)
                    elif point >= 0:
                        focus[point] = [pos]
                finder = len(max(focus.values(), key=len))
                if finder > match-k:
                    test_z, z = match_analyzer(focus, k)
                    if test_z >= match:
                        return z
        if finder < 2 and counter > rounds:
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


def tiny_handler(line, tiny_ls, tiny):
    '''
    perform smith-waterman alignment on reduced substring of adapters
    '''
    for tiny_motif in tiny_ls:
        z = smith_waterman(tiny_motif, line[-len(tiny_motif):], tiny)
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
            match, mismatch, gap_ref, gap_query = 0, 0, 0, 0
            if ref_base == query_base:
                match = 3 + matrix[i][j]
            else:
                mismatch = -1 + matrix[i][j]
                gap_ref = -12 + matrix[i + 1][j]
                gap_query = -12 + matrix[i][j + 1]
            matrix[i + 1][j + 1] = max(match, mismatch, gap_ref, gap_query)

    exp_scores = [(i/idx)+idx/ref_len if idx >= tiny else 0 for idx, i in enumerate(matrix[-1])]

    if max(exp_scores) > 2.5:
        return exp_scores.index(max(exp_scores))
    return None


if __name__ == "__main__":
   porifera_main()
