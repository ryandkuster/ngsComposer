import argparse
import gzip
import os
import shutil
import sys

from helpers.gzip_handling import gzip_test


def anemone_main():
    '''
    standalone, command line entry point to anemone using stdin
    '''
    in1 = args.r1
    out1 = os.path.basename(in1)
    try:
        in2 = args.r2 if args.r2 else False
        out2 = os.path.basename(in2)
    except TypeError:
        in2 = False
        out2 = False
    bcs_file = args.c
    mismatch = args.m if args.m else 0
    if args.o is None:
        proj = os.path.dirname(os.path.abspath(in1))
    elif os.path.exists(args.o) is True:
        proj = os.path.abspath(args.o)
    else:
        sys.exit('directory not found at ' + os.path.abspath(args.o))
    anemone_init(in1, in2, out1, out2, mismatch,
                 bcs_file, proj)


def anemone_comp(in1_ls, in2_ls, mismatch, bcs_dict, curr, in1):
    '''
    composer entry point to anemone
    '''
    if in1 in in2_ls:
        return

    out1 = os.path.basename(in1)
    compressed = gzip_test(in1)
    if compressed:
        base, _ = os.path.splitext(out1)
        if base in bcs_dict.keys():
            out1 = base
    else:
        if out1 + '.gz' in bcs_dict.keys():
            bcs_dict[out1] = bcs_dict.pop(out1 + '.gz')

    try:
        in2 = in2_ls[in1_ls.index(in1)]
        out2 = os.path.basename(in2)
    except (IndexError, ValueError) as e:
        in2 = False
        out2 = False

    '''
    the following copies files not found in index.txt
    '''
    try:
        bcs_file = bcs_dict[out1]
    except KeyError:
        shutil.copy(in1, curr)
        try:
            shutil.copy(in2, curr)
        except TypeError:
            pass
        return

    subdir = os.path.join(curr, os.path.basename(in1))
    os.mkdir(subdir)
    anemone_init(in1, in2, out1, out2, mismatch,
                 bcs_file, subdir)


def anemone_init(in1, in2, out1, out2, mismatch, bcs_file, proj):
    '''
    extract bcs from bcs_file and detect dual-indexing
    call paired-end or single-end anemone
    '''
    names_mx, R1_bcs, R2_bcs, dual_index = bc_reader(bcs_file)
    test = bc_test(R1_bcs, names_mx, True)
    if test:
        print('redundant R1 barcodes detected')
    test = bc_test(R2_bcs, names_mx, False)
    if test:
        print('redundant R2 barcodes detected')
    of1_ls = [open(proj + '/temp_unknown.' + out1, 'w')]
    if in2:
        of2_ls = [open(proj + '/temp_unknown.' + out2, 'w')]
    for i, item in enumerate(R1_bcs):
        of1_ls.append(open(proj + '/' + str(i) + '.' + out1, 'w'))
        if in2:
            of2_ls.append(open(proj + '/' + str(i) + '.' + out2, 'w'))

    of1_dict, of2_dict = {}, {}
    if in2:
        of1_ls, of2_ls = anemone_open(in1, in2, out1, out2, of1_ls, of2_ls,
                                      mismatch, R1_bcs, proj, True)
        if dual_index is True:
            of1_master, of2_master = dual_indexer(in1, in2, R2_bcs, proj,
                                                  of1_ls, of2_ls, mismatch)
            for file1, file2 in zip(of2_master, of1_master):
                of1_dict, of2_dict = rename_files(file1, file2, dual_index,
                                                  names_mx, proj, of1_dict,
                                                  of2_dict)
        else:
            for file1, file2 in zip(of1_ls, of2_ls):
                of1_dict, of2_dict = rename_files(file1.name, file2.name,
                                                  dual_index, names_mx,
                                                  proj, of1_dict, of2_dict)
    else:
        anemone_single_open(in1, out1, of1_ls, mismatch, R1_bcs, proj,
                            True)
        for file1 in of1_ls[1:]:
            of1_dict, of2_dict = rename_files(file1.name, None, dual_index,
                                              names_mx, proj, of1_dict,
                                              of2_dict)
    concatenate_files(proj, of1_dict, of2_dict)


def bc_reader(bcs_file):
    '''
    open user-defined bc file and extract forward and reverse
    bcs create matrix to pull corresponding sample ids that match
    bcs
    '''
    with open(bcs_file) as f:
        bcs_mx = [line.rstrip().split() for line in f]
    R2_bcs = {i: item for i, item in enumerate(bcs_mx[0])}
    R1_bcs = {i: item[0] for i, item in enumerate(bcs_mx[1:])}
    dual_index = False if len(R2_bcs) == 1 else True
    if (len(bcs_mx[1]) - len(bcs_mx[0])) != 1:
        sys.exit('format of barcodes file ' + bcs_file + ' is not recognized')
    names_mx = [item[1:] for item in bcs_mx[1:]]
    return names_mx, R1_bcs, R2_bcs, dual_index


def bc_test(bcs, names_mx, r1):
    '''
    test for redundant barcodes
    '''
    dupe_ls = []
    for k1, v1 in bcs.items():
        for k2, v2 in bcs.items():
            if v2.startswith(v1) and k1 != k2:
                dupe_ls.extend((k1, k2))
    if r1:
        for i in dupe_ls:
            for j in range(len(names_mx[i])):
                names_mx[i][j] = 'redundant'
    else:
        for i in range(len(names_mx)):
            for j in dupe_ls:
                names_mx[i][j] = 'redundant'
    return dupe_ls


def dual_indexer(in1, in2, R2_bcs, proj, of1_ls, of2_ls, mismatch):
    '''
    create of1/2_di_ls to direct output for final iteration
    create of1/2_masters to keep track of ALL output files in directory
    '''
    of1_master = [proj + '/unknown.' + os.path.basename(in2)]
    of2_master = [proj + '/unknown.' + os.path.basename(in1)]
    for file1, file2 in zip(of1_ls[1:], of2_ls[1:]):
        in2 = file1.name
        in1 = file2.name  # the R2 reads become in1
        out1 = os.path.basename(in1)
        out2 = os.path.basename(in2)
        of1_di_ls = [open(proj + '/temp_unknown.' + out1, 'w')]
        of2_di_ls = [open(proj + '/temp_unknown.' + out2, 'w')]
        of1_master.append(proj + '/unknown.' + out1)
        of2_master.append(proj + '/unknown.' + out2)
        for i, item in enumerate(R2_bcs):
            of1_di_ls.append(open(proj + '/' + str(i) + '.' + out1, 'w'))
            of2_di_ls.append(open(proj + '/' + str(i) + '.' + out2, 'w'))
            of1_master.append(proj + '/' + str(i) + '.' + out1)
            of2_master.append(proj + '/' + str(i) + '.' + out2)
        anemone_open(in1, in2, out1, out2, of1_di_ls, of2_di_ls, mismatch,
                     R2_bcs, proj, True)
        os.remove(in1)
        os.remove(in2)
    return of1_master, of2_master


def rename_files(file1, file2, dual_index, names_mx, proj, of1_dict,
                 of2_dict):
    '''
    using names_mx, rename files with possibility of duplication
    '''
    for j, element in enumerate(os.path.basename(file1).split('.')):
        if j == 0:
            if dual_index is True:
                y = element
            else:
                x = element
                y = 0
        if j == 1 and dual_index is True:
            x = element

    if x == 'unknown' or y == 'unknown':
        sample_id = 'unknown'
        rename1 = file1
        rename2 = file2 if file2 else None
    else:
        sample_id = names_mx[int(x)][int(y)]
        rename1 = proj + '/' + sample_id + '.' + \
            os.path.basename(file1)
        rename2 = proj + '/' + sample_id + '.' + \
            os.path.basename(file2) if file2 else None
        os.rename(file1, rename1)
        if file2:
            os.rename(file2, rename2)
    if sample_id in of1_dict:
        of1_dict[sample_id].append(rename1)
        if file2:
            of2_dict[sample_id].append(rename2)
    else:
        of1_dict[sample_id] = [rename1]
        if file2:
            of2_dict[sample_id] = [rename2]
    return of1_dict, of2_dict


def concatenate_files(proj, of1_dict, of2_dict):
    '''
    combine files with identical sample ids
    '''
    for sample_id in of1_dict.keys():
        with open(proj + '/' + sample_id + '.R1.fastq', 'w') as o1:
            for i in of1_dict[sample_id]:
                with open(i) as obj1:
                    shutil.copyfileobj(obj1, o1)
                os.remove(i)
    try:
        for sample_id in of2_dict.keys():
            with open(proj + '/' + sample_id + '.R2.fastq', 'w') as o2:
                for i in of2_dict[sample_id]:
                    with open(i) as obj2:
                        shutil.copyfileobj(obj2, o2)
                    os.remove(i)
    except AttributeError:
        pass


def anemone_open(in1, in2, out1, out2, of1_ls, of2_ls, mismatch, bcs, proj,
                 round_one):
    '''
    create IO file object based on gzipped status for pe data
    '''
    compressed = gzip_test(in1)
    if compressed:
        with gzip.open(in1, 'rt') as f1, gzip.open(in2, 'rt') as f2:
            of1_ls, of2_ls = anemone(f1, f2, out1, out2, of1_ls, of2_ls,
                                     mismatch, bcs, proj, round_one)
    else:
        with open(in1) as f1, open(in2) as f2:
            of1_ls, of2_ls = anemone(f1, f2, out1, out2, of1_ls, of2_ls,
                                     mismatch, bcs, proj, round_one)
    return of1_ls, of2_ls


def anemone_single_open(in1, out1, of1_ls, mismatch, bcs, proj, round_one):
    '''
    create IO file object based on gzipped status for se data
    '''
    compressed = gzip_test(in1)
    if compressed:
        with gzip.open(in1, 'rt') as f1:
            anemone_single(f1, out1, of1_ls, mismatch, bcs, proj,
                           round_one)
    else:
        with open(in1) as f1:
            anemone_single(f1, out1, of1_ls, mismatch, bcs, proj,
                           round_one)


def anemone(f1, f2, out1, out2, of1_ls, of2_ls, mismatch, bcs, proj,
            round_one):
    '''
    use active 'in1' file to demultiplex in a number of ways
    '''
    y, entry1, entry2 = 0, "", ""
    for line1, line2 in zip(f1, f2):
        y += 1
        if y == 2:
            z, output_prefix = exact_matches(line1, bcs) if round_one is\
                    True else mismatches(line1, bcs, mismatch)
        if y == 2 or y == 4:
            line1 = line1[z:]
        entry1 = entry1 + line1
        entry2 = entry2 + line2
        if y == 4:
            of1_ls[output_prefix].write(entry1)
            of2_ls[output_prefix].write(entry2)
            y, entry1, entry2 = 0, "", ""
    if round_one is True:
        if mismatch > 0:
            second_pass(out1, out2, of1_ls, of2_ls, mismatch, bcs, proj)
        else:
            for x in of1_ls:
                x.close()
            for x in of2_ls:
                x.close()
            of1_ls[0] = open(proj + '/unknown.' + out1, 'w')
            of2_ls[0] = open(proj + '/unknown.' + out2, 'w')
            os.rename(proj + '/temp_unknown.' + out1,
                      proj + '/unknown.' + out1)
            os.rename(proj + '/temp_unknown.' + out2,
                      proj + '/unknown.' + out2)
    else:
        os.remove(proj + '/temp_unknown.' + out1)
        os.remove(proj + '/temp_unknown.' + out2)
        for x in of1_ls:
            x.close()
        for x in of2_ls:
            x.close()
    return of1_ls, of2_ls


def anemone_single(f1, out1, of1_ls, mismatch, bcs, proj, round_one):
    '''
    use active 'in1' file to demultiplex in a number of ways
    '''
    y, entry1 = 0, ""
    for line1 in f1:
        y += 1
        if y == 2:
            z, output_prefix = exact_matches(line1, bcs) if round_one is\
                    True else mismatches(line1, bcs, mismatch)
        if y == 2 or y == 4:
            line1 = line1[z:]
        entry1 = entry1 + line1
        if y == 4:
            of1_ls[output_prefix].write(entry1)
            y, entry1 = 0, ""
    if round_one is True:
        if mismatch > 0:
            second_pass(out1, False, of1_ls, False, mismatch, bcs, proj)
        else:
            for x in of1_ls:
                x.close()
            os.rename(proj + '/temp_unknown.' + out1,
                      proj + '/unknown.' + out1)
    else:
        os.remove(proj + '/temp_unknown.' + out1)
        for x in of1_ls:
            x.close()


def exact_matches(line1, bcs):
    '''
    write to file only bcs with 100 percent match in first pass
    '''
    for file_prefix, x in bcs.items():
        if line1.startswith(x):
            output_prefix = file_prefix + 1
            z = len(x)
            break
        else:
            z = 0
            output_prefix = 0
    return z, output_prefix


def mismatches(line1, bcs, mismatch):
    '''
    if mismatch > 0 write to file bcs with leniency
    '''
    z, multi, output_prefix = 0, 0, 0
    for file_prefix, x in bcs.items():
        hamm = 0
        for j in range(len(x)):
            if x[j] != line1[j]:
                hamm = hamm + 1
                if hamm > mismatch:
                    break
        if hamm <= mismatch:
            output_prefix = file_prefix + 1
            z = len(x)
            multi += 1
        if multi > 1:
            z = 0
            output_prefix = 0
            break
    return z, output_prefix


def second_pass(out1, out2, of1_ls, of2_ls, mismatch, bcs, proj):
    '''
    after the first round of precision demultiplexing, attempt matches
    of unknown reads
    '''
    of1_ls[0].close()
    of1_ls[0] = open(proj + '/unknown.' + out1, 'w')
    if out2:
        of2_ls[0].close()
        of2_ls[0] = open(proj + '/unknown.' + out2, 'w')
    in1 = proj + '/temp_unknown.' + out1
    if out2:
        in2 = proj + '/temp_unknown.' + out2
        anemone_open(in1, in2, out1, out2, of1_ls,
                     of2_ls, mismatch, bcs, proj, False)
    else:
        anemone_single_open(in1, out1, of1_ls, mismatch,
                            bcs, proj, False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='demultiplex reads')
    parser.add_argument('-r1', type=str, required=True, metavar='',
            help='the full or relative path to R1 fastq file')
    parser.add_argument('-r2', type=str, metavar='',
            help='the full or relative path to R2 fastq file (optional)')
    parser.add_argument('-c', type=str, required=True, metavar='',
            help='the full or relative path to barcodes index file')
    parser.add_argument('-m', type=int, metavar='',
            help='mismatch value for barcode hamming distance (integer)')
    parser.add_argument('-o', type=str, metavar='',
            help='the full path to output directory (optional)')
    args = parser.parse_args()
    anemone_main()
