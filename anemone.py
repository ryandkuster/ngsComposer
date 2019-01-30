import sys
import os
import shutil

def anemone_main():
    '''
    standalone, command line entry point to anemone using stdin
    '''
    mismatch = int(sys.argv[1]) # number of mismatches allowed
    bcs_file = sys.argv[2] # bcs file
    input1 = sys.argv[3] # R1 reads
    output1 = os.path.basename(input1)
    try:
        input2 = sys.argv[4] # R2 reads
        output2 = os.path.basename(input2)
    except:
        input2 = False
        output2 = False
    chunk = 3000000 # how many reads to process before writing
    proj_dir = os.getcwd()
    anemone_init(input1, input2, output1, output2, mismatch, chunk,
                 bcs_file, proj_dir)


def anemone_pipeline(
        input1_ls, input2_ls, mismatch, bcs_dict, proj_dir, input1):
    '''
    composer entry point to anemone
    '''
    proj_dir = proj_dir + '/' + os.path.basename(input1)
    os.mkdir(proj_dir)
    for k, v in bcs_dict.items():
        if k == os.path.basename(input1):
            bcs_file = v
    try:
        input2 = input2_ls[input1_ls.index(input1)]
        output2 = os.path.basename(input2)
    except:
        input2 = False
        output2 = False
    output1 = os.path.basename(input1)
    chunk = 3000000
    anemone_init(input1, input2, output1, output2, mismatch, chunk,
                 bcs_file, proj_dir)


def anemone_init(
        input1, input2, output1, output2, mismatch, chunk, bcs_file, proj_dir):
    '''
    extract bcs from bcs_file and detect dual-indexing
    call paired-end or single-end anemone
    '''
    bcs_matrix, R1_bcs, R2_bcs, dual_index = bc_reader(bcs_file)
    of1_ls = [open(proj_dir + '/temp_unknown.' + output1, 'w')]
    try:
        of2_ls = [open(proj_dir + '/temp_unknown.' + output2, 'w')]
    except:
        of2_ls = []
    for i, item in enumerate(R1_bcs):
        of1_ls.append(open(proj_dir + '/' + str(i) + '.' + output1, 'w'))
        try:
            of2_ls.append(open(proj_dir + '/' + str(i) + '.' + output2, 'w'))
        except:
            pass
    if input2:
        anemone(
            input1, input2, output1, output2, of1_ls,
            of2_ls, mismatch, 3000000, R1_bcs, proj_dir, True)
        if dual_index == True:
            of1_master = [proj_dir + '/unknown.' + os.path.basename(input2)]
            of2_master = [proj_dir + '/unknown.' + os.path.basename(input1)]
            of1_dict, of2_dict = dual_indexer(
                bcs_matrix, R2_bcs, proj_dir,
                of1_ls, of2_ls, mismatch,
                of1_master, of2_master)
        elif dual_index == False:
            of1_dict, of2_dict = {}, {}
            for file1, file2 in zip(of1_ls, of2_ls):
                file1 = file1.name
                file2 = file2.name
                for j, element in enumerate(os.path.basename(file1).split('.')):
                    if j == 0:
                        x = element
                if x == 'unknown':
                    sample_id = 'unknown'
                    rename1 = file1
                    rename2 = file2
                else:
                    sample_id = bcs_matrix[0][int(x)]
                    rename1 = proj_dir + '/' + sample_id + '.' + os.path.basename(file1)
                    rename2 = proj_dir + '/' + sample_id + '.' + os.path.basename(file2)
                    os.rename(file1, rename1)
                    os.rename(file2, rename2)
                if sample_id in of1_dict:
                    of1_dict[sample_id].append(rename1)
                    of2_dict[sample_id].append(rename2)
                else:
                    of1_dict[sample_id] = [rename1]
                    of2_dict[sample_id] = [rename2]
    elif input2 == False:
        of1_dict, of2_dict = {}, None
        anemone_single(
            input1, output1, of1_ls, mismatch,
            chunk, R1_bcs, proj_dir, True)
        for file1 in of1_ls[1:]:
            file1 = file1.name
            for j, element in enumerate(os.path.basename(file1).split('.')):
                if j == 0:
                    x = element
            if x == 'unknown':
                sample_id = 'unknown'
                rename1 = file1
                rename2 = file2
            else:
                sample_id = bcs_matrix[0][int(x)]
                rename1 = proj_dir + '/' + sample_id + '.' + os.path.basename(file1)
                os.rename(file1, rename1)
            if sample_id in of1_dict:
                of1_dict[sample_id].append(rename1)
            else:
                of1_dict[sample_id] = [rename1]
    concatenate_files(proj_dir, of1_dict, of2_dict)


def bc_reader(bcs_file):
    '''
    open user-defined bc file and extract forward and reverse
    bcs create matrix to pull corresponding sample ids that match
    bcs
    '''
    bcs_matrix = []
    with open(bcs_file) as f:
        line = f.readline()
        R1_bcs, R2_bcs = {}, {}
        for i, item in enumerate(line.split()):
            R2_bcs[item] = i
        if len(R2_bcs) == 1:
            dual_index = False
        else:
            dual_index = True
        bcs_matrix = [[0] * 1 for i in range(len(R2_bcs))]
        for i, line in enumerate(f):
            for j, item in enumerate(line.split()):
                if j == 0:
                    R1_bcs[item] = i
                if j > 0:
                    bcs_matrix[j - 1].append(item)
        for i, line in enumerate(bcs_matrix):
            del line[0]
    return bcs_matrix, R1_bcs, R2_bcs, dual_index


def dual_indexer(
        bcs_matrix, R2_bcs, proj_dir, of1_ls, of2_ls, mismatch, of1_master,
        of2_master):
    '''
    create of1/2_final_lss to direct output for final iteration
    create of1/2_masters to keep track of ALL output files in directory
    create of1/2_dicts finds common sample ids associated with outputs
    '''
    of1_dict, of2_dict = {}, {}
    for file1, file2 in zip(of1_ls[1:], of2_ls[1:]):
        input2 = file1.name
        input1 = file2.name # the of2 results become input1...
        output1 = os.path.basename(input1)
        output2 = os.path.basename(input2)
        of1_final_ls = [open(proj_dir + '/temp_unknown.' + output1, 'w')]
        of2_final_ls = [open(proj_dir + '/temp_unknown.' + output2, 'w')]
        of1_master.append(proj_dir + '/unknown.' + output1)
        of2_master.append(proj_dir + '/unknown.' + output2)
        for i, item in enumerate(R2_bcs):
            of1_final_ls.append(open(proj_dir + '/' + str(i) + '.' + output1, 'w'))
            of2_final_ls.append(open(proj_dir + '/' + str(i) + '.' + output2, 'w'))
            of1_master.append(proj_dir + '/' + str(i) + '.' + output1)
            of2_master.append(proj_dir + '/' + str(i) + '.' + output2)
        anemone(
            input1, input2, output1, output2, of1_final_ls,
            of2_final_ls, mismatch, 3000000, R2_bcs, proj_dir,
            True)
        os.remove(input1)
        os.remove(input2)

    for i, (file1, file2) in enumerate(zip(of2_master, of1_master)):
        for j, element in enumerate(os.path.basename(file1).split('.')):
            if j == 0:
                x = element
            if j == 1:
                y = element
        if x == 'unknown':
            sample_id = 'unknown'
            rename1 = file1
            rename2 = file2
        else:
            sample_id = bcs_matrix[int(x)][int(y)]
            rename1 = proj_dir + '/' + sample_id + '.' + os.path.basename(file1)
            rename2 = proj_dir + '/' + sample_id + '.' + os.path.basename(file2)
            os.rename(file1, rename1)
            os.rename(file2, rename2)
            of1_master[i] = rename1
            of2_master[i] = rename2
        if sample_id in of1_dict:
            of1_dict[sample_id].append(rename1)
            of2_dict[sample_id].append(rename2)
        else:
            of1_dict[sample_id] = [rename1]
            of2_dict[sample_id] = [rename2]
    return of1_dict, of2_dict


def concatenate_files(proj_dir, of1_dict, of2_dict):
    '''
    combine files with identical sample ids
    '''
    for sample_id in of1_dict.keys():
        with open(proj_dir + '/' + sample_id + '.R1.fastq', 'w') as o1:
            for i in of1_dict[sample_id]:
                with open (i) as obj1:
                    shutil.copyfileobj(obj1, o1)
                os.remove(i)
    try:
        for sample_id in of2_dict.keys():
            with open(proj_dir + '/' + sample_id + '.R2.fastq', 'w') as o2:
                for i in of2_dict[sample_id]:
                    with open (i) as obj2:
                        shutil.copyfileobj(obj2, o2)
                    os.remove(i)
    except AttributeError:
        pass


def anemone(
        input1, input2, output1, output2, of1_ls, of2_ls, mismatch, chunk, bcs,
        proj_dir, round_one):
    '''
    use active 'input1' file to demultiplex in a number of ways
    '''
    row_len = len(bcs) + 1
    matrix_one = matrix_maker(row_len)
    matrix_two = matrix_maker(row_len)
    i, y, entry1, entry2 = 0, 0, "", ""
    with open(input1) as f1, open(input2) as f2:
        for line1, line2 in zip(f1, f2):
            i += 1
            y += 1
            if y == 2 and round_one == True:
                z, output_prefix = exact_matches(line1, bcs)
            if y == 2 and round_one == False:
                z, output_prefix = mismatches(line1, bcs, mismatch)
            if y == 2 or y == 4:
                line1 = line1[z:]
            entry1 = entry1 + line1
            entry2 = entry2 + line2
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                matrix_two[output_prefix].append(entry2)
                y, entry1, entry2 = 0, "", ""
            if i == chunk:
                unload(matrix_one, of1_ls)
                unload(matrix_two, of2_ls)
                i = 0
                matrix_one = matrix_maker(row_len)
                matrix_two = matrix_maker(row_len)
        unload(matrix_one, of1_ls)
        unload(matrix_two, of2_ls)
    if round_one == True:
        if mismatch > 0:
            second_pass(
                output1, output2, of1_ls, of2_ls,
                mismatch, chunk, bcs, proj_dir)
        else:
            for x in of1_ls:
                x.close()
            for x in of2_ls:
                x.close()
            os.rename(proj_dir + '/temp_unknown.' + output1,
                      proj_dir + '/unknown.' + output1)
            os.rename(proj_dir + '/temp_unknown.' + output2,
                      proj_dir + '/unknown.' + output2)
    else:
        os.remove(proj_dir + '/temp_unknown.' + output1)
        os.remove(proj_dir + '/temp_unknown.' + output2)
        for x in of1_ls:
            x.close()
        for x in of2_ls:
            x.close()


def anemone_single(
        input1, output1, of1_ls, mismatch, chunk, bcs, proj_dir, round_one):
    '''
    use active 'input1' file to demultiplex in a number of ways
    '''
    row_len = len(bcs) + 1
    matrix_one = matrix_maker(row_len)
    i, y, entry1 = 0, 0, ""
    with open(input1) as f1:
        for line1 in f1:
            i += 1
            y += 1
            if y == 2 and round_one == True:
                z, output_prefix = exact_matches(line1, bcs)
            if y == 2 and round_one == False:
                z, output_prefix = mismatches(line1, bcs, mismatch)
            if y == 2 or y == 4:
                line1 = line1[z:]
            entry1 = entry1 + line1
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                y, entry1 = 0, ""
            if i == chunk:
                unload(matrix_one, of1_ls)
                i = 0
                matrix_one = matrix_maker(row_len)
        unload(matrix_one, of1_ls)
    if round_one == True:
        if mismatch > 0:
            second_pass(
                output1, False, of1_ls, False, mismatch,
                chunk, bcs, proj_dir)
        else:
            for x in of1_ls:
                x.close()
            os.rename(proj_dir + '/temp_unknown.' + output1,
                      proj_dir + '/unknown.' + output1)
    else:
        os.remove(proj_dir + '/temp_unknown.' + output1)
        for x in of1_ls:
            x.close()


def matrix_maker(row_len):
    '''
    create empty matrix
    '''
    matrix = [[0] * 1 for i in range(row_len)]
    for i in range(row_len):
        matrix[i][0] = ''
    return matrix


def exact_matches(line1, bcs):
    '''
    write to file only bcs with 100 percent match in first pass
    '''
    for x, file_prefix in bcs.items():
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
    for x, file_prefix in bcs.items():
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


def unload(matrix, ofs):
    '''
    write matrix out to file once chunk value is achieved
    '''
    for x, fo in enumerate(ofs):
        fo.write(str(''.join(matrix[x])))


def second_pass(
        output1, output2, of1_ls, of2_ls, mismatch, chunk, bcs, proj_dir):
    '''
    after the first round of precision demultiplexing, attempt matches
    of unknown reads
    '''
    of1_ls[0].close()
    of1_ls[0] = open(proj_dir + '/unknown.' + output1, 'w')
    if output2:
        of2_ls[0].close()
        of2_ls[0] = open(proj_dir + '/unknown.' + output2, 'w')
    input1 = proj_dir + '/temp_unknown.' + output1
    if output2:
        input2 = proj_dir + '/temp_unknown.' + output2
        anemone(input1, input2, output1, output2, of1_ls,
                of2_ls, mismatch, chunk, bcs, proj_dir, False)
    else:
        anemone_single(input1, output1, of1_ls, mismatch,
                       chunk, bcs, proj_dir, False)

if __name__ == '__main__':
    anemone_main()
