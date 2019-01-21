import sys
import os


def anemone_main():
    '''
    standalone, command line entry point to anemone using stdin
    '''
    mismatch = int(sys.argv[1]) # number of mismatches allowed
    barcodes_file = sys.argv[2] # barcodes file
    input1 = sys.argv[3] # R1 reads
    output1 = os.path.basename(input1)
    try:
        input2 = sys.argv[4] # R2 reads
        output2 = os.path.basename(input2)
    except:
        input2 = False
        output2 = False
    chunk = 3000000 # how many reads to process before writing
    project_dir = os.getcwd()
    anemone_init(input1, input2, output1, output2, mismatch, chunk,
                 barcodes_file, project_dir)


def anemone_pipeline(
        input1_list, input2_list, mismatch,
        barcodes_dict, project_dir, input1):
    '''
    composer entry point to anemone
    '''
    project_dir = project_dir + '/' + os.path.basename(input1)
    os.mkdir(project_dir)
    for k, v in barcodes_dict.items():
        if k == os.path.basename(input1):
            barcodes_file = v
    try:
        input2 = input2_list[input1_list.index(input1)]
        output2 = os.path.basename(input2)
    except:
        input2 = False
        output2 = False
    output1 = os.path.basename(input1)
    chunk = 3000000
    anemone_init(input1, input2, output1, output2, mismatch, chunk,
                 barcodes_file, project_dir)


def anemone_init(
        input1, input2, output1, output2, mismatch, chunk,
        barcodes_file, project_dir):
    '''
    extract barcodes from barcodes_file and detect dual-indexing
    call paired-end or single-end anemone
    '''
    barcodes_matrix, R1_barcodes, R2_barcodes, dual_index = barcode_reader(barcodes_file)
    outfile1_list = [open(project_dir + '/temp_unknown.' + output1, 'w')]
    try:
        outfile2_list = [open(project_dir + '/temp_unknown.' + output2, 'w')]
    except:
        outfile2_list = []
    for i, item in enumerate(R1_barcodes):
        outfile1_list.append(open(project_dir + '/' + str(i) + '.' + output1, 'w'))
        try:
            outfile2_list.append(open(project_dir + '/' + str(i) + '.' + output2, 'w'))
        except:
            pass
    if input2:
        anemone(
            input1, input2, output1, output2, outfile1_list,
            outfile2_list, mismatch, 3000000, R1_barcodes, project_dir, True)
        if dual_index == True:
            dual_indexer(
                barcodes_matrix, R2_barcodes, project_dir,
                outfile1_list, outfile2_list, mismatch)
        elif dual_index == False:
            for i, filename in enumerate(outfile1_list[1:]):
                os.rename(filename.name, project_dir + '/' + barcodes_matrix[0][i])                
    elif input2 == False:
        anemone_single(
            input1, output1, outfile1_list, mismatch,
            chunk, R1_barcodes, project_dir, True)


def barcode_reader(barcodes_file):
    '''
    open user-defined barcode file and extract forward and reverse
    barcodes create matrix to pull corresponding sample ids that match
    barcodes
    '''
    barcodes_matrix = []
    with open(barcodes_file) as f:
        line = f.readline()
        R1_barcodes, R2_barcodes = {}, {}
        for i, item in enumerate(line.split()):
            R2_barcodes[item] = i
        if len(R2_barcodes) == 1:
            dual_index = False
        else:
            dual_index = True
        barcodes_matrix = [[0] * 1 for i in range(len(R2_barcodes))]
        for i, line in enumerate(f):
            for j, item in enumerate(line.split()):
                if j == 0:
                    R1_barcodes[item] = i
                if j > 0:
                    barcodes_matrix[j - 1].append(item)
        for i, line in enumerate(barcodes_matrix):
            del line[0]
    return barcodes_matrix, R1_barcodes, R2_barcodes, dual_index


def dual_indexer(
        barcodes_matrix, R2_barcodes, project_dir,
        outfile1_list, outfile2_list, mismatch):
    outfile1_master = []
    outfile2_master = []
    for file1, file2 in zip(outfile1_list[1:], outfile2_list[1:]):
        input2 = file1.name
        input1 = file2.name # the outfile2 results become input1...
        output1 = os.path.basename(input1)
        output2 = os.path.basename(input2)
        outfile1_final_list = [open(project_dir + '/temp_unknown.' + output1, 'w')]
        outfile2_final_list = [open(project_dir + '/temp_unknown.' + output2, 'w')]
        for i, item in enumerate(R2_barcodes):
            outfile1_final_list.append(open(project_dir + '/' + str(i) + '.' + output1, 'w'))
            outfile2_final_list.append(open(project_dir + '/' + str(i) + '.' + output2, 'w'))
            outfile1_master.append(open(project_dir + '/' + str(i) + '.' + output1, 'w'))
            outfile2_master.append(open(project_dir + '/' + str(i) + '.' + output2, 'w'))
        anemone(
            input1, input2, output1, output2, outfile1_final_list,
            outfile2_final_list, mismatch, 3000000, R2_barcodes, project_dir,
            True)
        os.remove(input1)
        os.remove(input2)
        
    for file1, file2 in zip(outfile2_master, outfile1_master):
        for i, element in enumerate(os.path.basename(file1.name).split('.')):
            if i == 0:
                x = element
            if i == 1:
                y = element
            if i > 1:
                print(element)
        print(x, y)
        print(barcodes_matrix[int(x)][int(y)])

#
# KEEP THIS FUNCTION IN CASE BINNING BY SAMPLE ID BECOMES RELEVANT
#
# def dual_indexer(barcodes_matrix, R1_barcodes, R2_barcodes, project_dir, infile1_list, infile2_list):
    # '''
    # handles each R2 output from first round of demultiplexing as new input
    # creates new subdir and avoids sample id redundancies
    # '''
    # print('project dir is ' + project_dir)
    # for i, filename in enumerate(infile2_list[1:]):
        # #TODO account for input2...
        # input1 = filename.name
        # final_dir = project_dir + '/' + os.path.basename(input1)
        # print(final_dir)
        # pre_output, counter = [], 0
        # for index, key in enumerate(R2_barcodes):
            # if barcodes_matrix[index][i] not in pre_output:
                # pre_output.append(barcodes_matrix[index][i])
                # counter += 1
            # R2_barcodes[key] = counter
        # print(R2_barcodes)
        # print(pre_output)
        # #TODO create binned barcodes list and write to sep folders per input
        # for j in pre_output:
            # pass


def anemone(
        input1, input2, output1, output2, outfile1_list,
        outfile2_list, mismatch, chunk, barcodes, project_dir, round_one):
    '''
    use active 'input1' file to demultiplex in a number of ways
    '''
    row_len = len(barcodes) + 1
    matrix_one = matrix_maker(row_len)
    matrix_two = matrix_maker(row_len)
    i, y, entry1, entry2 = 0, 0, "", ""
    with open(input1) as f1, open(input2) as f2:
        for line1, line2 in zip(f1, f2):
            i += 1
            y += 1
            if y == 2 and round_one == True:
                z, output_prefix = exact_matches(line1, barcodes)
            if y == 2 and round_one == False:
                z, output_prefix = mismatches(line1, barcodes, mismatch)
            if y == 2 or y == 4:
                line1 = line1[z:]
            entry1 = entry1 + line1
            entry2 = entry2 + line2
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                matrix_two[output_prefix].append(entry2)
                y, entry1, entry2 = 0, "", ""
            if i == chunk:
                unload(matrix_one, outfile1_list)
                unload(matrix_two, outfile2_list)
                i = 0
                matrix_one = matrix_maker(row_len)
                matrix_two = matrix_maker(row_len)
        unload(matrix_one, outfile1_list)
        unload(matrix_two, outfile2_list)
    if round_one == True:
        if mismatch > 0:
            second_pass(
                output1, output2, outfile1_list, outfile2_list,
                mismatch, chunk, barcodes, project_dir)
        else:
            for x in outfile1_list:
                x.close()
            for x in outfile2_list:
                x.close()
            os.rename(project_dir + '/temp_unknown.' + output1,
                      project_dir + '/unknown.' + output1)
            os.rename(project_dir + '/temp_unknown.' + output2,
                      project_dir + '/unknown.' + output2)
    else:
        os.remove(project_dir + '/temp_unknown.' + output1)
        os.remove(project_dir + '/temp_unknown.' + output2)
        for x in outfile1_list:
            x.close()
        for x in outfile2_list:
            x.close()


def anemone_single(
        input1, output1, outfile1_list, mismatch, chunk,
        barcodes, project_dir, round_one):
    '''
    use active 'input1' file to demultiplex in a number of ways
    '''
    row_len = len(barcodes) + 1
    matrix_one = matrix_maker(row_len)
    i, y, entry1 = 0, 0, ""
    with open(input1) as f1:
        for line1 in f1:
            i += 1
            y += 1
            if y == 2 and round_one == True:
                z, output_prefix = exact_matches(line1, barcodes)
            if y == 2 and round_one == False:
                z, output_prefix = mismatches(line1, barcodes, mismatch)
            if y == 2 or y == 4:
                line1 = line1[z:]
            entry1 = entry1 + line1
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                y, entry1 = 0, ""
            if i == chunk:
                unload(matrix_one, outfile1_list)
                i = 0
                matrix_one = matrix_maker(row_len)
        unload(matrix_one, outfile1_list)
    if round_one == True:
        if mismatch > 0:
            second_pass(
                output1, False, outfile1_list, False, mismatch,
                chunk, barcodes, project_dir)
        else:
            for x in outfile1_list:
                x.close()
            os.rename(project_dir + '/temp_unknown.' + output1,
                      project_dir + '/unknown.' + output1)
    else:
        os.remove(project_dir + '/temp_unknown.' + output1)
        for x in outfile1_list:
            x.close()


def matrix_maker(row_len):
    '''
    create empty matrix
    '''
    matrix = [[0] * 1 for i in range(row_len)]
    for i in range(row_len):
        matrix[i][0] = ''
    return matrix


def exact_matches(line1, barcodes):
    '''
    write to file only barcodes with 100 percent match in first pass
    '''
    for x, file_prefix in barcodes.items():
        if line1.startswith(x): 
            output_prefix = file_prefix + 1
            z = len(x)
            break
        else:
            z = 0
            output_prefix = 0
    return z, output_prefix


def mismatches(line1, barcodes, mismatch):
    '''
    if mismatch > 0 write to file barcodes with leniency
    '''
    z, multi, output_prefix = 0, 0, 0
    for x, file_prefix in barcodes.items():
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


def unload(matrix, outfiles):
    '''
    write matrix out to file once chunk value is achieved
    '''
    for x, fo in enumerate(outfiles):
        fo.write(str(''.join(matrix[x])))


def second_pass(
        output1, output2, outfile1_list, outfile2_list, mismatch,
        chunk, barcodes, project_dir):
    '''
    after the first round of precision demultiplexing, attempt matches
    of unknown reads
    '''
    outfile1_list[0].close()
    outfile1_list[0] = open(project_dir + '/unknown.' + output1, 'w')
    if output2:
        outfile2_list[0].close()
        outfile2_list[0] = open(project_dir + '/unknown.' + output2, 'w')
    input1 = project_dir + '/temp_unknown.' + output1
    if output2:
        input2 = project_dir + '/temp_unknown.' + output2
        anemone(input1, input2, output1, output2, outfile1_list,
                outfile2_list, mismatch, chunk, barcodes, project_dir, False)
    else:
        anemone_single(input1, output1, outfile1_list, mismatch,
                       chunk, barcodes, project_dir, False)

if __name__ == '__main__':
    anemone_main()
