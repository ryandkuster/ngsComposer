import sys
import os


def anemone_main():
    mismatch = int(sys.argv[1]) # number of mismatches allowed
    barcodes_file = sys.argv[2] # barcodes file
    input1 = sys.argv[3] # R1 reads
    output1 = os.path.basename(input1)
    try:
        input2 = sys.argv[4] # R2 reads
        output2 = os.path.basename(input2)
        paired = True
    except:
        paired = False
        pass
    chunk = 3000000 # how many reads to process before writing
    project_dir = os.getcwd()
    barcodes_matrix, R1_barcodes, R2_barcodes = barcode_reader(barcodes_file)
    anemone_init(input1, input2, output1, output2, mismatch, chunk, barcodes_matrix, project_dir, paired)


def anemone_pipeline(input1_list, input2_list, paired, mismatch, barcodes_index, project_dir, input1):
    if paired == True:
        input2 = input2_list[input1_list.index(input1)]
        output2 = os.path.basename(input2)
    else:
        input2 = False
        output2 = False
    output1 = os.path.basename(input1)
    chunk = 3000000
    #TODO get barcodes_file based on barcodes_index matching input1
    barcodes_matrix, R1_barcodes, R2_barcodes = barcode_reader(barcodes_file)
    anemone_init(input1, input2, output1, output2, mismatch, chunk, barcodes_matrix, project_dir, paired)


def barcode_reader(barcodes_file):
    '''
    open user-defined barcode file and extract forward and reverse barcodes
    create matrix to pull corresponding sample ids that match barcodes
    '''
    barcodes_matrix = []
    with open(barcodes_file) as f:
        R2_barcodes = R2_barcodes_maker([], f)
        barcodes_matrix = array_maker(R2_barcodes)
        barcodes_matrix, R1_barcodes = barcodes_matrix_maker([], f, barcodes_matrix)
    return barcodes_matrix, R1_barcodes, R2_barcodes


def R2_barcodes_maker(R2_barcodes, f):
    '''
    populate list of R2 barcodes
    '''
    line = f.readline()
    for item in line.split():
        R2_barcodes.append(item)
    if len(R2_barcodes) == 1 and dual_index == True:
        sys.exit("expected barcodes file to have reverse barcodes")
    return R2_barcodes


def array_maker(R2_barcodes):
    '''
    create empty matrix with dimensions 1 by length R2 barcodes 
    '''
    barcodes_matrix = [[0] * 1 for i in range(len(R2_barcodes))]
    for i, x in enumerate(R2_barcodes):
        barcodes_matrix[i][0] = x
    return barcodes_matrix


def barcodes_matrix_maker(R1_barcodes, f, barcodes_matrix):
    '''
    create a matrix of user-defined sample ids and populate list of R1 barcodes
    '''
    for i, line in enumerate(f):
        for j, item in enumerate(line.split()):
            if j == 0:
                R1_barcodes.append(item)
            if j > 0:
                barcodes_matrix[j - 1].append(item)    
    return barcodes_matrix, R1_barcodes


def anemone_init_paired(input1, input2, output1, output2, mismatch, chunk, barcodes_matrix, project_dir):
    #TODO update this to reflect sample ids in matrix
    filename = os.path.basename(input1)
    for i, item in enumerate(barcodes_matrix):
        if item[0] == filename:
            sample_id = barcodes_matrix[i][1:]
            R1_barcodes = barcodes_matrix[0][1:]
            na_list = []
    for i, item in enumerate(sample_id):
        if item == 'na':
            na_list.append(i)
    for i in sorted(na_list, reverse=True):
        del sample_id[i]
        del R1_barcodes[i]
    row_len = len(R1_barcodes) + 1
    outfile1_list = [open(project_dir + '/temp_unknown_' + output1, 'w')]
    outfile2_list = [open(project_dir + '/temp_unknown_' + output2, 'w')]
    for id in sample_id:
        outfile1_list.append(open(project_dir + '/' + str(id) + '_' + output1, 'w'))
        outfile2_list.append(open(project_dir + '/' + str(id) + '_' + output2, 'w'))
    #TODO call anemone from here, followed by mismatcher, then second pass
    anemone(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, 3000000, R1_barcodes, project_dir, True)


def anemone(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir, round_one):
    #TODO get to work in standard fashion, then break into functions, test speed
    matrix_one = matrix_maker(row_len)
    matrix_two = matrix_maker(row_len)
    i, y, entry1, entry2 = 0, 0, "", ""
    with open(input1) as f1, open(input2) as f2:
        for line1 in f1:
            i += 1
            y += 1
            if y == 2 and round_one == True:
                for file_prefix, x in enumerate(barcodes):
                    if line1.startswith(x): 
                        output_prefix = file_prefix + 1
                        z = len(x)
                        break
                    else:
                        z = 0
                        output_prefix = 0
            if y == 2 and round_one == False:
                z, multi, output_prefix = 0, 0, 0
                for file_prefix, x in enumerate(barcodes):
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
            if y == 2 or y == 4:
                line1 = line1[z:]
            entry1 = entry1 + line1
            for line2 in f2:
                entry2 = entry2 + line2
                break
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                matrix_two[output_prefix].append(entry2)
                y, entry1, entry2 = 0, "", ""
            if i == chunk:
                unload(matrix_one, row_len, outfile1_list)
                unload(matrix_two, row_len, outfile2_list)
                i = 0
                matrix_one = matrix_maker(row_len)
                matrix_two = matrix_maker(row_len)
        unload(matrix_one, row_len, outfile1_list)
        unload(matrix_two, row_len, outfile2_list)
    if round_one == True:
        if mismatch > 0:
            mismatcher(row_len, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir)
        else:
            for x in outfile1_list:
                x.close()
            for x in outfile2_list:
                x.close()
            os.rename(project_dir + '/temp_unknown_' + output1, project_dir + '/unknown_' + output1)
            os.rename(project_dir + '/temp_unknown_' + output2, project_dir + '/unknown_' + output2)
    else:
        os.remove(project_dir + '/temp_unknown_' + output1)
        os.remove(project_dir + '/temp_unknown_' + output2)
        for x in outfile1_list:
            x.close()
        for x in outfile2_list:
            x.close()


def matrix_maker(row_len):
    '''
    create empty matrix
    '''
    matrix = [[0] * 1 for i in range(row_len)]
    for i in range(row_len):
        matrix[i][0] = ''
    return matrix


def unload(matrix, row_len, outfiles):
    '''
    write matrix out to file
    '''
    for x, fo in enumerate(outfiles):
        fo.write(str(''.join(matrix[x])))


def mismatcher(row_len, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir):
    '''
    after the first round of precision demultiplexing, attempt matches of unknown reads
    '''
    outfile1_list[0].close()
    if paired == True:
        outfile2_list[0].close()
    outfile1_list[0] = open(project_dir + '/unknown_' + output1, 'w')
    if paired == True:
        outfile2_list[0] = open(project_dir + '/unknown_' + output2, 'w')
    input1 = project_dir + '/temp_unknown_' + output1
    if paired == True:
        input2 = project_dir + '/temp_unknown_' + output2
    anemone(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir, False)


if __name__ == '__main__':
    anemone_main()
