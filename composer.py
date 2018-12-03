import sys
import os


def comp_main():
    mismatch = int(sys.argv[1]) # number of mismatches allowed
    barcode_file = sys.argv[2] # barcodes file
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
    with open(barcode_file) as f:
        barcodes = []
        for line in f:
            barcodes.append(line.rstrip())
    project_dir = os.getcwd()
    if paired == True:
        comp_init_paired(input1, input2, output1, output2, mismatch, chunk, barcodes, project_dir)
    if paired == False:
        comp_init_single(input1, output1, mismatch, chunk, barcodes, project_dir)


def comp_piper_paired(input1_list, input2_list, mismatch, barcodes_matrix, project_dir, input1):
    input2 = input2_list[input1_list.index(input1)]
    output1 = os.path.basename(input1)
    output2 = os.path.basename(input2)
    chunk = 3000000
    comp_init_paired(input1, input2, output1, output2, mismatch, chunk, barcodes_matrix, project_dir)


def comp_piper_single(mismatch, barcodes, project_dir, input1):
    output1 = os.path.basename(input1)
    chunk = 3000000
    comp_init_single(input1, output1, mismatch, chunk, barcodes, project_dir)


def comp_init_paired(input1, input2, output1, output2, mismatch, chunk, barcodes_matrix, project_dir):
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
    composer(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, 3000000, R1_barcodes, project_dir, True)


def comp_init_single(input1, output1, mismatch, chunk, barcodes, project_dir):
    row_len = len(barcodes) + 1
    outfile1_list = [open(project_dir + '/temp_unknown_' + output1, 'w')]
    for row in range(row_len - 1):
        outfile1_list.append(open(project_dir + '/' + str(row + 1) + '_' + output1, 'w'))
    composer_single(row_len, input1, output1, outfile1_list, mismatch, 3000000, barcodes, project_dir, True)


def composer(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir, round_one):
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


def composer_single(row_len, input1, output1, outfile1_list, mismatch, chunk, barcodes, project_dir, round_one):
    matrix_one = matrix_maker(row_len)
    i, y, entry1 = 0, 0, ""
    with open(input1) as f1:
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
            if y == 4:
                matrix_one[output_prefix].append(entry1)
                y, entry1 = 0, ""
            if i == chunk:
                unload(matrix_one, row_len, outfile1_list)
                i = 0
                matrix_one = matrix_maker(row_len)
        unload(matrix_one, row_len, outfile1_list)
    if round_one == True:
        if mismatch > 0:
            mismatcher_single(row_len, output1, outfile1_list, mismatch, chunk, barcodes, project_dir)
        else:
            for x in outfile1_list:
                x.close()
            os.rename(project_dir + '/temp_unknown_' + output1, project_dir + '/unknown_' + output1)
    else:
        os.remove(project_dir + '/temp_unknown_' + output1)
        for x in outfile1_list:
            x.close()


def matrix_maker(row_len):
    matrix = [[0] * 1 for i in range(row_len)]
    for i in range(row_len):
        matrix[i][0] = ''
    return matrix


def unload(matrix, row_len, outfiles):
    for x, fo in enumerate(outfiles):
        fo.write(str(''.join(matrix[x])))


def mismatcher(row_len, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir):
    outfile1_list[0].close()
    outfile2_list[0].close()
    outfile1_list[0] = open(project_dir + '/unknown_' + output1, 'w')
    outfile2_list[0] = open(project_dir + '/unknown_' + output2, 'w')
    input1 = project_dir + '/temp_unknown_' + output1
    input2 = project_dir + '/temp_unknown_' + output2
    composer(row_len, input1, input2, output1, output2, outfile1_list, outfile2_list, mismatch, chunk, barcodes, project_dir, False)


def mismatcher_single(row_len, output1, outfile1_list, mismatch, chunk, barcodes, project_dir):
    outfile1_list[0].close()
    outfile1_list[0] = open(project_dir + '/unknown_' + output1, 'w')
    input1 = project_dir + '/temp_unknown_' + output1
    composer_single(row_len, input1, output1, outfile1_list, mismatch, chunk, barcodes, project_dir, False)


if __name__ == '__main__':
    comp_main()
