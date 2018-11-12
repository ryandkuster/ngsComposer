import sys
import os

def comp_main():
    input1 = sys.argv[1] # R1 reads
    input2 = sys.argv[2] # R2 reads
    mismatch = int(sys.argv[3]) # number of mismatches allowed
    barcode_file = sys.argv[4] # barcodes file
    chunk = 3000000 # how many reads to process before writing
    with open(barcode_file) as f:
        barcodes = []
        for line in f:
            barcodes.append(line.rstrip())
    project_dir = os.path.dirname(os.path.abspath(input1))
    composer(input1, input2, mismatch, chunk, barcodes, project_dir)

def comp_piper(input1_list, input2_list, mismatch, R1_barcodes, project_dir, input1):
    input2 = input2_list[input1_list.index(input1)] 
    composer(input1, input2, mismatch, 3000000, R1_barcodes, project_dir)

def composer(input1, input2, mismatch, chunk, barcodes, project_dir):
    # consider with open(all outputs, 'w') here,then indent the next with open
    output1 = os.path.basename(input1)
    output2 = os.path.basename(input2)
    row_len = len(barcodes) + 1
    matrix_one = matrix_maker(row_len)
    matrix_two = matrix_maker(row_len)
    i, y, entry1, entry2 = 0, 0, "", ""
    outfile1_list, outfile2_list  = [], []
    for row in range(row_len):
        outfile1_list.append(open(project_dir + '/' + str(row + 1) + '_' + output1, 'w'))
        outfile2_list.append(open(project_dir + '/' + str(row + 1) + '_' + output2, 'w'))
    with open(input1) as f1, open(input2) as f2:
        for line1 in f1:
            i += 1
            y += 1
            if y == 2:
                for file_prefix, x in enumerate(barcodes):
                    hamm, z = 0, 0
                    if line1.startswith(x): 
                        output_prefix = file_prefix
                        z = len(x)
                        break
                    else:
                        output_prefix = -1
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
                unload(matrix_one, matrix_two, row_len, outfile1_list, outfile2_list)
                i = 0

                matrix_one = matrix_maker(row_len)
                matrix_two = matrix_maker(row_len)
        unload(matrix_one, matrix_two, row_len, outfile1_list, outfile2_list)
    # close the unmatched file here?
    for x in outfile1_list:
        x.close()
    for x in outfile2_list:
        x.close()

def matrix_maker(row_len):
    matrix = [[0] * 1 for i in range(row_len)]
    for i in range(row_len):
        matrix[i][0] = ''
    return matrix

def unload(matrix_one, matrix_two, row_len, outfile1_list, outfile2_list):
    for x, fo in enumerate(outfile1_list):
        fo.write(str(''.join(matrix_one[x])))
    for x, fo in enumerate(outfile2_list):
        fo.write(str(''.join(matrix_two[x])))


if __name__ == '__main__':
    comp_main()
