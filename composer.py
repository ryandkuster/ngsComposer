import sys

class composer:
    def matrix_maker(row_len):
        matrix = [[0] * 1 for i in range(row_len)]
        counter = 0
        for i in range(row_len):
            counter += 1
            matrix[i][0] = counter
        return matrix

    def unload(matrix_one, matrix_two, row_len, output1, output2, project_dir):
        for row in range(row_len):
            output_data1 = matrix_one[row]
            output_data2 = matrix_two[row]
            with open(project_dir + '/' + str(output_data1[0]) + '_' + output1, 'a') as outfile1:
                outfile1.write(str(''.join(output_data1[1:])))
            with open(project_dir + '/' + str(output_data2[0]) + '_' + output2, 'a') as outfile2:
                outfile2.write(str(''.join(output_data2[1:])))

    def line_writer(input1, input2, cutoff, chunk, barcodes, project_dir):
        output1 = input1
        output2 = input2
        input1 = project_dir + '/trimmer_' + input1
        input2 = project_dir + '/trimmer_' + input2
        with open(input1) as f1, open(input2) as f2:
            row_len = len(barcodes) + 1
            matrix_one = composer.matrix_maker(row_len)
            matrix_two = composer.matrix_maker(row_len)
            i, y, z, entry1, entry2 = 0, 0, 0, "", ""
            for line1 in f1:
                i += 1
                y += 1
                if y == 2:
                    for file_prefix, x in enumerate(barcodes):
                        hamm = 0
                        for j in range(len(x)):
                            if x[j] != line1[j]:
                                hamm = hamm + 1
                                break
                        if hamm == 0:
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
                    y, z, entry1, entry2 = 0, 0, "", ""
                if i == chunk:
                    composer.unload(matrix_one, matrix_two, row_len, output1, output2, project_dir)
                    i = 0
                    matrix_one = composer.matrix_maker(row_len)
                    matrix_two = composer.matrix_maker(row_len)
            composer.unload(matrix_one, matrix_two, row_len, output1, output2, project_dir)

if __name__ == '__main__':
    input1 = sys.argv[1] # R1 reads
    input2 = sys.argv[2] # R2 reads
    cutoff = int(sys.argv[3]) # number of mismatches allowed
    chunk = 4*(int(sys.argv[4])) # how many reads to process before writing, 20000 is good so far
    barcode_file = sys.argv[5] # barcodes file
    with open(barcode_file) as f:
        barcodes = []
        for line in f:
            barcodes.append(line.rstrip())
    composer.line_writer(input1, input2, cutoff, chunk, barcodes)
