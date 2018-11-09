import sys
import os
from conf import *
from project_test import file_type
from project_test import reader
from trimmer import trimmer
from composer import composer

print(cutoff)

def pipe_writer_forward1(input1, input2, cutoff, chunk, barcodes, project_dir):
    output1 = input1
    output2 = input2
    input1 = project_dir + '/trimmer_' + input1
    input2 = project_dir + '/trimmer_' + input2
    composer.line_writer(input1, input2, cutoff, chunk, barcodes,
                        project_dir, output1, output2)

if __name__ == '__main__':
    if os.path.exists(project_dir) == True:
        project_dir = os.path.abspath(project_dir)
    else:
        sys.exit('''
                project directory not found
                ''')
    R1_barcodes = reader.barcode_reader(project_dir, R1_barcodes, 'R1')
    R2_barcodes = reader.barcode_reader(project_dir, R2_barcodes, 'R2')
    fastq_list, gz, pairs = reader.fastq_reader(project_dir)

    # 'fastq_list' will be used to find files repeatedly

    if trim == True:
        for x, input_file in enumerate(fastq_list):
            output_file =  project_dir + '/trimmer_' + os.path.basename(input_file)
            trimmer.test(input_file, b, e, output_file, gz[x])


    input1_list, input2_list = [], []
    for values in pairs.values():
        if len(values) == 2:
            for filename in values:
                with open(filename) as f:
                    header = f.readline()
                    for pos, x in enumerate(header):
                        if x == ' ':
                            space_pos = pos
                    end = header[space_pos+1]
                    if int(end) == 1:
                        input1_list.append(os.path.basename(filename))
                    if int(end) == 2:
                        input2_list.append(os.path.basename(filename))
        else:
            print('is this a single end library?')


    print(input1_list)
    print(input2_list)
    if R1_barcodes == True:
        for x, input1 in enumerate(input1_list):
            input1 = input1
            input2 = input2_list[x]
            pipe_writer_forward1(input1, input2, cutoff, 200000, R1_barcodes, project_dir)
