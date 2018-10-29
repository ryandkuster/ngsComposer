import sys
import os
from project_test import file_type
from project_test import reader
from trimmer import trimmer
from composer import composer

# user must input full path to project folder
b = 6 # this will be found in input file
e = 0 # this will be found in input file
config_R1_barcodes = True # this will be found in input file
config_R2_barcodes = False # this will be found in input file

if __name__ == '__main__':
    project_dir = sys.argv[1]
    assert os.path.exists(project_dir),'project directory not found'
    config = reader.config_reader(project_dir)
    reader.barcode_reader(project_dir, config_R1_barcodes, 'R1')
    reader.barcode_reader(project_dir, config_R2_barcodes, 'R2')
    fastq_list, gz, pairs = reader.fastq_reader(project_dir)
    for x, input_file in enumerate(fastq_list):
        output_file = os.path.basename(input_file)
        trimmer.test(input_file, b, e, output_file, gz[x])
    print(pairs)
