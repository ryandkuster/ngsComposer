import sys
import os
from project_test import file_type
from trimmer import trimmer

# user must input full path to project folder
b = 6 # this will be found in input file
e = 0 # this will be found in input file

def config_reader(project_dir):
    try:
        with open(project_dir + 'config.txt') as f:
            for line in f:
                print(line)
    except FileNotFoundError:
        print(
              '''
              your project directory must contain a 
              configuration file named 'config.txt'
              '''
             )

def fastq_reader():
    fastq_list, gz, pairs = [], [], {}
    for filename in os.listdir(project_dir):
        if filename == 'config.txt':
            pass
        else:
            try:
                fastq_test, pairs = \
                file_type.is_fq(os.path.abspath(filename), pairs)
                gz.append(0)
            except UnicodeDecodeError:
                fastq_test, pairs = \
                file_type.is_gz(os.path.abspath(filename), pairs)
                gz.append(1)
            if fastq_test is None:
                raise TypeError
            fastq_list.append(os.path.abspath(filename))
    return fastq_list, gz, pairs

if __name__ == '__main__':
    project_dir = sys.argv[1]
    assert os.path.exists(project_dir),'project directory not found'
    config = config_reader(project_dir)
    fastq_list, gz, pairs = fastq_reader()
    for x, input_file in enumerate(fastq_list):
        output_file = os.path.basename(input_file)
        trimmer.test(input_file, b, e, output_file, gz[x])
    print(pairs)
