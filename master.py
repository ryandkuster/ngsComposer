import sys
import os
from project_test import file_type


if __name__ == '__main__':
    project_dir = sys.argv[1]
    assert os.path.exists(project_dir),'project directory not found'
    fastq_list = []
    for filename in os.listdir(project_dir):
        ext = filename[-5:]
        if ext == 'fastq':
            fastq_list.append(os.path.abspath(filename))
    print(fastq_list)
    for fastq in fastq_list:
        fastq_test = file_type.is_fastq(fastq)
        if fastq_test is None:
            raise TypeError
        else:
            print('you have happy, healthy fastqs')
    print('nothing was killed')
