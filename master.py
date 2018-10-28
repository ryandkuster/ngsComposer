import sys
import os
from project_test import file_type
from trimmer import trimmer

b = 45
e = 45


def fastq_finder():
    fastq_list, gz = [], []
    for filename in os.listdir(project_dir):
        try:
            fastq_test = file_type.is_fastq(os.path.abspath(filename))
            gz.append(0)
        except UnicodeDecodeError:
            fastq_test = file_type.is_gzipped(os.path.abspath(filename))
            gz.append(1)
        if fastq_test is None:
            raise TypeError
        fastq_list.append(filename)
    return fastq_list, gz

if __name__ == '__main__':
    project_dir = sys.argv[1]
    assert os.path.exists(project_dir),'project directory not found'
    fastq_list, gz = fastq_finder()
    
    print(fastq_list)
    print(gz)
    print(b)
    print(e)
    
    for x, filename in enumerate(fastq_list):
        trimmer.test(filename, b, e, gz[x])
