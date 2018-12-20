# provide the full path to your fastq files
project_dir = '/home/ryan/Testing'
# project_dir = '/media/sf_E_DRIVE/Analysis/Testing'
# project_dir = '/home/rkuster/Desktop/pipeline_test'

# if fastq files should be treated as paired ends, use 'True', else 'False'
paired = True

# choose number of subprocesses
threads = 1

# positions to trim from front of read before demultiplexing, leave 0 if no buffer sequence
front_trim = 6

# name of  file with barcodes to demultiplex forward reads
barcodes_index = 'index.txt'

# number of mismatches (hamming distance) allowed in barcodes
mismatch = 1

# list sequences immediately adjacent to barcodes
overhang_list = ['TCC','TCT']

# optionally remove each intermediate file
remove_intermediates = True

# optionally remove files that do not successfully match desired characteristics
remove_fail = True
