# provide the full path to your fastq files
project_dir = '/home/ryan/Testing'
# project_dir = '/media/sf_E_DRIVE/Analysis/Testing'

# if fastq files are only paired and should be treated as such, use 'True', else 'False'
paired = True

# choose number of subprocesses
threads = 1

# positions to trim from front of read before demultiplexing, leave 0 if no buffer sequence
front_trim = 0

# positions to trim from end of read before demultiplexing, leave 0 if no buffer sequence
back_trim = 0

# name of  file with barcodes to demultiplex forward reads
barcodes_file = 'key.txt'

# if forward and reverse reads each have corresponding sets of barcodes
dual_index = False

# number of mismatches (hamming distance) allowed in barcodes
mismatch = 0

# list sequences immediately adjacent to barcodes
overhang_list = ['TCC','TCT']

# optionally remove each intermediate file
remove_intermediates = True

# optionally remove files that do not successfully match desired characteristics
remove_fail = True
