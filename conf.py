# provide the full path to your fastq files (only place full path is used)
proj_dir = '/home/descartes/GitHub/composer/examples/project2'

# if fastq files should be treated as paired ends, use 'True', else 'False'
paired = True

# choose number of subprocesses
threads = 2

# create initial QC output for reads and pause for user
initial_qc = False

# positions to trim from front of read before demultiplexing, leave 0 if no buffer sequence
front_trim = 6

# name of  file with barcodes to demultiplex forward reads (use 'False' if not demultiplexing)
bcs_index = 'index.txt'

# number of mismatches (hamming distance) allowed in barcodes
mismatch = 1

# list sequences immediately adjacent to barcodes
overhang_ls = ['TCC','TCT']

# additional, non-genomic base not found in barcode sequence
non_genomic = 'T'

# quality score minimum (Phred value 0-40)
q_min = 30

# percentage of reads containing q_min or high quality scores
q_percent = 95

# pause at QC steps for manual input (overwrites q_min and q_percent)
manual_entry = True

# optionally remove each intermediate file
remove_intermediates = True

# optionally remove files that do not successfully match desired characteristics
remove_fail = True
