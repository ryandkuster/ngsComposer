# if fastq files should be treated as paired ends, use 'True', else 'False'
paired = True

# choose number of subprocesses that should run simultaneously
procs = 2

# create initial QC output
initial_qc = True

# perform qc step at each filtering stage (time-consuming, but informative)
walkthrough = True

# run from beginning to end without pausing at qc steps
walkaway = True

# positions to trim from front of read before demultiplexing, leave 0 if no buffer sequence
front_trim = 2

# name of  file with barcodes to demultiplex forward reads (use 'False' if not demultiplexing)
bcs_index = 'index.txt'

# number of mismatches (hamming distance) allowed in barcodes
mismatch = 1

# list sequences immediately adjacent to barcodes
bases_ls = ['A', 'C', 'G', 'T']

# additional, non-genomic base not found in barcode sequence (e.g. 'T' complementary to A-tailing library prep)
non_genomic = 'T'

# quality score minimum (Phred value 0-40)(use 'False' to skip)
q_min = 30

# percentage of reads >/= q_min quality scores (use 'False' to skip)
q_percent = 95

# optionally remove each transitional file folder to save space
rm_transit = True
