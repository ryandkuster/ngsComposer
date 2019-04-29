<a href="https://imgur.com/gYTN9Hv"><img src="https://i.imgur.com/gYTN9Hv.png" title="source: imgur.com" /></a>

# composer: a deep QC adventure!

Base-call error-filtering and read preprocessing pipeline designed by biologists

## Features

- Few dependencies (Python3 and R)
- Easy to learn
- Supports variable length barcodes and dual-indexing
- Trims buffer sequences and quality filters on a read-by-read basis
- Accepts project directory of multiple libraries
- Designed by biologists (please don't run away!)

## Installation

Clone or download the Git repository to your desired tool folder

```bash
$ git clone https://github.com/ryandkuster/composer.git
```

## Usage

### Basic usage

Set up your project directory containing the following files:
- 1_R1.fastq
- conf.py (see "Configuration" below for detailed instructions)

Optionally, paired end files can also be included in the project directory:
- 1_R2.fastq

From command line, run composer with the specified directory of your project
```bash
$ python3 composer.py <path_to_directory>
```

***

### Configuration
Using a text editor, save a file containing the following variables as a python file (inludes '.py' as file extension) and include it in your project directory:

```
# if fastq files should be treated as paired ends (True or False)
paired = True

# choose number of subprocesses that should run simultaneously (1 or greater)
procs = 2

# create initial QC output (True or False)
initial_qc = False

# perform qc step at each filtering stage (time-consuming, but informative) (True or False)
walkthrough = False

# run from beginning to end without pausing at qc steps (True or False)
walkaway = True

# positions to trim from front of read before demultiplexing, (0 if no buffer sequence)
front_trim = 2

# name of file associating forward reads with barcode array (False if not demultiplexing)
bcs_index = 'index.txt'

# number of mismatches allowed in barcodes (hamming distance)
mismatch = 1

# list sequences immediately adjacent to barcodes (False if not used)
bases_ls = ['TCGAG','AGCTT']

# additional, non-genomic base not found in barcode sequence (e.g. 'T' complementary to A-tailing library prep)
non_genomic = 'T'

# quality score minimum (Phred value 0-40)(use False to skip)
q_min = 30

# percentage of reads >/= q_min quality scores (use False to skip)
q_percent = 95

# trim read 3' ends based on q_min threshold (True or False)
end_trim = True

# optionally remove each transitional file folder to save space (True or False)
rm_transit = True
```

***

### Barcode demultiplexing
#### Barcodes file
Optionally, one or more barcode files may be included in the project directory for demultiplexing. The following files are required at minimum:
- barcodes_1.txt
- index.txt

The barcodes file is a tab or space delimited file including forward barcodes as rows and reverse barcodes as columns. For example, the following would be required for a dual-indexed library:

```
	A	C	G	T
A	sample1	sample5	sample6	sample10
C	sample2	sample5	sample7	sample10
G	sample3	sample5	sample8	sample10
T	sample4	sample5	sample9	sample10
```

If reverse barcodes do not require demultiplexing, the barcode file can be set up as follows with "NA" or any other text used as a header in the first row:

```
	NA
A	sample1
C	sample2
G	sample3
T	sample4
```

***

#### Index file
The index file is a tab delimited file required to associate the barcodes file with a specific library in your project directory. It must include the filename followed by the appropriate barcodes file:

```
1_R1.fastq  barcodes_1.txt
```

Alternatively, multiple barcoding schemes may be included to accomodate multiple libraries. For example:

```
1_R1.fastq  barcodes_1.txt
2_R1.fastq  barcodes_2.txt
3_R1.fastq  barcodes_3.txt
```

The name of your index file must match the exact text for 'bcs_index' in the conf.py file.

***

## License

<a href="https://github.com/ryandkuster/Pipeline/blob/master/LICENSE">Apache License Version 2.0</a>

<a href="https://imgur.com/uQ0kCRk"><img src="https://i.imgur.com/uQ0kCRk.png" title="source: imgur.com" /></a>
