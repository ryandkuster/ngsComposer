# INITIAL DEVELOPMENT RELEASE
Please use ngsComposer at your own discretion. Currently, the software is in an early testing release. If you find any issues or have suggestions to improve ngsComposer or its documentation, please contact rkuster@utk.edu or 
bolukolu@utk.edu.

<a><img src="https://i.imgur.com/hqfJVWJ.png" title="source: imgur.com" /></a>

# ngsComposer: empirically different

Base-call error-filtering and read preprocessing pipeline designed by biologists

## Features

- Full start-to-finish pipeline for many library types
- Few dependencies (Python3 and R)
- Easy to learn
- Supports variable length barcodes and dual-indexing
- Trims buffer sequences and quality filters on a read-by-read basis
- Accepts project directory of multiple libraries
- Designed by biologists (please don't run away!)

## Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Demultiplexing](#demultiplexing)
  - [Standalone Tools](#standalone)
- [Troubleshooting](#troubleshooting)
- [Versioning](#versioning)
- [License](#license)


## Installation

Currently, ngsComposer is only available for unix-based systems (i.e. mac and linux).

Clone or download the Git repository to your desired folder

```bash
$ git clone https://github.com/ryandkuster/ngsComposer.git
```

Dependencies:
- Python3 version 3.5 or above
- R version compatible with the following required packackes:
	ggplot2

For help troubleshooting installation, see the troubleshooting section

## Usage

### Basic usage

Set up your project directory containing the following files:
- fastq file(s)
- conf.py (see "Configuration" below for detailed instructions on creating this file)

<a><img src="https://imgur.com/SYiOxfR.png" title="source: imgur.com" width=400 /></a>

Optionally, paired-end files and even multiple separate libraries can be included in the project directory.

From command line, run ngsComposer with the specified directory of your project
```bash
$ python3 <path_to_composer_directory>/composer.py -i <path_to_project_directory>
```

If this is the first time running the pipeline, you may need to wait for R to install the appropriate packages and dependencies.

Several **example datasets** are included in the "examples" directory. Users are encouraged to examine and run these small projects to assist in understanding pipeline functionality.

***

### Overview
The order of steps in the ngsComposer pipeline are outlined in the following figure:

<a><img src="https://i.imgur.com/99BsbsJ.png" title="source: imgur.com" /></a>

The steps implemented are first specified in a configuration file.

***

### Configuration
Using a text editor, save a file containing any of the following variables as a python script called 'conf.py' (includes '.py' as file extension) and include it in your project directory.

|Variable        |Usage           |Input |
|:-------------|:-------------|:-------------|
|paired|if fastq files should be treated as paired ends|True or False|
|procs|choose maximum number of subprocesses that can run simultaneously|integer|
|alt_dir|optional, full path to empty directory for file storage|e.g. '/path/to/dir/' (quotes required)|
|walkaway|run from beginning to end without pausing at qc steps (use True to run without interruption, use False to initiate user-controlled walkthrough mode: requires all_qc)|True or False|
|rm_transit|optional, remove each transitional file folder to save space|True or False|
|initial_qc|create initial QC output|True or False|
|all_qc|perform qc step at each filtering stage (use 'full' to produce visualizations for every file, use 'summary' for a summarized version of the R1, R2, and single reads)|'full' or 'summary' (quotes required)|
|front_trim|positions to trim from front of read before demultiplexing, leave 0 if no buffer sequence|integer|
|mismatch|number of mismatches (hamming distance) allowed in barcodes (must include 'index.txt' and barcodes file(s) in project directory; see "Demultiplexing" below)|integer|
|R1_bases_ls|list expected sequence motifs immediately adjacent to barcodes (e.g. restriction sites)|e.g. ['TCC', 'TCT'] (quotes, commas, and brackets required)|
|R2_bases_ls|list expected sequence motifs immediately adjacent to barcodes (e.g. restriction sites)|e.g. ['TCC', 'TCT'] (quotes, commas, and brackets required)|
|non_genomic|number of non-genomic bases not found in barcode sequence (e.g. 'T' complementary to A-tailing library prep)|integer|
|end_score|end-trim once entire window >= this Q score|integer between 0 and 40|
|window|size of window to test for >= end_trim|integer within read length|
|min_len|minimum read length to retain for end-trimming and adapter removal|integer > 0|
|q_min|Q score minimum (Phred value 0-40) applied to q_percent variable|integer between 0 and 40|
|q_percent|percentage of reads >= q_min Q scores|number between 0 and 100|
|adapter_match|number of base matches to identify adapters (requires 'adapters.txt')|integer (recommend 12)|
|p64|defaults to phred+33, use True if using phred+64 qscores|True or False|
|r_dir|optional, full path to R package installation directory (defaults to the usual local R installation library paths|e.g. '/path/to/dir/' (quotes required)|
|compress|gzip compress files after each step in the pipeline to save space (defaults to True)|True or False|


An example configuration file may look like this:

**conf.py**

```
paired = True
procs = 2
alt_dir = '/home/user/project'
walkaway = False
rm_transit = True
initial_qc = False
all_qc = 'summary'
front_trim = 6
mismatch = 1
R1_bases_ls = ['TCC', 'TCT']
R2_bases_ls = ['TCC', 'TCT']
non_genomic = 1
end_score = 30
window = 10
min_len = 50
q_min = 30
q_percent = 95
```
*In the above example, a paired library will be expected (**paired = True**) and the maximum number of subprocesses spawned will be 2 (**procs = 2**).  The empty directory '/home/user/project' will contain all resulting filtered data (**alt_dir = '/home/user/project'**). The pipeline will pause after relevant steps (**walkaway = False**) so users can view qc plots and have the option of modifying or bypassing the step. To save disk space, transitional directories will be removed (**rm_transit = True**) and only the final filtered data and any qc stats created in the pipeline will remain. No qc statistics will be created on the raw data in this example (**initial_qc = False**), but for all subsequent steps a summarized version will be created that collapses all R1, R2, and/or single-end reads to provide a helpful overview of the results of a given filtering step (**all_qc = 'summary'**).*

*A buffer sequence of length 6 (**front_trim = 6**) will be trimmed before demultiplexing, which will allow mismatch at a hamming distance of 1 (**mismatch = 1**).*

*In this case, samples were double-digested with AluI and HaeIII and A-tailed before adapter ligation (**R1_bases_ls = ['TCC', 'TCT']** and **R2_bases_ls = ['TCC', 'TCT']**). Only reads containing these motifs will pass to subsequent steps. As the T complement from A-tailing introduces an artificial residue not present in the specimen sequenced, it can simultaneously be removed alongside motif detection (**non_genomic = 1**).*

*Automatic end-trimming will be performed based on Q score. Here, groups of bases are considered within a moving window of 10 bases at a time (**window = 10**) until that window consists only of the desired Q score at or above 30 (**end_score = 30**). It is at this point that the read is trimmed. Reads that are less than 50 bp will be discarded (**min_len = 50**)*

*Only reads that have a Q score of 30 (**q_min = 30**) acrosss at least 95 percent of the read (**q_percent = 95**) will pass to subsequent steps. If a R1 read or an R2 read passes while its partner fails, it will be placed into a single-end read subfolder and the failing read will be discarded.*


Alternatively, a configuration file may only need to include necessary components for a run:

**conf.py**

```
paired = True
procs = 8
mismatch = 1
q_min = 30
q_percent = 95
```
*This example will demultiplex paired-end data using a mismatch value of one followed by a threshold filter for reads comprised of base-calls at or above 30 across at least 95 percent of the read. A maximum of 8 subprocesses will be called.*

***

### Demultiplexing
#### Barcodes file(s)
Optionally, one or more barcode files may be included in the project directory for demultiplexing. The following files are required at minimum:
- barcodes_1.txt
- index.txt

<a><img src="https://imgur.com/VlLCeY4.png" title="source: imgur.com" width=400 /></a>

Naming conventions: "index.txt" is required, the barcodes file can be named as desired (see "Index file for directing multiple barcode files")

The barcodes file is a tab or space delimited file with no spaces in sample names (or, copy directly from your favorite spreadsheet program into a text file). Forward barcodes begin each row and reverse barcodes begin each column with the desired sample names indicated in the interior of the matrix. For example, the following would be required for a dual-indexed library:

**barcodes_1.txt**
```
	A	C	G	T
A	sample1	sample5	sample6	sample10
C	sample2	sample5	sample7	sample10
G	sample3	sample5	sample8	sample10
T	sample4	sample5	sample9	sample10
```
*Note that in the example above the reverse barcode "C" corresponds with multiple identical sample names (sample5). While not common practice, ngsComposer accomodates repeated sample names and concatenates accordingly.*

If reverse barcodes do not require demultiplexing, the barcode file can be set up as follows with "NA" or any other text used as a header in the first row:

**barcodes_1.txt**
```
	NA
A	sample1
C	sample2
G	sample3
T	sample4
```

#### Index file for directing multiple barcode files
The index file is a tab delimited file required to associate the barcodes file with a specific library in your project directory. It must include the filename of the forward read (R1) followed by the appropriate barcodes file. Reverse reads (R2), if present, will automatically be detected and are not indicated in this file.

**index.txt**
```
1_R1.fastq  barcodes_1.txt
```

Alternatively, multiple barcoding schemes may be included to accomodate multiple libraries. For example:

**index.txt**
```
1_R1.fastq  1_bcs.txt
2_R1.fastq  2_bcs.txt
```

<a><img src="https://imgur.com/mnrrviL.png" title="source: imgur.com" width=600 /></a>

*In this example, sample "1_R1.fastq" and "1_R2.fastq" correspond with "1_bcs.txt" and "2_R1.fastq" and "2_R2.fastq" correspond with "2_bcs.txt"*

***

### Adapters
#### Adapters file(s)
Optionally, 'adapters.R2.txt' and 'adapters.R1.txt' may be included in the project directory for recognition and removal of adapters. The 'adapters.R2.txt' file contains the adapters expected to appear in the R1 reads. Adapter sequences should be newline-separated and be in 5' to 3' orientation. If libraries are barcoded, users are encouraged to provide adapter sequences that contain the corresponding barcodes expected in the opposing end of the read's adapter.

<a><img src="https://i.imgur.com/3In8TX0.png" title="source: imgur.com" width=300 /></a>


**adapters.R2.txt**
```
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGCTCAGTTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTATCTGACCT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTATATGAGACG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTTATGGAAT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTAATCTCGTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCGCGATGTT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAGAGCACTAG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCCTTGATC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTACTCAGTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCGTCTGACT
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGAACATACGG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCTATGACTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTAATGGCAAG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTGCCGCTTC
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCGGCAATGGA
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGCCGTAACCG
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAACCATTCTC 
```
*Each of the above sample adapters is presented in 5' to 3' orientation and shares a common 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' adapter sequence followed by expected barcodes. Adapter sequences may also include restriction motifs for greater detection, but these sequences will also be removed. Porifera.py creates all reverse-complements before detection.*

<a><img src="https://i.imgur.com/mYBQWIv.png" title="source: imgur.com" width=600 /></a>

*When paired end data is used, as above, 'adapters.R1.txt' and 'adapters.R2.txt' must be provided. Adapters are tested for the inclusion of barcodes and only those combinations of R1/R2 barcodes leading to a given sample will be used to search for adapters quickly and with a lower false positive rate.*

### Standalone
All tools available in the ngsComposer pipeline can be called individually from the command line. Please see the <a href="https://github.com/ryandkuster/composer/tree/master/tools">ngsComposer Standalone Tools page</a> for usage.

## Troubleshooting

### Python installation
To view Python version, from the terminal type:

```bash
$ python3 --version
```

If python3 is not found, you can try one of the python3 releases from the Python Software Foundation <a href="https://www.python.org/downloads/">downloads page</a>.

Alternatively, a package manager is an easy way to install Python from the terminal. For Ubuntu, Python can be installed directly using apt (replace 'X' with an existing version in the apt repository):

```bash
$ sudo apt-get update
$ sudo apt-get install python3.X
```

### R installation
To view R version, from the terminal type:

```bash
$ R --version
```

To install the newest version of R, see the releases available at the Comprehensive R Archive Network <a href="https://cran.r-project.org/">downloads page</a>.

For Ubuntu, R can be installed directly using apt:

```bash
$ sudo apt update
$ sudo apt install r-base
```

...or with homebrew on Mac using:

```
brew install r
```

## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License

<a href="https://github.com/ryandkuster/composer/blob/master/LICENSE">Apache License Version 2.0</a>

<a><img src="https://i.imgur.com/uQ0kCRk.png" title="source: imgur.com" /></a>
