<a href="https://imgur.com/gYTN9Hv"><img src="https://i.imgur.com/gYTN9Hv.png" title="source: imgur.com" /></a>

# composer: a deep QC aventure!

Base-calling, error-filtering read preprocessing designed by biologists

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
- 1_R2.fastq (optional)
- barcodes_1.txt
- index.txt

***

barcodes_1.txt (tab or space delimited)
```
	A	C	G	T
A	sample1	sample5	sample6	sample10
C	sample2	sample5	sample7	sample10
G	sample3	sample5	sample8	sample10
T	sample4	sample5	sample9	sample10
```

***

index.txt
```
1_R1.fastq  barcodes_1.txt
```

***

From command line, run composer with the specified directory of your project
```bash
$ python3 composer.py <path_to_directory>
```

## License

<a href="https://github.com/ryandkuster/Pipeline/blob/master/LICENSE">Apache License Version 2.0</a>

<a href="https://imgur.com/uQ0kCRk"><img src="https://i.imgur.com/uQ0kCRk.png" title="source: imgur.com" /></a>
