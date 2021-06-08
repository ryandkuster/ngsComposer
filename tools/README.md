# ngsComposer standalone tools

A guide to using the individual tools in the ngsComposer pipeline

See <a href="https://github.com/ryandkuster/composer/blob/master/README.md">main README page</a> for installation, pipeline usage, troubleshooting, and license information.

# Contents
- [Crinoid - QC stats](#crinoid)
- [Scallop - trimming](#scallop)
- [Anemone - demultiplexing](#anemone)
- [Rotifer - motif parsing](#rotifer)
- [Krill - filtering](#krill)
- [Porifera - adapter removal](#porifera)

## Crinoid - nucleotide and Q score summary visualizations

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 or R2 fastq file|
|-o|the full path to output directory (optional)|
|-a|create visualizations on existing qscore and nucleotide raw data (optional, requires -r1 input)|
|-s|type True for phred 64 samples (optional, default phred 33)|
|-t|number of subprocesses (optional, default 1)|

Example:
```bash
$ python3 crinoid.py -r1 1_R1.fastq
```

Crinoid creates three visualizations spanning the entire read length. The first shows the per-base nucleotide frequency. There are two qscore images, both display the per-base Phred scores as a boxplot and one (labelled with '_outliers') also displays the locations of outliers (orange dots) as well as the mean and standard error (gray diamonds and lines).

## Scallop - end-trimming

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 or R2 fastq file|
|-f|number of bases to remove from beginning of read (integer)|
|-b|final position to keep within a read (integer)|
|-e|end-trim where entire window >= score (integer)|
|-w|use with 'e', size of window consisting of >= e (integer)|
|-l|use with 'e' & 'w', minimum read length to keep (integer)|
|-o|the full path to output directory (optional)|

Example:
```bash
$ python3 scallop.py -r1 1_R1.fastq -f 6
```

or

```bash
$ python3 scallop.py -r1 1_R1.fastq -w 10 -e 30 -l 50
```

The output files are automatically named with "trimmed" prefix (e.g. "trimmed.1_R1.fastq")

## Anemone - demultiplexing of single-end or paired-end barcoded libraries

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 fastq file|
|-r2|the full or relative path to R2 fastq file (optional)|
|-c|the full or relative path to barcodes index file|
|-m|mismatch value for barcode hamming distance (integer)|
|-o|the full path to output directory (optional)|


Example:
```bash
$ python3 anemone.py -r1 1_R1.fastq -r2 1_R2.fastq -m 0 -c barcodes_1.txt
```

The barcodes file is a tab or space delimited file with no spaces in sample names (or, copy directly from your favorite spreadsheet program into a text file). Forward barcodes begin each row and reverse barcodes begin each column with the desired sample names indicated in the interior of the matrix. For example, the following would be required for a dual-indexed library:

Example barcodes file:

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


## Rotifer - retain only reads beginning with expected RE cut site motifs

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 fastq file|
|-r2|the full or relative path to R2 fastq file (optional)|
|-m1|space separated list of motifs expected in R1 (no brackets)|
|-m2|space separated list of motifs expected in R2 (no brackets)|
|-n|number of non-genomic bases to remove from read start (integer)|
|-o|the full path to output directory (optional)|

Example:
```bash
$ python3 rotifer.py -r1 1_R1.fastq -r2 1_R2.fastq -m1 TCT TCC -m2 TCT TCC
```

In the event that input data is paired-end, output files will indicate when pairing has been retained with the "pe" prefix. Reads that have no pair will have the "se" prefix. In the instance that a pe and se naming scheme are both applied to a file name, the leftmost prefix will indicate if a file is paired or single/unpaired (e.g. "**se**.pe.1_R1.fastq" will be unpaired).

## Krill - quality threshold filtering of reads with a percent of bases at or above a defined Q score

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 fastq file|
|-r2|the full or relative path to R2 fastq file (optional)|
|-q|the minimum qscore required at a position (integer)|
|-p|the percent frequency minimum qscore must occur per read (integer)|
|-o|the full path to output directory (optional)|

Example:
```bash
$ python3 krill.py -r1 1_R1.fastq -r2 1_R2.fastq -q 30 -p 95
```

As with Krill, paired-end output files will indicate when pairing has been retained with the "pe" prefix. Reads that have no pair will have the "se" prefix. In the instance that a pe and se naming scheme are both applied to a file name, the leftmost prefix will indicate if a file is paired or single/unpaired (e.g. "**se**.pe.1_R1.fastq" will be unpaired).

## Porifera - adapter removal

### Note: this tool is still being optimized to detect adapters with minimal overlap, but works well to detect adapters in general

|Variable|Usage|
|:--|:--|
|-r1|the full or relative path to R1 fastq file|
|-r2|the full or relative path to R2 fastq file|
|-a1|the full or relative path to adapter sequences file|
|-a2|the full or relative path to adapter sequences file|
|-k|size of adapter k-mer used in search (default 8)|
|-m|minimum matching value for to accept adapter match (integer, default 12)|
|-r|walk adapter length by k-mer this many rounds before skipping non-matched read|
|-o|the full path to output directory (optional)|


Example:
```bash
$ python3 porifera.py -r1 1_R1.fastq -a1 adapters.txt -m 12 -k 8 -r 1 
```

Example adapter file:

**adapters.txt**
```
CGCTCAGTTC
TATCTGACCT
ATATGAGACG
CTTATGGAAT
TAATCTCGTC
GCGCGATGTT
AGAGCACTAG
TGCCTTGATC
CTACTCAGTC
TCGTCTGACT
GAACATACGG
CCTATGACTC
TAATGGCAAG
GTGCCGCTTC
CGGCAATGGA
GCCGTAACCG
AACCATTCTC 
```
*Each of the above sample adapters is presented in 5' to 3' orientation and may have restriction motifs added if desired. Porifera.py creates all reverse-complements before alignment*

<a href="https://imgur.com/uQ0kCRk"><img src="https://i.imgur.com/uQ0kCRk.png" title="source: imgur.com" /></a>
