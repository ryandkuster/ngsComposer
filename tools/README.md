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

## Crinoid

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 or R2 fastq file|
|-o|O|the full path to output directory (optional)|
|-a|A|create visualizations on existing qscore and nucleotide raw data (optional, requires -r1 input)|
|-s|S|type True for phred 64 samples (optional, default phred 33)|

Example:
```bash
$ python3 crinoid.py -r1 1_R1.fastq
```

## Scallop

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 or R2 fastq file|
|-f|F|number of bases to remove from beginning of read (integer)|
|-b|B|final position to keep within a read (integer)|
|-o|O|the full path to output directory (optional)|

Example:
```bash
$ python3 scallop.py -r1 1_R1.fastq -f 6
```

## Anemone

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 fastq file|
|-r2|R2|the full or relative path to R2 fastq file (optional)|
|-c|C|the full or relative path to barcodes index file|
|-m|M|mismatch value for barcode hamming distance (integer)|
|-o|O|the full path to output directory (optional)|


Example:
```bash
$ python3 anemone.py -r1 1_R1.fastq -r2 1_R2.fastq -m 0 -c barcodes.txt
```

## Rotifer

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 fastq file|
|-r2|R2|the full or relative path to R2 fastq file (optional)|
|-m1|M1|space separated list of motifs expected in R1 (no brackets)|
|-m2|M2|space separated list of motifs expected in R2 (no brackets)|
|-n|N|number of non-genomic bases to remove from read start (integer)|
|-o|O|the full path to output directory (optional)|

Example:
```bash
$ python3 rotifer.py -r1 1_R1.fastq -r2 1_R2.fastq -m1 TCT TCC -m2 TCT TCC
```

## Krill

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 fastq file|
|-r2|R2|the full or relative path to R2 fastq file (optional)|
|-q|Q|the minimum qscore required at a position (integer)|
|-p|P|the percent frequency minimum qscore must occur per read (integer)|
|-o|O|the full path to output directory (optional)|

Example:
```bash
$ python3 krill.py -r1 1_R1.fastq -r2 1_R2.fastq -q 30 -p 95
```

## Porifera

|Variable|Alternate Variable|Usage|
|:--|:--|:--|
|-r1|R1|the full or relative path to R1 or R2 fastq file|
|-a|A|the full or relative path to adapter sequences file|
|-n|N|number of bases from end of read to begin adapter alignment|
|-m|M|mismatch value for adapter hamming distance (integer)|
|-o|O|the full path to output directory (optional)|


Example:
```bash
$ python3 porifera.py -r1 1_R1.fastq -a adapters.txt -n 18 -m 3
```

<a href="https://imgur.com/uQ0kCRk"><img src="https://i.imgur.com/uQ0kCRk.png" title="source: imgur.com" /></a>
