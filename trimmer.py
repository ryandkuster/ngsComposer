import sys
import os
import gzip

class trimmer:
    def test(input_file, b, e, output_file, gzipped):
        if gzipped == 0:
            with open(input_file) as f, \
                 open('trimmer_' + output_file, 'w') as o:
                trimmer.trim(f, b, e, o)
        if gzipped == 1:
            output_file = output_file.strip('.gz')
            with gzip.open(input_file, 'rt') as f, \
                 open('trimmer_' + output_file, 'w') as o:
                trimmer.trim(f, b, e, o)
    
    def trim(f, b, e, o):
	    i = 0
	    for line in f:
		    i += 1
		    line2or4 = i%2
		    if line2or4 == 0:
			    line = line[b:-(e+1)] + "\n"
			    o.write(line)
		    else:
			    o.write(line)

if __name__ == '__main__':
    input_file = sys.argv[1] # fastq file
    b = int(sys.argv[2]) # bases to remove from beginning of sequence
    e = int(sys.argv[3]) # bases to remove from end of sequence
    output_file = sys.argv[4]
    gzipped = 0 # change to '1' for gzipped file
    trimmer.test(input_file, b, e, output_file, gzipped)
