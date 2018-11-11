import sys
import os

def trimmer(trim, fastq):
	   with open(fastq) as f, open('trimmed_' + os.path.basename(fastq), 'w') as o:
	    for i, line in enumerate(f):
		    if i%2 == 0:
			    line = line[b:-(e+1)] + "\n"
			    o.write(line)
		    else:
			    o.write(line)

if __name__ == '__main__':
    input_file = sys.argv[1] # fastq file
    front_trim = int(sys.argv[2]) # bases to remove from beginning of sequence
    back_trim = int(sys.argv[3]) # bases to remove from end of sequence
    trimmer(front_trim, back_trim, input_file)
