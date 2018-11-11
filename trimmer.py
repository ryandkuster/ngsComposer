import sys
import os

def trimmer(front_trim, back_trim, project_dir, fastq):
    with open(fastq) as f, \
        open(project_dir + '/trimmed_' + os.path.basename(fastq), 'w') as o:
        for i, line in enumerate(f):
            i += 1
            if i%2 == 0:
                line = line[front_trim:-(back_trim+1)] + "\n"
                o.write(line)
            else:
                o.write(line)

if __name__ == '__main__':
    fastq = sys.argv[1] # fastq file
    front_trim = int(sys.argv[2]) # bases to remove from beginning of sequence
    back_trim = int(sys.argv[3]) # bases to remove from end of sequence
    project_dir = os.path.abspath(fastq)
    trimmer(front_trim, back_trim, project_dir, fastq)
