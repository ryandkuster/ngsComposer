import sys
import os

def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    input1 = sys.argv[1] # name of input fastq file
    front_trim = int(sys.argv[2]) # bases to remove from beginning of sequence
    back_trim = int(sys.argv[3]) # bases to remove from end of sequence
    output1 = int(sys.argv[4]) # name of output file
    project_dir = os.path.abspath(input1)
    scallop(input1, front_trim, back_trim, project_dir, output1)


def scallop(input1, front_trim, back_trim, project_dir, output1):
    '''
    composer entry point to scallop
    '''
    with open(input1) as f, \
        open(project_dir + '/' + output1, 'w') as o:
        for i, line in enumerate(f):
            i += 1
            if i%2 == 0:
                line = line[front_trim:-(back_trim+1)] + "\n"
                o.write(line)
            else:
                o.write(line)


def scallop_pipeline(front_trim, back_trim, project_dir, input1):
    output1 = os.path.basename(input1)
    scallop(input1, front_trim, back_trim, project_dir, output1)


if __name__ == '__main__':
    scallop_main()