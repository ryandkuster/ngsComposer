import sys
import os


def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    in1 = sys.argv[1]
    front_trim = int(sys.argv[2])
    back_trim = int(sys.argv[3])
    out1 = sys.argv[4]
    proj_dir = os.path.abspath(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


def scallop_comp(front_trim, back_trim, proj_dir, in1):
    '''
    composer entry point to scallop
    '''
    out1 = proj_dir + '/' + os.path.basename(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


def scallop(in1, front_trim, back_trim, proj_dir, out1):
    '''
    trim defined base numbers from the front or from the end of reads
    '''
    if back_trim == 0:
        back_trim = None
    else:
        back_trim = - back_trim
    with open(in1) as f, open(out1, 'w') as o:
        for i, line in enumerate(f):
            i += 1
            if i % 2 == 0:
                line = line.rstrip()[front_trim:back_trim] + "\n"
                o.write(line)
            else:
                o.write(line)


def scallop_fixed():
    '''
    trim defined slice of reads, producing constant read length
    '''
    #TODO open qscores raw file
    #TODO find point in raw file where read length is not variable
    #TODO calculate median and corrresponding IQRs
    #TODO find point at which q30 is acceptable x number of times (1 for now)
    #TODO run essentially identical function as 'scallop' without negative indexing




if __name__ == '__main__':
    scallop_main()
