import sys
import os


def scallop_main():
    '''
    standalone, command line entry point to scallop using stdin
    '''
    in1 = sys.argv[1]
    front_trim = int(sys.argv[2])
    back_trim = int(sys.argv[3])
    out1 = int(sys.argv[4])
    proj_dir = os.path.abspath(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


def scallop(in1, front_trim, back_trim, proj_dir, out1):
    '''
    composer entry point to scallop
    '''
    if back_trim == 0:
        back_trim = None
    else:
        back_trim = - back_trim
    with open(in1) as f, open(proj_dir + '/' + out1, 'w') as o:
        for i, line in enumerate(f):
            i += 1
            if i % 2 == 0:
                line = line.rstrip()[front_trim:back_trim] + "\n"
                o.write(line)
            else:
                o.write(line)


def scallop_comp(front_trim, back_trim, proj_dir, in1):
    out1 = os.path.basename(in1)
    scallop(in1, front_trim, back_trim, proj_dir, out1)


if __name__ == '__main__':
    scallop_main()
