import sys

def overhang():
    fastq = sys.argv[1] # name of input file
    out1 = sys.argv[2] # name of output file for reads starting with desired sequence(s)
    out2 = sys.argv[3] # name of output file for reads not starting with desired sequence(s)
    with open(fastq) as s:
        line1 = s.readline()
        line2 = s.readline()
    seq_len = len(line2) - 1
    with open(fastq) as f, open(str(seq_len) + "_" + out1, 'a') as o1, open(str(seq_len) + "_" + out2, 'a') as o2:
        y, entry = 0, ""
        for line in f:
            y += 1
            if y == 2:
                if line.startswith("TCC") or line.startswith("TCT"):
                    overhang = 1
                else:
                    overhang = 0
            entry = entry + line
            if y == 4 and overhang == 1:
                o1.write(entry)
                y, overhang, entry = 0, 0, ""
            if y == 4 and overhang == 0:
                o2.write(entry)
                y, overhang, entry = 0, 0, ""

if __name__ == "__main__":
    overhang()
