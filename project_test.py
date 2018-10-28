import gzip

class file_type:
    def is_fastq(fastq):
        with open(fastq) as f:
            i = 0
            for line in f:
                i += 1
                if i == 1 and line[0] != '@':
                    return
                if i == 3 and line[0] != '+':
                    return
                if i == 5:
                    if line[0] != '@':
                        return
                    else:
                        return True

    def is_gzipped(filename):
        with gzip.open(filename, 'rt') as f:
            i = 0
            for line in f:
                i += 1
                if i == 1 and line[0] != '@':
                    return
                if i == 3 and line[0] != '+':
                    return
                if i == 5:
                    if line[0] != '@':
                        return
                    else:
                        return True

    def is_barcode(file):
        pass
        
class fastq_type:
    def phred(file):
        pass
