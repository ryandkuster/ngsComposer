import gzip

class file_type:
    def is_fq(filename, pairs):
        with open(filename) as f:
            i = 0
            for line in f:
                i += 1
                if i == 1:
                    if line[0] != '@':
                        return
                    else:
                        pairs = file_type.is_barcode(filename, line, pairs)
                if i == 3 and line[0] != '+':
                    return
                if i == 5:
                    if line[0] != '@':
                        return
                    else:
                        return True, pairs

    def is_gz(filename, pairs):
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
                        return True, pairs

    def is_barcode(filename, line, pairs):
            if line in pairs:
                pairs[line].append(filename)
            else:
                pairs[line] = [filename]
            return pairs
        
    def is_paired(filename):
        pass

    def phred(file):
        pass
