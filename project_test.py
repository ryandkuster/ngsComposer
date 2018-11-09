import sys
import os
import gzip

class reader:
    def barcode_reader(project_dir, config_barcodes, prefix):
        if config_barcodes == True:
            barcode_list = []
            try:
                with open(project_dir + '/' + prefix + '_barcodes.txt') as f:
                    for line in f:
                        barcode_list.append(line.rstrip())
                    return barcode_list
            except FileNotFoundError:
                sys.exit('''
                based on your configuration file your project
                directory must contain a newline-separated list
                of barcodes named ''' + prefix + '''_barcodes.txt
                    ''')

    def fastq_reader(project_dir):
        fastq_list, gz, pairs = [], [], {}
        for filename in os.listdir(project_dir):
            if filename == 'config.txt':
                pass
            elif filename == 'R1_barcodes.txt':
                pass
            elif filename == 'R2_barcodes.txt':
                pass
            else:
                try:
                    fastq_test, pairs = \
                    file_type.is_fq(project_dir + '/' + filename, pairs)
                    gz.append(0)
                except UnicodeDecodeError:
                    fastq_test, pairs = \
                    file_type.is_gz(project_dir + '/' + filename, pairs)
                    gz.append(1)
                if fastq_test is None:
                    raise TypeError
                fastq_list.append(project_dir + '/' + filename)
        return fastq_list, gz, pairs

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
                        pairs = file_type.is_paired(filename, line, pairs)
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
                if i == 1:
                    if line[0] != '@':
                        return
                    else:
                        pairs = file_type.is_paired(filename, line, pairs)
                if i == 3 and line[0] != '+':
                    return
                if i == 5:
                    if line[0] != '@':
                        return
                    else:
                        return True, pairs

    def is_paired(filename, line, pairs):
        for pos, x in enumerate(line):
            if x == ' ':
                space_pos = pos
        header = line[:space_pos]
        if header in pairs:
            pairs[header].append(filename)
        else:
            pairs[header] = [filename]
        return pairs

    def phred(file):
        pass
