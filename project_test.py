import sys
import os
import gzip

class reader:
    def config_reader(project_dir):
        try:
            with open(project_dir + 'config.txt') as f:
                for line in f:
                    print(line)
        except FileNotFoundError:
            print('''
                your project directory must contain a 
                configuration file named config.txt
                ''')

    def barcode_reader(project_dir, config_barcodes, prefix):
        if config_barcodes == True:
            try:
                with open(project_dir + prefix + '_barcodes.txt') as f:
                    for line in f:
                        print(line)
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
                    file_type.is_fq(os.path.abspath(filename), pairs)
                    gz.append(0)
                except UnicodeDecodeError:
                    fastq_test, pairs = \
                    file_type.is_gz(os.path.abspath(filename), pairs)
                    gz.append(1)
                if fastq_test is None:
                    raise TypeError
                fastq_list.append(os.path.abspath(filename))
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
        end = line[space_pos+1]
        if header in pairs:
            pairs[header].append(filename)
        else:
            pairs[header] = [filename]
        return pairs

    def phred(file):
        pass
