sep = '#' * 50

nem_title = (sep + '\n' + ' ' * 22 + 'anemone' + '\n' + sep)

anemone_qc = (nem_title + '\ncheck qc file(s) in demultiplexed directory\n\n' +
        'continue without changes to configuration settings? (y/n)\n')

anemone_up = ('\nupdate mismatch value and rerun anemone demultiplexing ' +
        'step? (y/n)\n')

anemone_in = ('\nnew mismatch value? (int)\n')

rot_title = (sep + '\n' + ' ' * 22 + 'rotifer' + '\n' + sep)

rotifer_qc = (rot_title + '\ncheck qc file(s) in parsed directory\n\n' +
        'continue without changes to configuration settings? (y/n)\n')

rotifer_up = ('\nupdate motif lists and rerun rotifer parsing step? (y/n)\n')

rotifer_in1 = ('\nnew R1_bases_ls? (space separated list, or press enter to ' +
        'continue the pipeline disregarding these motifs\n')

rotifer_in2 = ('\nnew R2_bases_ls? (space separated list, or press enter to ' +
        'continue the pipeline disregarding these motifs\n')

confirm = ['Y', 'y', 'Yes', 'yes', 'YES']
