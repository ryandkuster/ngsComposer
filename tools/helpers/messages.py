sep = '#' * 50

crin_title = (sep + '\n' + ' ' * 15 + 'crinoid - qc stats' + '\n' + sep)
scal_title1 = (sep + '\n' + ' ' * 15 + 'scallop - front-trimming' + '\n' + sep)
scal_title2 = (sep + '\n' + ' ' * 15 + 'scallop - end-trimming' + '\n' + sep)
nem_title = (sep + '\n' + ' ' * 15 + 'anemone - demultiplexing' + '\n' + sep)
rot_title = (sep + '\n' + ' ' * 15 + 'rotifer - motif detection' + '\n' + sep)
kril_title = (sep + '\n' + ' ' * 15 + 'krill   - filtering' + '\n' + sep)


anemone_qc = ('\ncheck qc file(s) in demultiplexed directory\n\n' +
              'continue without changes to configuration settings? (y/n)\n')
anemone_up = ('\nupdate mismatch value and rerun anemone demultiplexing ' +
              'step? (y/n)\n')
anemone_in = ('\nnew mismatch value? (int)\n')


rotifer_qc = ('\ncheck qc file(s) in parsed directory\n\n' +
              'continue without changes to configuration settings? (y/n)\n')
rotifer_up = ('\nupdate motif lists and rerun rotifer parsing step? (y/n)\n')
rotifer_in1 = ('\nnew R1_bases_ls? (space separated list, or press enter to ' +
               'continue the pipeline disregarding these motifs\n')
rotifer_in2 = ('\nnew R2_bases_ls? (space separated list, or press enter to ' +
               'continue the pipeline disregarding these motifs\n')


scallop_qc = ('\ncheck qc file(s) in end_trimmed directory\n\n' +
              'continue without changes to configuration settings? (y/n)\n')
scallop_up = ('\nupdate threshold value and rerun end-trimming step? (y/n)\n')
scallop_in = ('\nnew end_trim value? (int, or press enter to continue the ' +
              'pipeline disregarding end-trimming\n')





confirm = ['Y', 'y', 'Yes', 'yes', 'YES']
