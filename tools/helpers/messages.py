class font:
    green='\033[32m'
    lightgreen = '\033[92m'
    blue='\033[34m'
    lightblue = '\033[94m'
    cyan = '\033[36m'
    lightcyan = '\033[96m'
    yellow = '\033[93m'
    bold = '\u001b[1m'
    reset = '\u001b[0m'


initialize1 = (font.yellow + '\n\nproject directory not found\n' + font.reset)


conf_confirm1 = (font.yellow + '\n\npaired must be True or False\n' + font.reset)
conf_confirm2 = (font.yellow + '\n\nprocs must be an integer >= 1\n' + font.reset)
conf_confirm3 = (font.yellow + '\n\npath to alt_dir not found\n' + font.reset)
conf_confirm4 = (font.yellow + '\n\ninitial_qc must be True or False\n' + font.reset)
conf_confirm5 = (font.yellow + '\n\nall_qc must be False, \'full\', or \'summary\' (quotes required)\n' + font.reset)
conf_confirm6 = (font.yellow + '\n\nwalkaway must be True or False\n' + font.reset)
conf_confirm7 = (font.yellow + '\n\nall_qc must be defined for walkthrough mode\n' + font.reset)
conf_confirm8 = (font.yellow + '\n\nfront_trim must be an integer\n' + font.reset)
conf_confirm9 = (font.yellow + '\n\nfront_trim can\'t be less than 1\n' + font.reset)
conf_confirm10 = (font.yellow + '\n\nmismatch defined, \'index.txt\' expected in project directory\n' + font.reset)
conf_confirm11 = (font.yellow + '\n\nmismatch must be an integer >= 0\n' + font.reset)
conf_confirm12 = (font.yellow + '\n\nR1_bases_ls must be a list of motifs (quotes and brackets required)\n' + font.reset)
conf_confirm13 = (font.yellow + '\n\nR2_bases_ls must be a list of motifs (quotes and brackets required)\n' + font.reset)
conf_confirm14 = (font.yellow + '\n\nnon_genomic must be an integer >= 1\n' + font.reset)
conf_confirm15 = (font.yellow + '\n\nboth auto_trim and trim_mode variables must be defined\n' + font.reset)
conf_confirm16 = (font.yellow + '\n\nauto_trim must be an integer between 0 and 42\n' + font.reset)
conf_confirm17 = (font.yellow + '\n\ntrim_mode must be \'whisker\', \'quartile\', \'median\', or \'mean\' (quotes required)\n' + font.reset)
conf_confirm18 = (font.yellow + '\n\nboth q_min and q_percent variables must be defined\n' + font.reset)
conf_confirm19 = (font.yellow + '\n\nq_min must be an integer between 0 and 42\n' + font.reset)
conf_confirm20 = (font.yellow + '\n\nq_percent must be an integer between 0 and 100\n' + font.reset)
conf_confirm21 = (font.yellow + '\n\nrm_dirs must be True or False\n' + font.reset)
conf_confirm22 = (font.yellow + '\n\nmin_start must be an integer\n' + font.reset)
conf_confirm23 = (font.yellow + '\n\n\'adapters.R1.txt\' (adapters and barcodes in corresponding reverse end of read), \'adapters.R2.txt\' and adapter_match must be present\n' + font.reset)
conf_confirm24 = (font.yellow + '\n\nadapter_match must be an integer >= 10 for pipeline\n' + font.reset)
conf_confirm25 = (font.yellow + '\n\nphred64 must be True or False\n' + font.reset)


r_packages1 = (font.yellow + '\n\nplease install latest version of R\n' + font.reset)


fastq_test1 = (font.yellow + ' was not expected in project directory\n' + font.reset)


adapters_test1 = (font.yellow + '\nbarcodes not found in adapters file(s)\n' + font.reset +
                 '\nplease select from the following options:\n\n' +
                 font.lightblue + ' 1' + font.reset + ' - continue without including barcodes in adapters\n' +
                 font.lightblue + ' 2' + font.reset + ' - exit\n' +
                 font.lightblue + 'number selection > ' + font.reset)


dir_size1 = (font.yellow + '\n\nan estimated ' + font.reset)
dir_size2 = (font.yellow +' bytes are required to process, consider rm_transit or alt_dir variables\n' + font.reset +
             '\nplease select from the following options:\n\n' +
             font.lightblue + ' 1' + font.reset + ' - continue\n' +
             font.lightblue + ' 2' + font.reset + ' - exit\n\n' +
             font.lightblue + 'number selection > ' + font.reset)


nucleotide_test1 = (font.yellow + '\nempty list, sequence motif(s) expected\n' + font.reset)
nucleotide_test2 = (font.yellow + '\nbarcodes and motifs must be upper-case only\n' + font.reset)
nucleotide_test3 = (font.yellow + '\nbarcodes and motifs must consist of A, C, G, or T only\n' + font.reset)


dec1 = ('\n' + '#' * 50 + '\n' + ' ' * 15  + font.lightblue)
dec2 = (font.reset + '\n' + '#' * 50)
crin_title = (dec1 + 'crinoid  - qc stats' + dec2)
scal_title1 = (dec1 + 'scallop  - trimming' + dec2)
scal_title2 = (dec1 + 'scallop  - auto-trimming' + dec2)
nem_title = (dec1 + 'anemone  - demultiplexing' + dec2)
rot_title = (dec1 + 'rotifer  - motif detection' + dec2)
kril_title = (dec1 + 'krill    - filtering' + dec2)
porf_title = (dec1 + 'porifera - adapter removal' + dec2)


walkthrough1 = ('\nplease select from the following options:\n\n' +
                font.lightblue + ' 1' + font.reset + ' - accept results of current step and continue\n' +
                font.lightblue + ' 2' + font.reset + ' - modify parameters and rerun current step\n' +
                font.lightblue + ' 3' + font.reset + ' - remove results and bypass current step\n' +
                font.lightblue + ' 4' + font.reset + ' - toggle walkaway mode\n' +
                font.lightblue + ' 5' + font.reset + ' - exit ngsComposer\n\n' +
                font.lightblue + 'number selection > ' + font.reset)
walkthrough2 = 'enabled'
walkthrough3 = (font.yellow + 'about to be canceled' + font.reset)
walkthrough4 = ('\nplease select from the following options:\n\n' +
                font.lightblue + ' 1' + font.reset + ' - continue with modified parameters\n' +
                font.lightblue + ' 2' + font.reset + ' - modify parameters\n\n' +
                font.lightblue + 'number selection > ' + font.reset)
