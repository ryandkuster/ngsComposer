class font_mod:
    green='\033[32m'
    lightgreen = '\033[92m'
    blue='\033[34m'
    lightblue = '\033[94m'
    cyan = '\033[36m'
    lightcyan = '\033[96m'
    yellow = '\033[93m'
    bold = '\u001b[1m'
    reset = '\u001b[0m'

conf_confirm1 = (font_mod.yellow + '\n\npaired must be True or False\n'  + font_mod.reset)
conf_confirm2 = (font_mod.yellow + '\n\nprocs must be an integer >= 1\n'  + font_mod.reset)
conf_confirm3 = (font_mod.yellow + '\n\npath to alt_dir not found\n'  + font_mod.reset)
conf_confirm4 = (font_mod.yellow + '\n\ninitial_qc must be True or False\n'  + font_mod.reset)
conf_confirm5 = (font_mod.yellow + '\n\nall_qc must be False, \'full\', or \'summary\' (quotes required)\n'  + font_mod.reset)
conf_confirm6 = (font_mod.yellow + '\n\nwalkaway must be True or False\n'  + font_mod.reset)
conf_confirm7 = (font_mod.yellow + '\n\nall_qc must be defined for walkthrough mode\n'  + font_mod.reset)
conf_confirm8 = (font_mod.yellow + '\n\nfront_trim must be an integer\n'  + font_mod.reset)
conf_confirm9 = (font_mod.yellow + '\n\nfront_trim can\'t be less than 1\n'  + font_mod.reset)
conf_confirm10 = (font_mod.yellow + '\n\nmismatch defined, conf_confirmpected \'indconf_confirm.txt\' in project directory\n' + font_mod.reset)
conf_confirm11 = (font_mod.yellow + '\n\nmismatch must be an integer >= 0\n'  + font_mod.reset)
conf_confirm12 = (font_mod.yellow + '\n\nR1_bases_ls must be a list of motifs (quotes and brackets required)\n'  + font_mod.reset)
conf_confirm13 = (font_mod.yellow + '\n\nR2_bases_ls must be a list of motifs (quotes and brackets required)\n'  + font_mod.reset)
conf_confirm14 = (font_mod.yellow + '\n\nnon_genomic must be an integer >= 1\n'  + font_mod.reset)
conf_confirm15 = (font_mod.yellow + '\n\nboth auto_trim and trim_mode variables must be defined\n'  + font_mod.reset)
conf_confirm16 = (font_mod.yellow + '\n\nauto_trim must be an integer between 0 and 42\n'  + font_mod.reset)
conf_confirm17 = (font_mod.yellow + '\n\ntrim_mode must be \'whisker\', \'quartile\', \'median\', or \'mean\' (quotes required)\n'  + font_mod.reset)
conf_confirm18 = (font_mod.yellow + '\n\nboth q_min and q_percent variables must be defined\n'  + font_mod.reset)
conf_confirm19 = (font_mod.yellow + '\n\nq_min must be an integer between 0 and 42\n'  + font_mod.reset)
conf_confirm20 = (font_mod.yellow + '\n\nq_percent must be an integer between 0 and 100\n'  + font_mod.reset)
conf_confirm21 = (font_mod.yellow + '\n\nrm_dirs must be True or False\n'  + font_mod.reset)
conf_confirm22 = (font_mod.yellow + '\n\nmin_start must be an integer\n'  + font_mod.reset)
conf_confirm23 = (font_mod.yellow + '\n\n\'adapters.txt\', min_start, and adapter_mismatch must be present\n'  + font_mod.reset)
conf_confirm24 = (font_mod.yellow + '\n\nadapter_mismatch must be an integer\n'  + font_mod.reset)
conf_confirm25 = (font_mod.yellow + '\n\nphred64 must be True or False\n'  + font_mod.reset)


is_fq1 = (font_mod.yellow + '\n\nproject directory should not contain subdirectories\n'  + font_mod.reset)


paired_vars = (font_mod.yellow + '\nunexpected paired libraries found\n' + font_mod.reset +
               '\nplease select from the following options:\n\n' +
               font_mod.lightblue + ' 1' + font_mod.reset + ' - continue treating files as single-end libraries?\n' +
               font_mod.lightblue + ' 2' + font_mod.reset + ' - continue treating files as paired-end libraries?\n' +
               font_mod.lightblue + ' 3' + font_mod.reset + ' - exit ngs-composer\n\n' +
               font_mod.lightblue + 'number selection > ' + font_mod.reset)


nucleotide_test1 = (font_mod.yellow + '\nempty list, sequence motif(s) expected\n' + font_mod.reset)
nucleotide_test2 = (font_mod.yellow + '\nbarcodes and motifs must be upper-case only\n' + font_mod.reset)
nucleotide_test3 = (font_mod.yellow + '\nbarcodes and motifs must consist of A, C, G, or T only\n' + font_mod.reset)


sep = '#' * 50
crin_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'crinoid  - qc stats' + font_mod.reset + '\n' + sep)
scal_title1 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop  - trimming' + font_mod.reset + '\n' + sep)
scal_title2 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop  - auto-trimming' + font_mod.reset + '\n' + sep)
nem_title =   ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'anemone  - demultiplexing' + font_mod.reset + '\n' + sep)
rot_title =   ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'rotifer  - motif detection' + font_mod.reset + '\n' + sep)
kril_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'krill    - filtering' + font_mod.reset + '\n' + sep)
porf_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'porifera - adapter removal' + font_mod.reset + '\n' + sep)


walkthrough1 = ('\nplease select from the following options:\n\n' +
         font_mod.lightblue + ' 1' + font_mod.reset + ' - accept results of current step and continue\n' +
         font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters and rerun current step\n' +
         font_mod.lightblue + ' 3' + font_mod.reset + ' - remove results and bypass current step\n' +
         font_mod.lightblue + ' 4' + font_mod.reset + ' - toggle walkaway mode\n' +
         font_mod.lightblue + ' 5' + font_mod.reset + ' - exit ngs-composer\n\n' +
         font_mod.lightblue + 'number selection > ' + font_mod.reset)

walkthrough2 = 'enabled'

walkthrough3 = (font_mod.yellow + 'about to be canceled' + font_mod.reset)

walkthrough4 = ('\nplease select from the following options:\n\n' +
         font_mod.lightblue + ' 1' + font_mod.reset + ' - continue with modified parameters\n' +
         font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters\n\n' +
         font_mod.lightblue + 'number selection > ' + font_mod.reset)

trim_assist1 = ('\nplease select from the following options:\n\n' +
                font_mod.lightblue + ' 1' + font_mod.reset + ' - accept results of current step and continue\n' +
                font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters and rerun current step\n' +
                font_mod.lightblue + ' 3' + font_mod.reset + ' - remove results and bypass current step\n' +
                font_mod.lightblue + ' 4' + font_mod.reset + ' - toggle walkaway mode\n' +
                font_mod.lightblue + ' 5' + font_mod.reset + ' - exit ngs-composer\n\n' +
                font_mod.lightblue + 'number selection > ' + font_mod.reset)
