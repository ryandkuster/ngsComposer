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

ex1 = (font_mod.yellow + 'paired must be True or False'  + font_mod.reset)
ex2 = (font_mod.yellow + 'procs must be an integer >= 1'  + font_mod.reset)
ex3 = (font_mod.yellow + 'path to alt_dir not found'  + font_mod.reset)
ex4 = (font_mod.yellow + 'initial_qc must be True or False'  + font_mod.reset)
ex5 = (font_mod.yellow + 'all_qc must be False, \'full\', or \'summary\' (quotes required)'  + font_mod.reset)
ex6 = (font_mod.yellow + 'walkaway must be True or False'  + font_mod.reset)
ex7 = (font_mod.yellow + 'all_qc must be defined for walkthrough mode'  + font_mod.reset)
ex8 = (font_mod.yellow + 'front_trim must be an integer'  + font_mod.reset)
ex9 = (font_mod.yellow + 'front_trim can\'t be less than 1'  + font_mod.reset)
ex10 = (font_mod.yellow + 'mismatch defined, expected \'index.txt\' in project directory' + font_mod.reset)
ex11 = (font_mod.yellow + 'mismatch must be an integer >= 0'  + font_mod.reset)
ex12 = (font_mod.yellow + 'R1_bases_ls must be a list of motifs (quotes and brackets required)'  + font_mod.reset)
ex13 = (font_mod.yellow + 'R2_bases_ls must be a list of motifs (quotes and brackets required)'  + font_mod.reset)
ex14 = (font_mod.yellow + 'non_genomic must be an integer >= 1'  + font_mod.reset)
ex15 = (font_mod.yellow + ''  + font_mod.reset)
ex16 = (font_mod.yellow + ''  + font_mod.reset)
ex17 = (font_mod.yellow + ''  + font_mod.reset)
ex18 = (font_mod.yellow + ''  + font_mod.reset)
ex19 = (font_mod.yellow + ''  + font_mod.reset)
ex20 = (font_mod.yellow + ''  + font_mod.reset)




conf_end = (font_mod.yellow + '\nplease use either manual or auto end trim' + font_mod.reset)
q_vars = (font_mod.yellow + '\nboth q_min and q_percent variables must be defined' + font_mod.reset)
trim_vars = (font_mod.yellow + '\nboth auto_trim and trim_mode variables must be defined' + font_mod.reset)

paired_vars = (font_mod.yellow + '\nunexpected paired libraries found\n' + font_mod.reset +
               '\nplease select from the following options:\n\n' +
               font_mod.lightblue + ' 1' + font_mod.reset + ' - continue treating files as single-end libraries?\n' +
               font_mod.lightblue + ' 2' + font_mod.reset + ' - continue treating files as paired-end libraries?\n' +
               font_mod.lightblue + ' 3' + font_mod.reset + ' - exit ngs-composer\n\n' +
               font_mod.lightblue + 'number selection > ' + font_mod.reset)

nucs1 = (font_mod.yellow + '\nbarcodes and motifs must be upper-case only' + font_mod.reset)
nucs2 = (font_mod.yellow + '\nbarcodes and motifs must consist of A, C, G, or T only' + font_mod.reset)


sep = '#' * 50
crin_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'crinoid  - qc stats' + font_mod.reset + '\n' + sep)
scal_title1 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop  - trimming' + font_mod.reset + '\n' + sep)
scal_title2 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop  - auto-trimming' + font_mod.reset + '\n' + sep)
nem_title =   ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'anemone  - demultiplexing' + font_mod.reset + '\n' + sep)
rot_title =   ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'rotifer  - motif detection' + font_mod.reset + '\n' + sep)
kril_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'krill    - filtering' + font_mod.reset + '\n' + sep)
porf_title =  ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'porifera - adapter removal' + font_mod.reset + '\n' + sep)


walk1 = ('\nplease select from the following options:\n\n' +
         font_mod.lightblue + ' 1' + font_mod.reset + ' - accept results of current step and continue\n' +
         font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters and rerun current step\n' +
         font_mod.lightblue + ' 3' + font_mod.reset + ' - remove results and bypass current step\n' +
         font_mod.lightblue + ' 4' + font_mod.reset + ' - toggle walkaway mode\n' +
         font_mod.lightblue + ' 5' + font_mod.reset + ' - exit ngs-composer\n\n' +
         font_mod.lightblue + 'number selection > ' + font_mod.reset)

walk2 = 'enabled'

walk3 = (font_mod.yellow + 'about to be canceled' + font_mod.reset)

walk4 = ('\nplease select from the following options:\n\n' +
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
