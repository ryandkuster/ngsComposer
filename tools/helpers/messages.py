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


conf_end = (font_mod.yellow + '\nplease use either manual or auto end trim' + font_mod.reset)
q_vars = (font_mod.yellow + '\nboth q_min and q_percent variables must be defined' + font_mod.reset)
nucs1 = (font_mod.yellow + '\nbarcodes and motifs must be upper-case only' + font_mod.reset)
nucs2 = (font_mod.yellow + '\nbarcodes and motifs must consist of A, C, G, or T only' + font_mod.reset)


sep = '#' * 50
crin_title = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'crinoid - qc stats' + font_mod.reset + '\n' + sep)
scal_title1 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop - trimming' + font_mod.reset + '\n' + sep)
scal_title2 = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'scallop - auto-trimming' + font_mod.reset + '\n' + sep)
nem_title = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'anemone - demultiplexing' + font_mod.reset + '\n' + sep)
rot_title = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'rotifer - motif detection' + font_mod.reset + '\n' + sep)
kril_title = ('\n' + sep + '\n' + ' ' * 15 + font_mod.lightblue + 'krill   - filtering' + font_mod.reset + '\n' + sep)


walk1 = ('\nplease select from the following options:\n\n' +
         font_mod.lightblue + ' 1' + font_mod.reset + ' - accept results of current step and continue\n' +
         font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters and rerun current step\n' +
         font_mod.lightblue + ' 3' + font_mod.reset + ' - remove results and bypass current step\n' +
         font_mod.lightblue + ' 4' + font_mod.reset + ' - toggle walkaway mode\n' +
         font_mod.lightblue + ' 5' + font_mod.reset + ' - exit ngs-composer\n\n' +
         font_mod.lightblue + 'number selection > ' + font_mod.reset)

walk2 = 'enabled'

walk3 = (font_mod.yellow + 'about to be canceled' + font_mod.reset)

trim_assist1 = ('\nplease select from the following options:\n\n' +
                font_mod.lightblue + ' 1' + font_mod.reset + ' - accept results of current step and continue\n' +
                font_mod.lightblue + ' 2' + font_mod.reset + ' - modify parameters and rerun current step\n' +
                font_mod.lightblue + ' 3' + font_mod.reset + ' - remove results and bypass current step\n' +
                font_mod.lightblue + ' 4' + font_mod.reset + ' - toggle walkaway mode\n' +
                font_mod.lightblue + ' 5' + font_mod.reset + ' - exit ngs-composer\n\n' +
                font_mod.lightblue + 'number selection > ' + font_mod.reset)
