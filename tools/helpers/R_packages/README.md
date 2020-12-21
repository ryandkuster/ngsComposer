### Optional folder for writing necessary R packages  
Typically, ngsComposer tools that call the Rscript 'qc_plots.R' require no direction to the R packages library (uses default lib.loc). In some instances, it may be useful to manually provide a 'r_dir' variable:  
- ngsComposer pipeline (composer.py): use 'r_dir' to define the R packages directory in the 'conf.py' file (e.g., r_dir = '.../ngsComposer/tools/helpers/R_packages'  
- QC plotting tool (crinoid.py): use the -g argument (e.g., -g '.../ngsComposer/tools/helpers/R_packages'
