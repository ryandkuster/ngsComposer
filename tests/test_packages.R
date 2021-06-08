args <- commandArgs(trailingOnly = TRUE)

r_version = R.Version()

if (as.numeric(r_version[[6]]) > 3) {
  invisible()
} else {
  if (as.numeric(r_version[[6]]) == 3 & as.numeric(r_version[[7]]) > 5) {
    invisible()
  } else {
    r_v = paste(r_version[[6]], '.', r_version[[7]], sep='')
    print(paste('your version of R, ', r_v, ', needs to be 3.5 or higher for ngsComposer crinoid.py and composer.py usage', sep=''))
    quit(status=1)}
}

.libPaths(args[1])

if (suppressMessages(!require(ggplot2))) {
   print('attempting to install ggplot2')
   install.packages('ggplot2', quiet=TRUE, repos='http://cran.us.r-project.org')
   }
 suppressMessages(library(ggplot2))

if (suppressMessages(!require(quantreg))) {
   print('attempting to install quantreg')
   install.packages('quantreg', quiet=TRUE, repos='http://cran.us.r-project.org')
   }
 suppressMessages(library(quantreg))
