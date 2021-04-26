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

if (suppressMessages(!require(ggplot2, lib.loc = .libPaths()[1]))) {
  print(paste('Installing dependency (ggplot2) into', .libPaths()[1]))
  install.packages('ggplot2', dependencies=TRUE, quiet=TRUE, repos='http://cran.us.r-project.org')
  }
suppressMessages(library(ggplot2))
