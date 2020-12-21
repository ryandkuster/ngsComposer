args <- commandArgs(trailingOnly = TRUE)

if (args[1] == 'None') {
  if (suppressMessages(!require(ggplot2))) install.packages('ggplot2', dependencies=TRUE, repos='http://cran.us.r-project.org')
  suppressMessages(library(ggplot2))
} else {
  if (suppressMessages(!require(ggplot2, lib.loc = args[1]))) install.packages('ggplot2', dependencies=TRUE, repos='http://cran.us.r-project.org', lib = args[1])
  suppressMessages(library(ggplot2))
}
