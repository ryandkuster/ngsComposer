library(ggplot2)
library(reshape)

x <- commandArgs(trailingOnly = TRUE)
# qscore_files = dir(path=x, pattern="qscores*")
# for (i in qscore_files){
#   write(paste(x, "/", i, sep = ""), file = "/home/ryan/testing/composer_test/output.txt",
#         append = TRUE, sep = " ")
# }

qscore_files = dir(path=x, pattern="qscores*")
for (i in qscore_files){
  qscore_path <- paste(x, "/", i, sep = "")
  qscores <- read.table(qscore_path, header=F, sep=",")
  names(qscores) <- gsub("V","", names(qscores))
  sum <- format(round(as.numeric(sum(qscores[1,])), 1), nsmall=0, big.mark=",")
  xyz <- stack(qscores)
  for (i in levels(xyz$ind)){
    j = strtoi(i) - 1
    print(j)
    levels(xyz$ind)[levels(xyz$ind) == i] <- toString(j)
  }
  
  qc <- reshape(qscores, direction="long", varying=list(c(1:ncol(qscores))), v.names="count")
  qc$score <- qc$time-1
  colnames(qc)[colnames(qc)=="id"] <- "position"
  var <- c("position","score","count")
  qc <- qc[var]
  boxplot <- ggplot(qc[which(qc$count>0),], aes(x = position, y=score, weight=count, group=position), stat='identity')+
    geom_boxplot(fill="white", colour="cornflowerblue", lwd = 0.25,
                 outlier.alpha = 0.2, alpha=0.2, outlier.colour="tomato", outlier.size=1.25)+
    scale_x_continuous(expand=c(0,0), limit=c(0,nrow(qscores)+1), breaks=seq(0,nrow(qscores)+1,by=5),
                       sec.axis = dup_axis(name = NULL)) +
    scale_y_continuous(expand=c(0,0), limit=c(0,41), breaks=seq(0,41, by=2), sec.axis = dup_axis(name = NULL)) +
    theme_classic()+
    xlab(paste("Read Position (Total Number of reads = ",sum,")", sep="")) +
    ylab("Quality Scores (phred+33)")
  ggsave(filename= paste(qscore_path, ".tiff"), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw")
  gc()
  
}

# nuc_files = dir(path=x, pattern="nucleotides*")
# for (i in nuc_files){
#   write(paste(x, "/", i, sep = ""), file = "/home/ryan/testing/composer_test/output.txt",
#         append = TRUE, sep = " ")
# }