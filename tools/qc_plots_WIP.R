library(ggplot2)

qscores <- read.table(file = '/home/ryan/Desktop/examples/project2/qscores.1_R2.fastq', header=F, sep=",")
names(qscores) <- gsub("V","", names(qscores))
sum <- format(round(as.numeric(sum(qscores[1,])), 1), nsmall=0, big.mark=",")
qc <- reshape(qscores, direction="long", varying=list(c(1:ncol(qscores))), v.names="count")
qc$score <- qc$time-1
colnames(qc)[colnames(qc)=="id"] <- "position"
var <- c("position","score","count")
qc <- qc[var]
qc <- qc[which(qc$count>0),]
qc$mean <- "NA"
qc$sd <- "NA"


boxplot <- ggplot(qc, aes(x = position, y = score, weight = count, group = position)) +
  geom_hline(aes(yintercept = 40), lty=1, size=1, color = "white") +
  geom_hline(aes(yintercept = 30), lty=1, size=1, color = "white") +
  geom_hline(aes(yintercept = 20), lty=1, size=1, color = "white") +
  geom_hline(aes(yintercept = 10), lty=1, size=1, color = "white") +
  geom_boxplot(fill = "white", color = "#568EA6", outlier.shape = 19, outlier.size = 1, outlier.alpha = .25) +
  scale_x_continuous(expand=c(0,0), limit=c(0, nrow(qscores)+1), breaks=seq(0, nrow(qscores)+1,by=5), sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(expand=c(0,0), limit=c(0,43), breaks=seq(0,43, by=2), sec.axis = dup_axis(name = NULL)) +
  theme_classic()+
  xlab(paste("Read Position\n(Total Number of reads = ",sum,")", sep="")) +
  ylab("Quality Scores (phred+33)") +
  theme(panel.background = element_rect(fill = "#f0f2f1"))
boxplot
suppressMessages(ggsave(filename= paste(composer_in[2], "_outliers.tiff", sep = ""), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw"))
invisible(gc())

