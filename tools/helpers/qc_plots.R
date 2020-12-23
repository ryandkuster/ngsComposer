options(warn=-1)

composer_in <- commandArgs(trailingOnly = TRUE)

.libPaths(composer_in[3])
suppressMessages(library(ggplot2))

qscores <- read.table(composer_in[2], header=F, sep=",")

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

if (nrow(qscores) <= 250) {
  scale = 5
} else {
  scale = round(nrow(qscores)/40,-1)
}


# create qscore distribution visualization
boxplot <- ggplot(qc, aes(x = position, y = score, weight = count, group = position)) +
  geom_hline(aes(yintercept = 40), lty=1, alpha=0.05, size=0.35) +
  geom_hline(aes(yintercept = 30), lty=1, alpha=0.05, size=0.35)+
  geom_hline(aes(yintercept = 20), lty=1, alpha=0.05, size=0.35)+
  geom_hline(aes(yintercept = 10), lty=1, alpha=0.05, size=0.35) +
  geom_boxplot(fill = "white", color = "cornflowerblue", outlier.alpha = 0.25, alpha=0.25, outlier.colour="tomato", outlier.size=1.25) +
  scale_x_continuous(expand=c(0,0), limit=c(0, nrow(qscores)+1), breaks=seq(0, nrow(qscores)+1,by=scale), sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(expand=c(0,0), limit=c(0,43), breaks=seq(0,43, by=2), sec.axis = dup_axis(name = NULL)) +
  theme_classic()+
  xlab(paste("Read Position\n(Total Number of reads = ",sum,")", sep="")) +
  ylab("Quality Scores (phred+33)")
suppressMessages(ggsave(filename= paste(composer_in[2], ".tiff", sep = ""), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw"))
invisible(gc())


# create nucleotide distribution visualization
nucs <- read.table(composer_in[1], header=F, sep=",")

names(nucs) <- c("A", "C", "G", "T", "N")

nucsm <- reshape(nucs, direction="long", varying=1:5, times=names(nucs)[1:5], v.names="count")
names(nucsm)[1] <- "nucleotide"
colnames(nucsm)[colnames(nucsm)=="id"] <- "position"

nucsm$nucleotide <- factor(nucsm$nucleotide, levels = c("A","C","G","T","N"))

barplot <- ggplot(nucsm, aes(x = position, y = count, fill = nucleotide)) +
  geom_bar(position = "fill", width = .75, stat = "identity")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title=element_blank())+
  scale_x_continuous(expand=c(0,0), limit=c(0,nrow(nucs)+1), breaks=seq(0,nrow(nucs)+1,by=scale), sec.axis = dup_axis(name = NULL))+
  scale_y_continuous(expand=c(0,0), sec.axis = dup_axis(name = NULL)) +
  scale_fill_manual("legend", values = c("A" = "gold2",
                                        "C" = "tomato2",
                                        "G" = "cornflowerblue",
                                        "T" = "palegreen2",
                                        "N" = "pink"))+
  xlab(paste("Read Position\n(Total Number of reads = ",sum,")", sep="")) +
  ylab("Proportion")
ggsave(filename= paste(composer_in[1], ".tiff", sep = ""), plot=barplot, width=15, height= 5, dpi=600, compression = "lzw")
invisible(gc())
