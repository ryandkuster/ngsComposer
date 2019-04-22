options(warn=-1)
library(ggplot2)
library(reshape)


x <- commandArgs(trailingOnly = TRUE)


qscores <- read.table(x[2], header=F, sep=",")
names(qscores) <- gsub("V","", names(qscores))
sum <- format(round(as.numeric(sum(qscores[1,])), 1), nsmall=0, big.mark=",")
xyz <- stack(qscores)
for (i in levels(xyz$ind)){
j = strtoi(i) - 1
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
ggsave(filename= paste(x[2], ".tiff", sep = ""), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw")
invisible(gc())


nucs <- read.table(x[1], header=F, sep=",")
names(nucs) <- c("A", "C", "G", "T", "N")
# names(nucs) <- gsub("V","", names(nucs))
nucsm <- melt(cbind(nucs, ind = rownames(nucs)), id.vars = c('ind'))
names(nucsm)[2] <- "Nucleotide"
nucsm$ind=as.numeric(levels(nucsm$ind))[nucsm$ind]
barplot <- ggplot(nucsm, aes(x = ind, y = value, fill = Nucleotide)) + 
geom_bar(position = "fill", width = .75, stat = "identity")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title=element_blank())+
scale_x_continuous(expand=c(0,0), limit=c(0,nrow(nucs)), breaks=seq(0,nrow(nucs),by=5))+
scale_fill_manual("legend", values = c("A" = "gold2",
                                       "C" = "tomato2",
                                       "G" = "cornflowerblue",
                                       "T" = "palegreen2",
                                       "N" = "pink"))+
xlab(paste("Read Position (Total Number of reads = ",sum,")", sep="")) +
ylab("Proportion")
ggsave(filename= paste(x[1], ".tiff", sep = ""), plot=barplot, width=15, height= 5, dpi=600, compression = "lzw")
invisible(gc())
