options(warn=-1)

suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
suppressMessages(library(Hmisc, warn.conflicts = FALSE))

composer_in <- commandArgs(trailingOnly = TRUE)

qscores <- read.table(composer_in[2], header=F, sep=",")
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
qc <- qc[which(qc$count>0),]
qc$mean <- "NA"
qc$sd <- "NA"

#function across all positions to calculate mean value for each position
#source code for weighted_se: https://rdrr.io/cran/diagis/src/R/ses.R
for (i in c(1:nrow(qscores))) {
  wt <- (qscores[i,])/sum(qscores[i,])
  x <- c(0:40)
  wm <- wtd.mean(x, wt, normwt=FALSE, na.rm=TRUE)
  wsd <- sqrt(wtd.var(x, wt, normwt=FALSE, na.rm=TRUE, method=c('ML')))
  qc$mean[qc$position==i] <- wm
  qc$sd[qc$position==i] <- wsd
}
qc[,4] <- as.numeric(as.character(qc[,4]))
qc[,5] <- as.numeric(as.character(qc[,5]))
qc$minSD <- round((qc$mean - qc$sd), digits=0)
qc <- as.data.frame(qc)
qc <- qc[order(qc$position),]


# create qscore distribution visualization (no outliers)
boxplot <- ggplot(qc, aes(x = position, y=score, weight=count, group=position), stat='identity')+
  geom_boxplot(fill="white", colour="cornflowerblue", lwd = 0.35, outlier.shape = NA)+
  geom_hline(aes(yintercept = 30), lty=1, alpha=0.05, size=0.35)+
  geom_hline(aes(yintercept = 20), lty=1, alpha=0.05, size=0.35)+
  stat_summary(aes(x = position, y = mean, group = position),
               shape=23, size=0.2, color="black", alpha=0.25)+
  scale_x_continuous(expand=c(0,0), limit=c(0,nrow(qscores)+1), breaks=seq(0,nrow(qscores)+1,by=5),
                     sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(expand=c(0,0), limit=c(0,43), breaks=seq(0,43, by=2), sec.axis = dup_axis(name = NULL)) +
  theme_classic()+
  xlab(paste("Read Position (Total Number of reads = ",sum,")", sep="")) +
  ylab("Quality Scores (phred+33)")
suppressMessages(ggsave(filename= paste(composer_in[2], ".tiff", sep = ""), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw"))
invisible(gc())


# create nucleotide distribution visualization (with outliers and sd)
boxplot <- ggplot(qc, aes(x = position, y=score, weight=count, group=position), stat='identity')+
  geom_boxplot(fill="white", colour="cornflowerblue", lwd = 0.35,
               outlier.alpha = 0.25, alpha=0.25, outlier.colour="tomato", outlier.size=1.25)+
  geom_hline(aes(yintercept = 30), lty=1, alpha=0.05, size=0.35)+
  geom_hline(aes(yintercept = 20), lty=1, alpha=0.05, size=0.35)+
  stat_summary(aes(x = position, y=mean, group=position),
               shape=23, size=0.2, color="black", alpha=0.25)+
  geom_linerange(aes(x = position, ymin=round(mean-sd,0), ymax=round(mean+sd,0), group=position),
                 color="grey", alpha=0.25, size=0.5)+
  scale_x_continuous(expand=c(0,0), limit=c(0,nrow(qscores)+1), breaks=seq(0,nrow(qscores)+1,by=5),
                     sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(expand=c(0,0), limit=c(0,43), breaks=seq(0,43, by=2), sec.axis = dup_axis(name = NULL)) +
  theme_classic()+
  xlab(paste("Read Position (Total Number of reads = ",sum,")", sep="")) +
  ylab("Quality Scores (phred+33)")
suppressMessages(ggsave(filename= paste(composer_in[2], "_outliers.tiff", sep = ""), plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw"))
invisible(gc())


# create nucleotide distribution visualization
nucs <- read.table(composer_in[1], header=F, sep=",")
names(nucs) <- c("A", "C", "G", "T", "N")
nucsm <- melt(cbind(nucs, ind = rownames(nucs)), id.vars = c('ind'))
names(nucsm)[2] <- "Nucleotide"
nucsm$ind=as.numeric(levels(nucsm$ind))[nucsm$ind]
barplot <- ggplot(nucsm, aes(x = ind, y = value, fill = Nucleotide)) + 
geom_bar(position = "fill", width = .75, stat = "identity")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"),
      legend.title=element_blank())+
scale_x_continuous(expand=c(0,0), limit=c(0,nrow(nucs)+1), breaks=seq(0,nrow(nucs)+1,by=5))+
scale_fill_manual("legend", values = c("A" = "gold2",
                                       "C" = "tomato2",
                                       "G" = "cornflowerblue",
                                       "T" = "palegreen2",
                                       "N" = "pink"))+
xlab(paste("Read Position (Total Number of reads = ",sum,")", sep="")) +
ylab("Proportion")
ggsave(filename= paste(composer_in[1], ".tiff", sep = ""), plot=barplot, width=15, height= 5, dpi=600, compression = "lzw")
invisible(gc())
