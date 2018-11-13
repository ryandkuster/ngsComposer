library(plyr)
library(car)
library(MASS)
library(data.table)
library(ggplot2)
library(reshape)

setwd ("/home/rkuster/Desktop/QC_plotting")

xyz <- read.table("quals.txt", header=F, sep=",")
names(xyz) <- gsub("V","", names(xyz))
boxplot <- ggplot(stack(xyz), aes(x = ind, y=values), stat='identity')+
  geom_boxplot(fill="white", colour="cornflowerblue",
               outlier.alpha = 0.25, outlier.colour="tomato", outlier.size=2.0)+
  coord_flip()+
  scale_x_discrete((xyz$ind)+1, limits=c(1:40))+
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x=element_text(colour="cornflowerblue", size=18),
        axis.text.y=element_text(colour="cornflowerblue", size=18),
        axis.title=element_text(size=24)) +
  xlab("Quality Scores") +
  ylab("Frequency")
ggsave(filename="Frequency_QC_Score.tiff", plot=boxplot, width=15, height= 15, dpi=600, compression = "lzw")
gc

#Now transform data
qc <- reshape(xyz, direction="long", varying=list(c(1:ncol(xyz))), v.names="count")
qc$score <- qc$time-1
qc <- rename(qc, c("id"="position"))
qc <- subset(qc, select=c("position","score","count"))
boxplot <- ggplot(qc[which(qc$count>0),], aes(x = position, y=score, weight=count, group=position), stat='identity')+
  geom_boxplot(fill="white", colour="cornflowerblue",
               outlier.alpha = 0.2, outlier.colour="tomato", outlier.size=2)+
   scale_x_continuous(expand=c(0,0), limit=c(0,nrow(xyz)), breaks=seq(0,nrow(xyz),by=1)) +
  scale_y_continuous(expand=c(0,0), limit=c(0,40), breaks=seq(0,40, by=2)) +
  theme(axis.text.x=element_text(colour="cornflowerblue", size=8),
        axis.text.y=element_text(colour="cornflowerblue", size=8),
        axis.title=element_text(size=12)) +
  xlab("Read Position") +
  ylab("Quality Scores")
ggsave(filename="QC_plot.tiff", plot=boxplot, width=15, height= 5, dpi=600, compression = "lzw")
gc
