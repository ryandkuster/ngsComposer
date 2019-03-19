library(ggplot2)
library(reshape)


current_folder <- "C:\\Users\\rkuster\\vm_shared_data\\composer_qscores\\Nucs"
basefile <- "se.sample1.R2.fastq.txt"


setwd (current_folder)
input2 <- paste("nucleotides_", basefile, sep = '')
nucs <- read.table(input2, header=F, sep=",")
nuc_headers <- c("A", "C", "G", "T", "N")
names(nucs) <- nuc_headers


#Nucleotide distribution
#nucs <- as.data.frame(t(nucs))
names(nucs) <- gsub("V","", names(nucs))
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
  xlab("Read Position") +
  ylab("Proportion")
ggsave(filename= paste(input2, ".tiff"), plot=barplot, width=15, height= 5, dpi=600, compression = "lzw")
gc()