### Hong X. et al 2023
### produce figure 5 
### hongxiaoning@outlook.com
library("ggplot2")
library("ggpubr")
library("ggsci")
library(lemon)


#Figure5 b
SEM_ebtop<-function(x){
  return(mean(x)+(sd(x)/sqrt(length(x))))
}
SEM_ebbottom<-function(x){
  #return(mean(x)-(sd(x)/sqrt(length(x))))
  return(mean(x))
}
data <- read.table("loci.txt",header=T)
data$Chr<- factor(data$Chr, levels=rev(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY','chrM','UnMapped')))
data = data[data$Type %in% c("OB", "mPFC", "NAc", "CPu", "Hip", "Cb"), ]
data$Type = factor(data$Type, levels=c("OB", "mPFC", "NAc", "CPu", "Hip", "Cb"))
pdf("Figure5b.pdf",w=16,h=8)
ggplot(data,aes(x = Chr,y = Ratio *100))+
	geom_bar(stat = "summary", fun="mean",colour= "black", position=position_dodge(), alpha=0.7) +
	stat_summary(geom='errorbar',fun.max = SEM_ebtop,fun.min=SEM_ebbottom,colour = "black",  width = 0.25, position = position_dodge( .9)) +
	theme_classic() + ylab("eccDNA reads ratio (%)") + xlab("eccDNA reads in different loci") +
	coord_flip() +
	theme(plot.title = element_text(face="bold", color="black", size=16, vjust=0.5, hjust=0.5)) +
	theme(axis.line=element_line(size=0.5,colour="black"))+
	theme(axis.ticks=element_line(size=0.5,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
	theme(axis.text.x = element_text(face="bold", color="black", size=16, angle=90, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=16),
        axis.title.x = element_text(face="bold", color="black", size=24),
		    strip.text = element_text(color="black", size=24),
        axis.title.y = element_text(face="bold",color="black", size=24))+
	theme(strip.background = element_blank(), strip.text = element_text(face="bold", color="black", size=24), strip.placement = "outside") +
	facet_rep_wrap(~ data$Type,ncol=6, scales='free_x', repeat.tick.labels = 'right')
dev.off()

#Figure 5c&d&e
classic_theme <-   theme_classic()+theme(legend.title = element_blank())+theme(legend.position="NA")+
  theme(panel.grid = element_blank(),axis.title = element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.4))+
  #theme(panel.border = element_rect(size=1,colour="black")) +
  theme(axis.line=element_line(size=0.5,colour="black"))+
  theme(axis.ticks=element_line(size=0.5,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
  theme(axis.text.x = element_text(face="bold", color="black", size=12, angle=30, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=12),
        legend.text = element_text(face="bold", color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))
data<-read.table("Merge_Qualimap.txt", sep = "\t", h = T, check.names = T)
data$Status <-factor(data$Group, levels=c("Young","Aged"))
#Figure 5c
Clippedplot<-ggviolin(data,x="Group",y='Clipped_reads...', fill="Group", bxp.errorbar=T,  alpha=0.9, outlier.shape = NA, add = "boxplot", add.params = list(fill = "white")) +
  scale_fill_manual(values=c("#1c72a6","#ed7e00")) +
  theme_classic()+theme(legend.title = element_blank())+theme(legend.position="NA")+
  theme(panel.grid = element_blank(),axis.title = element_text(size=12))+
  xlab(label = NULL)+ylab(label = "Clipped reads(%)\n") +
  classic_theme
ggsave(Clippedplot,file="Figure5c_Clipped_reads.pdf",w=4,h=4)
#Figure 5d
GCplot<-ggboxplot(data,x="Type",y="GC_Percentage...", fill="Group", alpha=0.9, outlier.shape = NA) +
  geom_hline(yintercept = mean(data$GC_Percentage...), linetype=2) +
  scale_fill_manual(values=c("#1c72a6","#ed7e00")) +
  xlab(label = NULL)+ylab(label = "GC(%)\n") +
  classic_theme
ggsave(GCplot,file="Figure5d_GCBoxplot.pdf",w=6,h=3)
#Figure 5e
insertplot<-ggboxplot(data,x="Type",y="Median_insert_size", fill="Group", alpha=0.9, outlier.shape = NA) +
  geom_hline(yintercept = mean(data$Median_insert_size), linetype=2) +
  scale_fill_manual(values=c("#1c72a6","#ed7e00")) +
  xlab(label = NULL)+ylab(label = "Median insert size\n") +
  classic_theme
ggsave(insertplot,file="Figure5e_Median_insert_size.pdf",w=6,h=3)
