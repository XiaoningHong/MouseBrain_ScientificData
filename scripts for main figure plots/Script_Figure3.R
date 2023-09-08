### Hong X. et al 2023
### produce figure 3 
### hongxiaoning@outlook.com
library("ggplot2")
library("ggpubr")
library("ggsci")
library(ggridges)
suppressPackageStartupMessages(library(circlize))


#Figure3 A&b
for (group in c("Young", "Aged")) {
  data<-read.table(paste0(group, "_Genomic_coverage.txt"), sep = "\t", h = T, check.names = F)
  data$Group <-factor(data$Group, levels =rev(c("OB","mPFC","NAc","CPu","Hip","Cb")))
  data$Type <-factor(data$Type, levels=c("3UTR","5UTR","CpG","Exon","Gene2KU","Gene2KD","Intron"))
  group1= subset(data, data$Group=="OB")
  group2= subset(data, data$Group=="mPFC")
  group3= subset(data, data$Group=="NAc")
  group4= subset(data, data$Group=="CPu")
  group5= subset(data, data$Group=="Hip")
  group6= subset(data, data$Group=="Cb")
  data_mean = aggregate(Coverage ~ Type + Group, data, mean)
  pdf(paste0(group,"_Genomic_Figure3.pdf"),w=5,h=5)
  ggplot(data,aes(x=Type,y=Coverage,group=Group))+
  geom_boxplot(data=group1, aes(group=Type), colour="#E64B35FF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group2, aes(group=Type),color="#4DBBD5FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  geom_boxplot(data=group3, aes(group=Type), colour="#00A087FF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group4, aes(group=Type),color="#3C5488FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  geom_boxplot(data=group5, aes(group=Type), colour="#F39B7FFF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group6, aes(group=Type),color="#8491B4FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  geom_point(data=data_mean,aes(group=Type, color = Group, fill= Group), shape=24, size = 2.5) +
  geom_line(data=data_mean,aes(group=Group,color=Group)) +
  scale_fill_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))) +
	scale_color_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))) +
  coord_cartesian(ylim = c(0,3)) +
  geom_hline(yintercept = mean(data$Coverage), linetype=2) +
  ggtitle(group)+
  theme_classic()+
  theme(legend.title = element_blank(),
    legend.position = c(0.95, 0.5), legend.key.size=unit(10, "pt"),
    legend.text=element_text(face="bold", color="black", size=9))+#theme(legend.position="NA")+
  theme(panel.grid = element_blank(),axis.title = element_text(size = 9))+
  theme(plot.title = element_text(hjust = 0.4, face="bold", color="black", size=12))+
  xlab(label = NULL)+ylab(label = paste0("Normalized Genomic Coverage\n") )+
  theme(axis.ticks=element_line(size=1,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
  theme(axis.text.x = element_text(face="bold", color="black", size=9, angle=30, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=9),
        axis.title.y = element_text(face="bold",color="black", size=9)) 
  dev.off()
}

#Figure3 c&d
for (group in c("Young", "Aged")) {
  data<-read.table(paste0(group,"_Ratio.txt"), sep = "\t", h = T, check.names = F)
  data$Group <-factor(data$Group, levels =rev(c("OB","mPFC","NAc","CPu","Hip","Cb")))
  data$Type <-factor(data$Type, levels=rev(c("Transposable Elements","Coding Genes","SINE Repeats","LTR Repeats","LINE Repeats","cis-Regulatory Elements","Enhancer Elements","DNA Repeats","Simple Repeats","Satellite Repeats","Promotor Elements","RNA Repeats","MiroRNAs")))
  group1= subset(data, data$Group=="OB")
  group2= subset(data, data$Group=="mPFC")
  group3= subset(data, data$Group=="NAc")
  group4= subset(data, data$Group=="CPu")
  group5= subset(data, data$Group=="Hip")
  group6= subset(data, data$Group=="Cb")
  pdf(paste0(group,"_Ratio_Figure3.pdf"),w=10,h=5)
  ggplot(data,aes(x=Type,y=Ratio,group=Group))+
  geom_boxplot(data=group1, aes(group=Type), colour="#E64B35FF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group2, aes(group=Type),color="#4DBBD5FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  geom_boxplot(data=group3, aes(group=Type), colour="#00A087FF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group4, aes(group=Type),color="#3C5488FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  geom_boxplot(data=group5, aes(group=Type), colour="#F39B7FFF",alpha=0.5,size=0.5,width=0.2,outlier.shape = NA)+
  geom_boxplot(data=group6, aes(group=Type),color="#8491B4FF",alpha=0.5,size=0.5,width=0.2, outlier.shape = NA)+
  scale_fill_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))) +
	scale_color_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))) +
  coord_flip() +
  ggtitle(group)+
  scale_y_continuous(labels = scales::percent)+
  theme_classic()+
  theme(legend.title = element_blank(),
    legend.key.size=unit(10, "pt"), legend.text=element_text(face="bold", color="black", size=9))+#theme(legend.position="NA")+
  theme(panel.grid = element_blank(),axis.title = element_text(size = 9))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", color="black", size=12))+
  xlab(label = NULL)+ylab(label = paste0("Elements Ratio of eccDNA\n") )+
  theme(axis.ticks=element_line(size=1,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
  theme(axis.text.x = element_text(face="bold", color="black", size=9, angle=30, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=9),
        axis.title.y = element_text(face="bold",color="black", size=9))
  dev.off()
}
