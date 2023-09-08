### Hong X. et al 2023
### produce figure 2 
### hongxiaoning@outlook.com
library("ggplot2")
library("ggpubr")
library("ggsci")
library(ggridges)
suppressPackageStartupMessages(library(circlize))


#Figure2a
A_bed = read.table("Aged_brain.bed",header=F)
C_bed = read.table("Young_brain.bed",header=F)
pdf("Figure2a_circos.pdf")
circos.clear()
circos.par("start.degree" = 90)
circos.par("gap.degree" = c(rep(c(2), 20), 22))
circos.initializeWithIdeogram(species="mm10",track.height = 0.3 )
circos.trackHist(C_bed[,1], x = C_bed[,2], bin.size = 1e6, ylim=c(0,1000),
				     col = "#1c72a6", border = "#1c72a6", track.height = 0.15, bg.border = "transparent")
circos.yaxis(side="right")
circos.trackHist(A_bed[,1], x = A_bed[,2], bin.size = 1e6, ylim=c(0,1000),
				     col = "#ed7e00", border = "#ed7e00", track.height = 0.15, bg.border = "transparent")
circos.yaxis(side="right")
text(0, 0.1, "Young: 1,168,079", cex = 1, col = "#1c72a6")
text(0, -0.1, "Aged: 876,918", cex = 1, col = "#ed7e00")
dev.off()

#Figure2b&c
for (group in c("Young", "Aged")) {
  data<-read.table("Brain_eccDNA.group.EPM.stat", sep = "\t", h = T, check.names = F)
  data$Status <-factor(data$Status, levels=c("Young","Aged"))
  data = subset(data, Status==group)
  data$Group <- factor(data$Group, levels = c("OB","mPFC","NAc","CPu","Hip","Cb"))
  pdf(paste0(group,"_EPM_Figure2.pdf"),w=8,h=4)
  ggviolin(data,x="Group",y="EPM", fill="Group", outlier.shape = NA, add = c("boxplot", "jitter"),add.params = list(fill = "white",shape=21, size=1)) +
	  scale_fill_npg()+
      ggtitle(group)+
  theme_classic()+theme(legend.title = element_blank())+theme(legend.position="NA")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold", color="black"))+
  xlab(label = NULL)+ylab(label = paste0("eccDNA Number of Per Million Reads\n") )+
  theme(axis.ticks=element_line(size=0.5,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
  theme(axis.text.x = element_text(face="bold", color="black", size=14, angle=30, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=14),
        legend.text = element_text(face="bold", color="black", size=14),
        axis.title.x = element_text(face="bold", color="black", size=14),
        axis.title.y = element_text(face="bold",color="black", size=14))
  dev.off()
}

#Figure2d&e
for (group in c("Young", "Aged")) {
  data<-read.table(paste0(group, "_GC.txt"), sep = "\t", h = T, check.names = F)
  data$GC=data$GC*100
  Number_Mean = mean(data$GC)
  data$Type=factor(data$Type, levels=rev(c("OB", "mPFC", "NAc", "CPu", "Hip", "Cb")))
  pdf(paste0(group,"_GC_Figure2.pdf"),w=3,h=4)
  ggplot(data, aes(x = GC, y = Type, fill = Type)) +
    geom_density_ridges(scale = 1.5) +
    ggtitle(group)+
    scale_fill_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))) +
    theme_classic()+theme(legend.title = element_blank())+theme(legend.position="NA")+
    theme(panel.grid = element_blank(),axis.title = element_text(size = 12))+
    geom_vline(aes(xintercept = Number_Mean), colour="black", linetype="dashed")+
    xlab(label = "GC content(%)\n")+ylab(label = NULL)+
    theme(plot.title = element_text(hjust = 0.4, face="bold", color="black", size=12))+
    theme(axis.line=element_line(size=0.5,colour="black"))+
    theme(axis.ticks=element_line(size=0.5,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
    theme(axis.text.x = element_text(face="bold", color="black", size=12, angle=90, vjust=0.85, hjust=0.75),
          axis.text.y = element_text(face="bold", color="black", size=12),
          axis.title.x = element_text(face="bold", color="black", size=12),
          axis.title.y = element_text(face="bold",color="black", size=12))
  dev.off()
}

#Figure2f&g
densFindPeak <- function(x){
 td <- density(x)
 maxDens <- which.max(td$y)
 list(x=td$x[maxDens],y=td$y[maxDens])
}

for (group in c("Young", "Aged")) {
  data <- read.table(paste0(group, "_Length.txt"),header=T)
  data$Group <-factor(data$Group, levels =rev(c("OB","mPFC","NAc","CPu","Hip","Cb")))
  for (Type in data$Group) {
    data_Group = subset(data, data$Group == Type)
    limit1 <- data_Group[,3][data_Group[,3]>=0 & data_Group[,3]<=200]
    Peak1=densFindPeak(limit1)
    limit2 <- data_Group[,3][data_Group[,3]>=200 & data_Group[,3]<=500]
    Peak2=densFindPeak(limit2)
    limit3 <- data_Group[,3][data_Group[,3]>=500 & data_Group[,3]<=700]
    Peak3=densFindPeak(limit3)
    limit4 <- data_Group[,3][data_Group[,3]>=700 & data_Group[,3]<=900]
    Peak4=densFindPeak(limit4)
    limit5 <- data_Group[,3][data_Group[,3]>=900 & data_Group[,3]<=1100]
    Peak5=densFindPeak(limit5)
    pdf(paste0(group,"_Density_Figure2.pdf"),h=3,w=4)
    ggplot(data_Group,aes(x=Length,y=..density..))+
	  geom_density(aes(color = Group,fill=Group), alpha=1)+
    annotate("text",x=Peak1[[1]], y= 0.0018,label=paste("",round(Peak1[[1]], 0)),face="bold",color="black",,size=5) +
	  annotate("text",x=Peak2[[1]], y= 0.0045,label=paste("",round(Peak2[[1]], 0)),face="bold",color="black",,size=5) +
	  annotate("text",x=Peak3[[1]], y= 0.0012,label=paste("",round(Peak3[[1]], 0)),face="bold",color="black",,size=5) +
	  annotate("text",x=Peak4[[1]], y= 0.001,label=paste("",round(Peak4[[1]], 0)),face="bold",color="black",,size=5) +
	  annotate("text",x=Peak5[[1]], y= 0.0005,label=paste("",round(Peak5[[1]], 0)),face="bold",color="black",,size=5) +
    xlab("The length distribution of eccDNA")+ylab("Density")+
	  theme_classic() +
	  theme(legend.title = element_blank(), legend.key.size=unit(10, "pt"), legend.text=element_text(face="bold", color="black", size=21)) +
	  theme(legend.position = c(0.8, 0.5)) + #theme(legend.position="NA")+
	  scale_fill_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))[6]) +
	  scale_color_manual(values=rev(c(pal_npg("nrc", alpha = 1)(6)))[6]) +
	  scale_x_continuous(limits = c(0, 2000)) +
	  theme(axis.line=element_line(size=0.5,colour="black"))+
	  theme(axis.ticks=element_line(size=0.5,colour="black"),axis.ticks.length=unit(0.5,"lines"))+
	  theme(axis.text.x = element_text(face="bold", color="black", size=12, angle=90, vjust=0.85, hjust=0.75),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))
    dev.off()
  }
}
