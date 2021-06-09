#########################################################################

#fastML
#Bat and human and Pangolin polymorphism

#########################################################################
#Human Bat Astr
#g	5812	19.41
#a	9011	30.10
#t	9680	32.33
#c	5436	18.16

setwd("/Dell/Dell13/shankj/projects/Cov/Plot/Fig6/fastML/")
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig6/fastML/Human_Bat_SNP.txt",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

plot<-as.data.frame(Count %>% filter(Sample %in% c("RATG13")) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))


plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
plot

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/RaTG13_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(breaks = c(seq(0,250,50)),limits = c(0,250)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()


fisher.test(matrix(c(8,11,(5812-8),(5436-11)),nrow=2))

#########################################################################
#Human

plot<-as.data.frame(Count %>% filter(Sample %in% c("WUHAN-HU-1")) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))

plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/Hu-1_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()

fisher.test(matrix(c(7,9,(5812-7),(5436-9)),nrow=2))

#########################################################################

#Human Bat pangolin Astr
#g	5761	19.24
#a	9079	30.32
#t	9716	32.45
#c	5383	17.98


Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig6/fastML/Human_Bat_GD_GX_SNP.txt",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)
unique(Count$Sample)
#GD Pangolin

plot<-as.data.frame(Count %>% filter(Sample %in% c("GUANGDONG")) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
head(SNP)

plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/GD_Pangolin_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(breaks = c(seq(0,400,100)),limits = c(0,400)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()


fisher.test(matrix(c(37,53,(5761-37),(5383-53)),nrow=2))

#GX Pangolin

plot<-as.data.frame(Count %>% filter(Sample %in% c("GUANGXI_P5L")) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))


plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/GX_Pangolin_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(breaks = c(seq(0,600,100)),limits = c(0,600)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()


fisher.test(matrix(c(103,121,(5761-103),(5383-121)),nrow=2))
#########################################################################
#Human and Bat Astr

plot<-as.data.frame(Count %>% filter(Sample %in% c("N2")) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
head(SNP)

plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/Human_Bat_Astr_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()

fisher.test(matrix(c(13,10,(5761-13),(5383-10)),nrow=2))
############################################################################
#CirSeq

Ebola_EpoNi <- read.csv("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ebola/Ebola_EpoNi.csv",header = T)
Ebola_EpoNi_2<-summarySE(Ebola_EpoNi, measurevar="Norm", groupvars=c("Sample"))
Ebola_EpoNi_2$Recourse<-rep("Bat",7)
Ebola_293T <- read.csv("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ebola/Ebola_293T.csv",header = T)
Ebola_293T_2<-summarySE(Ebola_293T, measurevar="Norm", groupvars=c("Sample"))
Ebola_293T_2$Recourse<-rep("Human",7)


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/Ebola_CirSeq.pdf",width = 3,height = 3)

ggplot(Ebola, aes(x= as.factor(Sample), y = Norm,fill=Recourse)) +
  geom_bar(position="dodge", stat="identity", colour="black",width = 0.7)+
  geom_errorbar(aes(ymin = (`Norm` - se), ymax = (`Norm` + se)), 
                position = position_dodge(0.7),width=0.4,col='black',size=.5 )+  
  geom_hline(yintercept=1,col="black",linetype="dashed")+
  labs(x='Passage',y='The ratio between C > A and G > U\n mutation rate normalized by passage1')+
  my_theme2+guides(fill=F)
dev.off()
