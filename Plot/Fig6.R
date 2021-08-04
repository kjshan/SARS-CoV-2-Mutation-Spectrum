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

#G>U vs C>A
fisher.test(matrix(c(8,11,(5812-8),(5436-11)),nrow=2))
#C>U vs G>A
fisher.test(matrix(c(134,54,(5436-134),(5812-54)),nrow=2))

#G>U vs C>A
fisher.test(matrix(c(8,11,(5812),(5436)),nrow=2))
#C>U vs G>A
fisher.test(matrix(c(134,54,(5436),(5812)),nrow=2))
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
#G>U vs C>A
fisher.test(matrix(c(7,9,(5812-7),(5436-9)),nrow=2))
#C>U vs G>A
fisher.test(matrix(c(125,33,(5436-125),(5812-33)),nrow=2))

#G>U vs C>A
fisher.test(matrix(c(7,9,(5812),(5436)),nrow=2))
#C>U vs G>A
fisher.test(matrix(c(125,33,(5436),(5812)),nrow=2))

#########################################################################

#Correlation between RATG13 and SARS-CoV-2
Mid<-Count %>% group_by(Sample,SNP)  %>% dplyr::summarise(count=n())
Mid<-spread(Mid, Sample, count)


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/HuVSRatG13_2.pdf",width = 3,height = 3,useDingbats = F)
Mid$SNP<-factor(Mid$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Mid, aes(x=`WUHAN-HU-1`,y=RATG13,col=SNP)) + 
  geom_point(size=3)+
  scale_color_manual(values = Self1)+
  scale_y_continuous(breaks = c(seq(0,250,50)),limits = c(0,250)) +
  labs(x='Hu-1',y='RatG13')+my_theme2+guides(col=F)

dev.off()

cor.test(log10(Mid$RATG13),log10(Mid$`WUHAN-HU-1`))

SCV2_Poly<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/VeroSNP_Polymorphism.txt",header = T,sep="\t")
Poly_MCAR<-merge(Mid,SCV2_Poly,by="SNP")
head(Poly_MCAR)

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/HuVPoly_2.pdf",width = 3,height = 3,useDingbats = F)
Poly_MCAR$SNP<-factor(Poly_MCAR$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Poly_MCAR, aes(x=`WUHAN-HU-1`,y=Polymorphism,col=SNP)) + 
  geom_point(size=3)+
  scale_color_manual(values = Self1)+
  scale_y_continuous(breaks = c(seq(0,2500,500)),limits = c(0,2700)) +
  labs(x='Hu-1',y='Polymorphism')+my_theme2+guides(col=F)

dev.off()

cor.test(log10(Poly_MCAR$Polymorphism),log10(Mid$`WUHAN-HU-1`))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig6/RaTG13VPoly_1.pdf",width = 3,height = 3,useDingbats = F)
Poly_MCAR$SNP<-factor(Poly_MCAR$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Poly_MCAR, aes(x=RATG13,y=Polymorphism,col=SNP)) + 
  geom_point(size=3)+
  scale_color_manual(values = Self1)+#scale_x_continuous(breaks = c(seq(0,250,50)),limits = c(0,250)) +
  #scale_y_continuous(breaks = c(seq(0,2500,500)),limits = c(0,2700)) +
  labs(x='Polymorphism in Human',y='MARC to RATG13')#+my_theme2+guides(col=F)

dev.off()

cor.test(log10(Poly_MCAR$Polymorphism),log10(Mid$RATG13))

#########################################################################