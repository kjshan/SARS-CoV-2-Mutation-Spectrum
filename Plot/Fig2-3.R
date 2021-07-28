#########################################################################

#SARS-CoV-2
#polymorphism 
#Ref:Hu-1
#number:34852

#########################################################################
library(RColorBrewer)
brewer.pal(12, name="Paired")
Self<-c("#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#E31A1C","#FB9A99","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA" )

#SARS-CoV-2 genome
#A      7801
#C      4723
#G      5093
#T      8354


setwd("/Dell/Dell13/shankj/projects/Cov/Plot/Fig2/")
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/SARS_CoV_2_GASID/SCV2_Align.SNP",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

#每个位点的SNP只算1次
SNP<-as.data.frame(Count %>% group_by(SNP,Site)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
sum(plot$count)
#plot<- as.data.frame(SNP %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
#SNP<- as.data.frame(Count %>% group_by(SNP,Sample)  %>% dplyr::summarise(count=n()))


plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("SCV2_polymorphism_Hu-1.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme
dev.off()

fisher.test(matrix(c(877,127,(5093),(4723)),nrow=2))


#########################################################################

#Fig2D
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/A549.txt",header = T,sep="\t")
head(Mutation)
Mutation$SNP<-factor(Mutation$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210716/Fig2D.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=Mutation ,aes(x=log10(SCV2_count),y=log10(SARS_CoV2_Polymorphism),col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+scale_y_continuous(limits = c(1.5,3.5)) +scale_x_continuous(limits = c(0,2)) +
  scale_color_manual(values = Self1)+
  labs(x='SARS-CoV-2 mutation  in Vero',y='SARS-CoV-2 mutation polymorphism') +
  my_theme2+guides(col=F,size=F)
dev.off()
cor.test(log10(Mutation$SCV2_count),log10(Mutation$SARS_CoV2_Polymorphism))


#########################################################################

#MERS-CoV
#polymorphism 
#Ref:NC_019843
#number:258

#########################################################################
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/MERS/MERS_Astr_Align.SNP2",sep="\t",header = F,stringsAsFactors = F)

#Astr
#g	6303	20.90
#a	7937	26.32
#t	9829	32.60
#c	6085	20.18

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

#每个位点的SNP只算1???
SNP<-as.data.frame(Count %>% group_by(SNP,Site)  %>% dplyr::summarise(count=n()))
head(SNP)
#plot<- as.data.frame(SNP %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))

plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("MERS_polymorphism.pdf",width = 2,height = 2)
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/MERS_polymorphism_Astr.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  scale_y_continuous(breaks = c(seq(0,300,50)),limits = c(0,300)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()

fisher.test(matrix(c(83,	25,(6303-83),(6085-25)),nrow=2))
#########################################################################

#four-fold degenerate site


Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/Codon/SCV2_NC_045512.SNP",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("CDS","Site","SNP","Sample")
head(Count)

#每个位点的SNP只算1次
SNP<-as.data.frame(Count %>% group_by(SNP,Site,CDS)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
#SNP<- as.data.frame(Count %>% group_by(SNP,Sample)  %>% dplyr::summarise(count=n()))


plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("SCV2_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
dev.off()

# Ref Count
#   A  1774
#   C   830
#   G   395
#   T  3195
fisher.test(matrix(c(418,606,(3*395-418),(3*830-606)),nrow=2))

#########################################################################
setwd("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/")
#tmp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Genome/",pattern = "*.txt")
tmp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Astr/",pattern = "*.SNP2")

Filter<-data.frame()
for (i in 1:length(tmp)) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Astr/",tmp[i]),header = F,sep="\t",stringsAsFactors = F)
  #Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Genome/",tmp[i]),header = F,sep="\t",stringsAsFactors = F)
  name=sub(".SNP.txt","",tmp[i])
  Mutation$sample<-rep(name,nrow(Mutation))
  Filter<-rbind(Filter,Mutation)
  rm(Mutation)
}
head(Filter)
colnames(Filter)<-c("Site","SNP","Variant","Segment")

#每个位点的SNP只算1次
SNP<-as.data.frame(Filter %>% group_by(SNP,Site,Segment)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))

plot$SNP <- c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C")
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#Ref Count
#C	3057	24.03
#T	4270	33.57
#A	2938	23.10
#G	2456	19.31

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210625/H1N1_polymorphism_Astr.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(breaks = c(seq(0,900,300)),limits = c(0,900)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)
dev.off()
fisher.test(matrix(c(156,	85,(2456-156),(3057-85)),nrow=2))

fisher.test(matrix(c(855,629,(3057-855),(2456-629)),nrow=2))

#########################################################################

#Ebola virus mutations in 293T cells

#########################################################################

Ebola_293T <- as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ebola/293T.txt",stringsAsFactors = F))
Ebola_293T$Type <- rep(c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C"),7)
Ebola_293T$Type<-factor(Ebola_293T$Type,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
head(Ebola_293T)
Ebola_293T_plot<-summarySE(Ebola_293T, measurevar="MLE", groupvars=c("Type"))

pdf("Ebola_293T.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=Ebola_293T, aes(x=Type,y=MLE,fill=Type)) + 
  geom_boxplot()+
  scale_fill_manual(values = Self)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='Mutation frequency',title = "Ebola (-) from 293T" )+ 
  #facet_grid(~Sample)+
  #scale_y_continuous(limits = c(0,300)) +
  my_theme
dev.off()