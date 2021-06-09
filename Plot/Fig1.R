###########################################################

#SARS-CoV-2 mutation spectrum

###########################################################

#color and themes
library(RColorBrewer)
brewer.pal(12, name="Paired")
Self1<-c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99" )
Self<-c("#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#E31A1C","#FB9A99","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA","#AAAAAA" )
Self3<-c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928" )

theme_set(theme_bw()+theme(panel.grid=element_blank(),panel.border=element_rect(size=1,color="black")))
my_theme<-theme(axis.line.x=element_line(size=0,color="black"),axis.line.y=element_line(size=0,color="black"),
                axis.ticks=element_line(size=0.5,color="black"),axis.ticks.length=unit(0.05,"inches"),
                axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
                axis.text.x = element_text(angle = 45,hjust = 1,size=8,color="black"),
                axis.text.y =  element_text(size=10,color="black"),
                strip.text.x = element_text(size=10,face = "bold"),
                strip.background = element_rect(color = "black",size=1),
                legend.position = 1,
                legend.text = element_text(size=10),legend.title = element_text(size=10))


###########################################################                

setwd("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/")

###########################################################

#SARS-CoV-2 mutation the in negative strand

###########################################################
#1.parameter's cut-off
Ctrl<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Vero_Total_Barcode_SNP.csv",stringsAsFactors = F))
Gaussi<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Guassi.txt",sep = "\t",header = F)
head(Ctrl)
head(Gaussi)
colnames(Gaussi)<-c("Pos","Mismatch_Number","Total_Number")
Gaussi$Frac<-Gaussi$Mismatch_Number/Gaussi$Total_Number
Ctrl<-merge(Ctrl,Gaussi,by=c("Pos"))
write.csv(Ctrl,"/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Vero_Total_Barcode_SNP_Mismatch.csv",row.names = F,quote = F)

Ctrl<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Vero_Total_Barcode_SNP_Mismatch.csv",stringsAsFactors = F))

#FigS2A
Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),#The four konwn polymorphisms sites, which were different between our SARS-CoV-2 reference genome and BetaCoV/Korea/KCDC03/2020
               #Dis >= 15,
               #UMI_Alt_no_PCR_reads>=2,
               #UMI_ref_reads==0,
               #UMI_reads<=20,
               #All_Alt_reads/Total_Number <=  0.002,
               !UMI %in% c("26257-26283")) #Known indel in BetaCoV/Korea/KCDC03/2020 
head(Filter)

plot<-as.data.frame(Filter) %>% group_by(SNP)  %>% dplyr::summarise(count=n())
max(plot$count)

pdf("Control_noCutoff_1.pdf",width = 3,height = 3)
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+
  labs(x='',y='Mutation count')+ 
  scale_y_continuous(breaks = seq(0,24000,6000),limits = c(0,24000))+
  my_theme
dev.off()

#FigS2B
Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),
               #Dis >= 15,
               UMI_Alt_no_PCR_reads>=2,
               UMI_ref_reads==0,
               #UMI_reads<=20,
               # All_Alt_reads/Total_Number <=  0.002,
               !UMI %in% c("26257-26283")) #,
#UMI_Alt_no_PCR_reads==1 
head(Filter)

plot<-as.data.frame(Filter) %>% group_by(SNP)  %>% dplyr::summarise(count=n())
max(plot$count)
pdf("Control_Cutoff_C1and2.pdf",width = 3,height = 3)
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+
  labs(x='',y='Mutation count')+ my_theme
#my_theme
dev.off()

#FigS2C
Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),
               Dis >= 15,
               UMI_Alt_no_PCR_reads>=2,
               UMI_ref_reads==0,
               #UMI_reads<=20,
               # All_Alt_reads/Total_Number <=  0.002,
               !UMI %in% c("26257-26283")) #,
#UMI_Alt_no_PCR_reads==1 
head(Filter)

plot<-as.data.frame(Filter) %>% group_by(SNP)  %>% dplyr::summarise(count=n())
max(plot$count)
pdf("Control_Cutoff_C1and2and3.pdf",width = 3,height = 3)
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+scale_y_continuous(limits = c(0,250))+
  labs(x='',y='Mutation count')+ my_theme
#my_theme
dev.off()

###########################################################

#To discard potential polymorphisms
#cut-off 0.2%

###########################################################

Gaussi<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Guassi.txt",sep = "\t",header = F)

head(Gaussi)
colnames(Gaussi)<-c("Pos","Mismatch_Number","Total_Number")
Gaussi$Frac<-Gaussi$Mismatch_Number/Gaussi$Total_Number

#FigS2D
ggplot(Gaussi) +
  #geom_vline(xintercept = log2(20),color="red",linetype="dotted")+
  labs(x='log10(mismatch read number/total read number)')+ 
  #scale_x_continuous(limits = c(0,0.01)) +
  geom_density(aes(x =log10(Frac), y =..count..))

Gaussi<-filter(Gaussi, Frac<0.9,Frac>0)

library(mixtools)
mid<-mixtools::normalmixEM(log10(Gaussi$Frac), arbvar = T, epsilon = 1e-03)

mid$lambda
mid$mu
mid$sigma

pnorm(mean = -2.356281,sd=0.5475604,lower.tail = T,q=log10(0.002))*0.06472484
pnorm(mean = -2.966511,sd=0.2419629,lower.tail = T,q=log10(0.002))*0.93527516

0.01719789/(0.01719789+0.80955)

ggplot(Gaussi) +
  #geom_vline(xintercept = log2(20),color="red",linetype="dotted")+
  labs(x='log10(mismatch read number/total read number)')+ 
  #scale_x_continuous(limits = c(0,0.01)) +
  geom_density(aes(x =log10(Frac)))+my_theme

pdf("Gaussi.pdf",width = 3,height = 3)
data.frame(x = mid$x) %>%
  ggplot() +
  geom_histogram(aes(x, y =..density..), binwidth = .2, colour = "grey30", fill = "grey", lwd = .5) +
  stat_function(geom = "line", fun = plot_mix_comps, # here is the function
                args = list(mid$mu[1], mid$sigma[1], lam = mid$lambda[1]),
                colour = "red", lwd = 1) +
  stat_function(geom = "line", fun = plot_mix_comps, # and here again because k = 2
                args = list(mid$mu[2], mid$sigma[2], lam = mid$lambda[2]),
                colour = "blue", lwd = 1) +
  labs(x='log10(mismatch read number/total read number)',y="Density")+ my_theme+  
  geom_vline(xintercept=log10(0.002),col="black",linetype="dashed",lwd=1)
dev.off()                

#Fig1F
Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),
               Dis >= 15,
               UMI_Alt_no_PCR_reads>=2,
               UMI_ref_reads==0,
               #UMI_reads<=20,
               All_Alt_reads/Total_Number <=  0.002,
               !UMI %in% c("26257-26283")) #,
#UMI_Alt_no_PCR_reads==1 
head(Filter)

#mutation spectrum
plot<-as.data.frame(Filter) %>% group_by(SNP)  %>% dplyr::summarise(count=n())
sum(plot$count)
pdf("Fig1F.pdf",width = 3,height = 3)
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+
  labs(x='',y='Mutation count')+ my_theme
#my_theme
dev.off()

#Count
ts<-sum(51,15,35,21)
tv<-sum(8,4,6,6,9,24,11,7)
ts+tv
binom.test(ts,(ts+tv),4/12)
binom.test(24,tv,1/8)

#Freq
#Coverage_consensus_reads.txt
#Base count in Consensus reads
#A	3631127
#T	3212571
#C	2354397
#G	2175628
wilcox.test(c(9.65E-06,2.17E-05,6.90E-06,6.55E-06),
            c(2.21E-06,1.10E-06,2.55E-06,2.55E-06,4.14E-06,1.10E-05,3.43E-06,2.18E-06))

t.test(c(2.21E-06,1.10E-06,2.55E-06,2.55E-06,4.14E-06,1.10E-05,3.43E-06,2.18E-06), mu = 1.10E-05)

#Fisher exact test
fisher.test(matrix(c(22,6,(5863-22),(5492-6)),nrow=2))
#Potential_mutation_site.txt 
#  Var1 Freq
#    A 8929
#    C 5492
#    G 5863
#    T 9594

###########################################################

#Consenus VS Unconsensus
#Fig1D

###########################################################

#mismatch frequency
#Consensuse
Consensus<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/ConsensusVSunconsensus/Consensus.result",stringsAsFactors = F,header = F))
colnames(Consensus)<-c("Aver","SNP_Freq")
head(Consensus)

##Inconsensuse
Inconsensuse<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/ConsensusVSunconsensus/Unconsensus.result",stringsAsFactors = F,header = F))

colnames(Inconsensuse)<-c("BQ","SNP_Freq")
head(Inconsensuse)

pdf("ConsensusVSInconsensue_Mismatch_Freq.pdf",height = 3,width = 3,useDingbats = F)
ggplot()+
  geom_point(data=Consensus,aes(x=Aver,y=log10(SNP_Freq)),col="black")+
  geom_point(data=Inconsensuse,aes(x=BQ,y=log10(SNP_Freq)),col="red")+
  scale_y_continuous(breaks = seq(-5,0,1),limits = c(-5,0))+
  my_theme
dev.off()

###################################################################

#Small indel

###################################################################

test<-read.table("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/Polymorphism_Consensuse.Indel.txt",header = F,sep="\t")
colnames(test)<-c("Fraction","Barcode","Barcode_All_Read_pair_number","Covered_reads_number","Covered_readpair_number","Pos","Indel","Alt","Indel_len","Dis","Reads_Pos","Reads","Reads_PCR","non_PCR_number","SNP_read_number")


Filter<-filter(test,Fraction<=0.002,abs(Dis)>=15)
write.csv(Filter,"Indel.csv",row.names = F,quote = F)

Count<-Filter %>% group_by(Pos,Indel,Alt,Indel_len) %>% dplyr::summarise(count=n())
nrow(Count)
ggplot()+
  geom_histogram(data=filter(Count,count==1),aes(x=Indel_len),binwidth=1)+ 
  #scale_x_continuous(limits = c(0,220)) +
  labs(x='Indel length(Insertion:+ \t Deletion:-)')+
  #my_theme2+guides(fill=F)
  my_theme1

#Distance 
ggplot()+
  geom_histogram(data=Filter,aes(x=Dis),binwidth=1)+ 
  scale_x_continuous(limits = c(-260,170)) +
  labs(x='Distance to junction site')+#my_theme1
  my_theme2+guides(fill=F)

#read distance
as.vector(test$Reads_Pos)
mid<-as.data.frame(na.omit(as.numeric(unlist(map(test$Reads_Pos,~str_split(.,','))))))
colnames(mid)<-"ReadPosition"
ggplot()+
  geom_histogram(data=mid,aes(x=ReadPosition),binwidth=1)+ 
  #scale_x_continuous(limits = c(0,220)) +
  labs(x='Read Position')+my_theme
