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

#Ctrl<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Vero_Total_Barcode_SNP_Mismatch.csv",stringsAsFactors = F))

#FigS2A
Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),#The four konwn polymorphisms sites, which were different between our SARS-CoV-2 reference genome and BetaCoV/Korea/KCDC03/2020
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
               UMI_Alt_no_PCR_reads>=2,
               UMI_ref_reads==0,
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
#FigS2D
#cut-off 0.2%

Gaussi<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/Guassi.txt",sep = "\t",header = F)

head(Gaussi)
colnames(Gaussi)<-c("Pos","Mismatch_Number","Total_Number")
Gaussi$Frac<-Gaussi$Mismatch_Number/Gaussi$Total_Number
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


pnorm(mean = -2.423157,sd=0.5649736,lower.tail = T,q=log10(0.002))*0.06368319
pnorm(mean = -3.044717 ,sd=0.2420202,lower.tail = T,q=log10(0.002))*0.93631681

0.01991428/(0.01991428+0.8646311)

ggplot(Gaussi) +
  labs(x='log10(mismatch read number/total read number)')+ 
  geom_density(aes(x =log10(Frac)))+my_theme

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210609/FigS3_E-2.pdf",width = 3,height = 3)
data.frame(x = mid$x) %>%
  ggplot() +
  geom_histogram(aes(x, y =..density..), binwidth = .2, colour = "grey30", fill = "grey", lwd = .5) +
  stat_function(geom = "line", fun = plot_mix_comps, # here is the function
                args = list(mid$mu[1], mid$sigma[1], lam = mid$lambda[1]),
                colour = "red", lwd = 1) +
  stat_function(geom = "line", fun = plot_mix_comps, # and here again because k = 2
                args = list(mid$mu[2], mid$sigma[2], lam = mid$lambda[2]),
                colour = "blue", lwd = 1) +
  labs(x='log10(mismatch read number/total read number)',y="Density")+ my_theme+#2+guides(fill=F)+
  geom_vline(xintercept=log10(0.002),col="black",linetype="dashed",lwd=1)
dev.off()                
head(Ctrl)


###########################################################

#De novo mutations in Vero
#Fig1F

###########################################################

Filter<-filter(Ctrl, 
               !Pos %in% c(4402,5062,8782,28144),
               Dis >= 15,
               UMI_Alt_no_PCR_reads>=2,
               UMI_ref_reads==0,
               All_Alt_reads/Total_Number <=  0.002,
               !UMI %in% c("26257-26283")) #,
#UMI_Alt_no_PCR_reads==1 
head(Filter)
CT_Pos<-filter(Filter,SNP=="C > T")$Pos
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
ts<-c(9.65E-06,2.17E-05,6.90E-06,6.55E-06)
tv<-c(2.21E-06,1.10E-06,2.55E-06,2.55E-06,4.14E-06,1.10E-05,3.43E-06,2.18E-06)
t.test(log10(ts),log10(tv))

wilcox.test(c(9.65E-06,2.17E-05,6.90E-06,6.55E-06),
       c(2.21E-06,1.10E-06,2.55E-06,2.55E-06,4.14E-06,1.10E-05,3.43E-06,2.18E-06))

t.test(c(2.21E-06,1.10E-06,2.55E-06,2.55E-06,4.14E-06,1.10E-05,3.43E-06,2.18E-06), mu = 1.10E-05)

#C>U VS G>A
fisher.test(matrix(c(51,15,2354397,2185628),nrow = 2))
# G>U VS C>A
fisher.test(matrix(c(24,6,2185628,2354397),nrow = 2))


#Fisher exact test, only use the sites covered by junction read pairs
fisher.test(matrix(c(22,6,(5863-22),(5492-6)),nrow=2))
#Potential_mutation_site.txt 
#  Var1 Freq
#    A 8929
#    C 5492
#    G 5863
#    T 9594

###########################################################

#Consenus VS Unconsensus
#FigS3A

###########################################################

#mismatch frequency
#Consensuse
Consensus<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/ConsensusVSunconsensus/Consensus.result2",stringsAsFactors = F,header = F))
colnames(Consensus)<-c("Aver","SNP","Cover","SNP_Freq")
head(Consensus)

##Inconsensuse
Inconsensuse<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/ConsensusVSunconsensus/Unconsensus.result2",stringsAsFactors = F,header = F))

colnames(Inconsensuse)<-c("BQ","SNP","Cover","SNP_Freq")
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
#FigS3D-E

###################################################################

test<-read.table("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/Polymorphism_Consensuse.Indel.txt",header = F,sep="\t")
colnames(test)<-c("Fraction","Barcode","Barcode_All_Read_pair_number","Covered_reads_number","Covered_readpair_number","Pos","Indel","Alt","Indel_len","Dis","Reads_Pos","Reads","Reads_PCR","non_PCR_number","SNP_read_number")


Filter<-filter(test,Fraction<=0.002,abs(Dis)>=15)
write.csv(Count,"/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Vero_JCR_default/JCR_indel/Indel.csv",row.names = F,quote = F)

Count<-Filter %>% group_by(Pos,Indel,Alt,Indel_len) %>% dplyr::summarise(count=n())
nrow(Count)
ggplot()+
  geom_histogram(data=filter(Count,count==1),aes(x=Indel_len),binwidth=1)+ 
  labs(x='Indel length(Insertion:+ \t Deletion:-)')+
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

###################################################################

#A549

###################################################################

tmp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Read2/",pattern = "*.UMI.mutation.countPvalue.csv")
tmp2 <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Read2/",pattern = "*.mismatch")
Filter<-data.frame()
for (i in 1:3) {
  Mismatch<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Read2/",tmp2[i]),header = T,stringsAsFactors = F)
  colnames(Mismatch)<-c("Pos","Mismatch_Number","Total_Number")
  Mutation<-read.csv(paste0("/Dell/Dell13/shankj/projects/Cov/SARS_CoV2/Read2/",tmp[i]),header = T,stringsAsFactors = F)
  name=sub(".UMI.mutation.countPvalue.csv","",tmp[i])
  #Mutation$Qvalue<-qvalue(Mutation$Pvalue)$qvalue
  Mid<-merge(Mutation,Mismatch,by=c("Pos"))
  Mid$sample<-rep(name,nrow(Mid))
  Filter<-rbind(Filter,Mid)
  rm(Mid)
  rm(Mutation)
  
}
#Filter$Qvalue<-qvalue(Filter$Pvalue)$qvalue
head(Filter)
nrow(Filter)


write.csv(Filter,"A549.csv",row.names = F,quote = F)
head(Mid)
#C8782T,T28144C,C18060T
Mid<-filter(Filter,Dis >= 15,
            UMI_Alt_no_PCR_reads>=2,
            UMI_ref_reads==0,
            !Pos %in% c(18060,8782,28144)) #Dis >= 15, Count <=5,UMI_Alt_Fraction==1
write.csv(Mid,"A549_Filter",row.names = F,quote = F)


plot<-as.data.frame(Mid %>% group_by(SNP)  %>% dplyr::summarise(count=n()))

sum(plot$count)
binom.test(29,sum(5+12+11+1+4+29+11+10),1/8)
binom.test(65,sum(65+14+5+8),1/4)


pdf("A549_1.pdf",width = 3,height = 3)

plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+scale_y_continuous(limits = c(0,80)) +
  labs(x='',y='Mutation count')+ my_theme2+guides(fill=F)
dev.off()

#    A 8637
#    C 5338
#    G 5677
#    T 9253
fisher.test(matrix(c(29,11,(5677-29),(5338-11)),nrow=2))


#mutation rate

#C	260205+283168+273443
#A	378420+410749+394759
#G	239705+262490+254775
#T	293026+316413+303995
#
fisher.test(matrix(c(29,11,(239705+262490+254775),(260205+283168+273443)),nrow=2))
fisher.test(matrix(c(65,5,(260205+283168+273443),(239705+262490+254775)),nrow=2))


Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig1/A549.txt",header = T,sep="\t")
head(Mutation)
Mutation$SNP<-factor(Mutation$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210519/SciAdv_VS_Vero_rate.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=Mutation ,aes(x=SCV2_count,y=SciAdv,col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+scale_y_continuous(limits = c(0,300)) +scale_x_continuous(limits = c(0,52)) +
  scale_color_manual(values = Self1)+
  labs(x='SARS-CoV-2 mutation rate in Vero',y='SARS-CoV-2 mutation from Sci. Adv.') +
  my_theme2+guides(col=F,size=F)
dev.off()

cor.test(Mutation$SCV2_count,Mutation$SciAdv,method = "s")


ggplot(data=Mutation,aes(x=A549,y=SARS_CoV2_Polymorphism,col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+#scale_y_continuous(limits = c(-6,-4)) +scale_x_continuous(limits = c(-6,-4.5)) +
  scale_color_manual(values = Self1)+
  labs(x='# of SARS-CoV-2 mutation rate from A549',y='# of SARS-CoV-2 polymorphisms') +
  my_theme2+guides(col=F,size=F)
cor.test(log10(Mutation$SCV2_count),log10(Mutation$A549))


ggplot(data=Mutation ,aes(x=SCV2_Rate,y=A549_Rate,col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+
  scale_color_manual(values = Self1)+
  scale_x_continuous(limits = c(0,2.5*10^-5)) +
  labs(x='log10(SARS-CoV-2 mutation rate from Vero)',y='log10(SARS-CoV-2 mutation rate from A549)') +
  my_theme2+guides(col=F,size=F)

cor.test(log10(Mutation$SCV2_Rate),log10(Mutation$A549_Rate))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210519/A549_VS_Vero_count.pdf",height = 3,width = 3,useDingbats = F)

ggplot(data=Mutation ,aes(x=SCV2_count,y=A549,col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+scale_y_continuous(limits = c(0,80)) +scale_x_continuous(limits = c(0,52)) +
  scale_color_manual(values = Self1)+
  labs(x='SARS-CoV-2 mutation count from Vero',y='SARS-CoV-2 mutation count from A549') +
  my_theme2+guides(col=F,size=F)
dev.off()
cor.test(log10(Mutation$SCV2_count),log10(Mutation$A549))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210519/A549_rate.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=Mutation, aes(x=factor(SNP),y=A549_Rate*10^5,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+#scale_y_continuous(limits = c(0,25)) +
  labs(x='',y='Mutation rate (x10^-5)')+ 
  my_theme2+guides(fill=F)
dev.off()
t.test(c(4.22323e-06,1.01358e-05,1.34669e-05,1.22427e-06,1.20425e-05,1.09477e-05), mu = 3.83106e-05)
t.test(c(1.18250e-05,6.60528e-06,8.75816e-06), mu = 7.95773e-05)


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210716/Fig2D-1.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=Mutation ,aes(x=log10(SCV2_count),y=log10(SARS_CoV2_Polymorphism),col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+scale_y_continuous(limits = c(1,3.5)) +scale_x_continuous(limits = c(0.5,2)) +
  scale_color_manual(values = Self1)+
  labs(x='SARS-CoV-2 mutation  in Vero',y='SARS-CoV-2 mutation polymorphism')# +
  #my_theme2+guides(col=F,size=F)
dev.off()
cor.test(log10(Mutation$SCV2_count),log10(Mutation$SARS_CoV2_Polymorphism))

############################################################################################

#Fig1C

############################################################################################

Junction<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/20210609/Junction_Site_mismatch_consensus_reads.txt",stringsAsFactors = F))
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210609/Consensus_Mismatch.pdf",height = 3,width = 3,useDingbats = F)
colnames(Junction)<-c("Dis","Mis","Cover")
ggplot(data=Junction, aes(x=(Dis),y=Mis/Cover)) + #geom_abline(intercept=0.5, slope=0,color="red")+
  geom_bar(position="dodge", stat="identity",width = 0.5)+
  #theme_classic()+theme(panel.background=element_rect(colour='black'))+
  labs(x='Position(Junction site=0) ',y='Mismatch number')+#barplot_theme +
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+  
  scale_x_continuous(breaks = seq(-30,30,15),limits = c(-30,30))+  
  scale_y_continuous(limits = c(0,0.025))+my_theme2
dev.off()


###############################################################################################

#FigS1C

###############################################################################################
Junction<-read.table("/Dell/Dell13/shankj/projects/Cov/SCV2SJ.out.tab",sep="\t",header = F,stringsAsFactors = F)
head(Junction)
colnames(Junction)<-c("Chr","start","end","strand","intron_motif","annotated","uniquely_mapped","multiple_mapped","Overhang")
#Junction2<-separate(Junction,UMI,into = c('start','end'),sep = '-')[,1:3]
Junction$sum<-Junction$uniquely_mapped+Junction$multiple_mapped
nrow(filter(Junction))

cutoff<-2^5
for (i in c(1:nrow(Junction))){
  if (Junction[i,10]> cutoff){
    Junction[i,10]<- cutoff
  }
}

pdf("FigS1C.pdf",height = 4,width = 5.5)
jpeg("FigS1C.jpg",units = "cm",height = 8, width = 11,res = 300)
jpeg("/Dell/Dell13/shankj/projects/Cov/Plot/20210702/FigS1C_nolegend.jpg",units = "cm",height = 8, width = 8,res = 300)
ggplot() + 
  geom_point(data=Junction, aes(x=as.numeric(start),y=-as.numeric(end),color=log2(sum)),shape=15,size=0.2)+
  geom_point(data=filter(Junction,start==3108,end==28158),aes(x=as.numeric(start),y=-as.numeric(end)),color="red",size=2,shape=23)+
  scale_color_continuous(type = "viridis")+my_theme2+
  ylab("End")+xlab("Start")+guides(col=F)
dev.off()

mid<-Junction %>% group_by(sum) %>% dplyr::summarise(count=n())

sum(filter(mid,count<=20)$count)
