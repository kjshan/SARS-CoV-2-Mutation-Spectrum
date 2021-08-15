theme_set(theme_bw()+theme(panel.grid=element_blank(),panel.border=element_rect(size=1,color="black")))
my_theme<-theme(axis.line.x=element_line(size=0,color="black"),axis.line.y=element_line(size=0,color="black"),
                axis.ticks=element_line(size=0.5,color="black"),axis.ticks.length=unit(0.05,"inches"),
                axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
                axis.text.x = element_text(angle = 45,hjust = 1,size=8,color="black"),
                axis.text.y =  element_text(size=10,color="black"),
                strip.text.x = element_text(size=10,face = "bold"),
                strip.background = element_rect(color = "black",size=1),
                legend.position = "none",
                legend.text = element_text(size=10),legend.title = element_text(size=10))

#########################################################################

setwd("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/")

#########################################################################

#PV https://www.jbc.org/content/289/43/29531.long#fn-group-1
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/PV.txt",sep="\t",header = T,stringsAsFactors = F)
Count$SNP<-factor(Count$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/PV_Mutation.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=Count, aes(x=factor(SNP),y=PV_JBC/100,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='Mutation freq. (x10-4)')+my_theme2+guides(fill=F)
dev.off()

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/PV_Mutation_Nature.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=Count, aes(x=factor(SNP),y=(10^PV_Nature)*10^5,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='Mutation freq. (x10-5)')+my_theme2+guides(fill=F)
dev.off()

pdf("Hela_PV.pdf",width = 2,height = 2)
ggplot(data=Count, aes(x=log10(Hela*0.000001),y=log10(PV*0.000001),col=SNP)) + 
  geom_point()+
  scale_color_manual(values = Self1)+
  labs(x=' log10(Hela Mutation freq.)',y='log10(PV Mutation freq.)')+my_theme
dev.off()

cor.test(Count$Hela,Count$PV,method = "spearman")


Ebola_293T <- as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ebola/293T.txt",stringsAsFactors = F))
Ebola_293T$Type <- rep(c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C"),7)
Ebola_293T$Type<-factor(Ebola_293T$Type,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
head(Ebola_293T)
Ebola_293T_plot<-summarySE(filter(Odds_ratio,Odds_ratio!=Inf), measurevar="Odds_ratio", groupvars=c("Chr","Group"))

pdf("Ebola_293T.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=Ebola_293T, aes(x=Type,y=MLE,fill=Type)) + 
  geom_boxplot()+
  scale_fill_manual(values = Self)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='Mutation frequency',title = "Ebola (-) from 293T" )+ 
  #facet_grid(~Sample)+
  #scale_y_continuous(limits = c(0,300)) +
  my_theme
dev.off()




#########################################################################

#MERS-CoV
#polymorphism 
#Ref:NC_019843
#number:258

#########################################################################
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/MERS/MERS_Astr_Align.SNP",sep="\t",header = F,stringsAsFactors = F)

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

fisher.test(matrix(c(73,	22,(6303-73),(6085-22)),nrow=2))


#########################################################################

#H1N1
#polymorphism 
#Ref:EPI_ISL_29801 (ssRNA - ) 

#########################################################################
setwd("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/")
tmp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Genome/",pattern = "*.txt")
tmp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/ChangShuo/FluA/flu_Seq/SNP/Astr/",pattern = "*.SNP")

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

pdf("H1N1_polymorphism_Astr.pdf",width = 2,height = 2)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme
dev.off()
fisher.test(matrix(c(122,	60,(2456-122),(3057-60)),nrow=2))

##########################################################################################


#Measles
#polymorphism 
#Ref:FJ161211 (ssRNA - ) 
#number:144

#########################################################################
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/Measles_FJ161211.SNP",sep="\t",header = F,stringsAsFactors = F)
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/Astr/Measles_Astr_Align.SNP",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

#每个位点的SNP只算1???
SNP<-as.data.frame(Count %>% group_by(SNP,Site)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
plot$SNP <- c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C")
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#Astr
#C	3731	23.47
#T	4646	29.23
#A	3729	23.46
#G	3788	23.83

pdf("Measles_polymorphism.pdf",width = 2,height = 2)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme
dev.off()

fisher.test(matrix(c(111,	79,(3788-111),(3731-79)),nrow=2))

##########################################################################################


#HPIV3
#polymorphism 
#Ref:FJ455842 (ssRNA - ) 
#number:338

#########################################################################
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig3/HPIV_3_FJ455842.SNP",sep="\t",header = F,stringsAsFactors = F)
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/Astr/HPIV3_Astr_Align.SNP",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

#每个位点的SNP只算1???
SNP<-as.data.frame(Count %>% group_by(SNP,Site)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
plot$SNP <- c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C")
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#Astr
#C	2839	18.36
#T	5940	38.42
#A	4088	26.44
#G	2595	16.78


pdf("HPIV3_polymorphism.pdf",width = 2,height = 2)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='SNV count')+my_theme
dev.off()

fisher.test(matrix(c(103,	56,(2595-103),(2839-56)),nrow=2))
#########################################################################################


#herpesvirus(paozhen) 
#polymorphism 
#Ref:NC_006273 (dsDNA,+ ) 
#number:299

#########################################################################################

#Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/Astr_Human_betaherpesvirus_5.SNP",sep="\t",header = F,stringsAsFactors = F)
#Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/USA_28/betaherpesvirus_Astr_Align.SNP",sep="\t",header = F,stringsAsFactors = F)
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/SNP2/",pattern = "*.SNP")
Count<-data.frame()
for(i in 1:length(temp)){
  Count1<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/SNP2/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  colnames(Count1)<-c("Site","SNP","Sample")
  name=sub(".SNP","",temp[i])
  Count1$Name<-rep(name,nrow(Count1))
  Count<-rbind(Count,Count1)
}


#每个位点的SNP只算1???
SNP<-as.data.frame(Count %>% group_by(SNP,Site,Name)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>=2)%>% group_by(SNP)  %>% dplyr::summarise(count=n()))
#plot$SNP <- c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C")
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#Astr
#g	68526	29.08
#a	50307	21.35
#t	49430	20.98
#c	67383	28.60

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/betaherpesvirus_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(limits = c(0,12000)) +
  labs(x='',y='SNV count')+my_theme+my_theme2+guides(fill=F)
dev.off()

fisher.test(matrix(c(650,	736,(68526-650),(67383-736)),nrow=2))
binom.test(650,(650+736),68526/(67383+68526))

#########################################################################################


#Vaccinia(niudou) 
#polymorphism 
#Ref:NC_006998 (dsDNA,+ ) 
#number:96

#########################################################################################

#Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/Astr_Human_betaherpesvirus_5.SNP",sep="\t",header = F,stringsAsFactors = F)
#Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA/USA_28/betaherpesvirus_Astr_Align.SNP",sep="\t",header = F,stringsAsFactors = F)
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA_Vaccinia/SNP/",pattern = "*.SNP")
Count<-data.frame()
for(i in 1:length(temp)){
  Count1<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/Genome/NegativeStrand/Kejia/dsDNA_Vaccinia/SNP/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  colnames(Count1)<-c("Site","SNP","Sample")
  name=sub(".SNP","",temp[i])
  Count1$Name<-rep(name,nrow(Count1))
  Count<-rbind(Count,Count1)
}


#每个位点的SNP只算1???
SNP<-as.data.frame(Count %>% group_by(SNP,Site,Name)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>=2)%>% group_by(SNP)  %>% dplyr::summarise(count=n()))
#plot$SNP <- c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C")
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#Astr
#g	32385	16.63
#a	64901	33.33
#t	64858	33.31
#c	32567	16.73

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Vaccinia_polymorphism.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+scale_y_continuous(breaks = c(seq(0,750,250)),limits = c(0,750)) +
  labs(x='',y='SNV count')+my_theme2+guides(fill=F)#my_theme
dev.off()

fisher.test(matrix(c(120,	124,(32385-120),(32567-124)),nrow=2))




########################################################################################


test <-
  array(c(122,60,(2456-122),(3057-60), #H1N1         
          111,79,(3788-111),(3731-79), #Mealse
          103,56,(2595-103),(2839-56), #HPIV3
          58,37,(4024-58),(3722-37), #Ebola         
          125,49,(2674-125),(2367-49), #RSV_A  
          68,17,(2725-68),(2400-17), #RSV_B         
          27,29,(3373-27),(3093-29), # Mumps_genotypeG
          73,	22,(6303-73),(6085-22),#MERS
          877,127,(5093-877),(4723-127)#SARS-CoV-2
  ),
  dim = c(2, 2, 9),
  dimnames = list(
    #Mutation = c("Virus", "TranscriptError"),
    C = c("C>A", "C"),
    G = c("G>T", "G"),
    Strain = c("H1N1", "Mealse","HPIV3","Ebola","RSV_A","RSV_B","Mumps_genotypeG","MERS","SARS-CoV-2")))
test
mantelhaen.test(test,exact = T)
###############################################################################################################

#Cirseq
#mRNA error count and rate

###############################################################################################################
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",pattern = "*.TranscriptionError")
plot_list=list()
Filter<-data.frame()
i=1
head(Mutation)
for (i in 1:7) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  
  colnames(Mutation)<-c("Strand","Chr","Site","Total_reads","SNP","SNP_Number","Fraction","Ref_Reads")
  Mutation$SNP_Ref<-Mutation$Ref_Reads/(Mutation$Total_reads)
  Mid<-as.data.frame(filter(Mutation, !Chr %in% c("MT","Y"),Fraction<=0.01,SNP_Ref>=0.99) %>% 
                       group_by(Chr,Site,SNP,SNP_Number,Total_reads)  %>% 
                       dplyr::summarise(count=n()) %>% 
                       filter(count==1))
  name=sub(".TranscriptionError","",temp[i])
  
  #ATCG<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/",name,".ATGC"),header = F,stringsAsFactors = F)
  #colnames(ATCG)<-c("base","Num")
  Mid$sample<-rep(name,nrow(Mid))
  Filter<-rbind(Filter,Mid)
  
  plot<-as.data.frame(filter(Mid) %>% group_by(SNP)  %>% dplyr::summarise(Count=n()))#sum(SNP_Number)
  plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
  
  plot_list[[temp[i]]]<-ggplot(data=plot, aes(x=factor(SNP),y=Count/10000,fill=SNP)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_fill_manual(values = Self)+
    labs(x='',y='mRNA error(x 10^-4)',title = paste0(name,"\t mRNA"))+my_theme
  
  rm(Mid)
  rm(Mutation)
}
do.call(grid.arrange, plot_list)

tail(Filter)
head(plotmid)
plotmid<-as.data.frame(Filter  %>% group_by(Chr,Site,SNP)  %>% dplyr::summarise(count=n(),SNP_count=sum(SNP_Number)))
head(Filter)
#filter somatic mutation (filter(plotmid,count==1)) 
#and add sample name to data frame
merge<-merge(filter(plotmid,count==1),Filter[,c(1,2,3,4,5,7)],by=c("Chr","Site","SNP"))
head(merge)

#mRNA error in each sample
plot<-as.data.frame(merge %>% group_by(SNP,sample)  %>% dplyr::summarise(Count=n(),SNP_Number=sum(SNP_count)))
#total mRNA error
plot<-as.data.frame(merge %>% group_by(SNP)  %>% dplyr::summarise(Count=n(),SNP_Number=sum(SNP_count)))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Cirseq_Hela.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=plot, aes(x=factor(SNP),y=Count/10^4,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+#facet_grid(~sample)+
  labs(x='',y='mRNA error(x10^4)')+my_theme
dev.off()


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Cirseq_Hela2.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=plot, aes(x=factor(SNP),y=Count/10^4,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+#facet_grid(~sample)+
  labs(x='',y='mRNA error(x10^4)')+my_theme2+guides(fill=F)
dev.off()

#############################################################
#mutation rate

temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",pattern = "*.ATGC")
Filter<-data.frame()
i=1
head(Mutation)
for (i in 1:7) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",temp[i]),sep="\t",header = T,stringsAsFactors = F)
  #only exists in 1 gene
  Count<-as.data.frame(Mutation %>% group_by(Chr,Site,Base,count)  %>% dplyr::summarise(Count=n())) %>% filter(Count==1)
  name=sub(".ATGC","",temp[i])
  Count$sample<-rep(name,nrow(Count))
  Filter<-rbind(Filter,Count)
}

head(Filter)
Count<-Filter %>% group_by(Base,sample)  %>% dplyr::summarise(Count=sum(count))

colnames(Count)<-c("Ref","sample","Ref_total")
head(Count)
head(plot)
plot<-plot %>% 
  separate(SNP,into = c('Ref','Alt'),sep = ' > ')

plot$SNP<-paste0(plot$Ref," > ",plot$Alt)
rate<-merge(plot,Count,by=c("Ref","sample"))
head(rate)
rate$rate<-rate$SNP_Number/rate$Ref_total
rate$SNP<-factor(rate$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=rate, aes(x=factor(SNP),y=rate*10^5,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+facet_grid(~sample)+
  labs(x='',y='mRNA error rate(x10^-5)')+my_theme

rate %>% group_by(SNP) %>% dplyr::summarise(Rate=median(rate))

#fisher test
Count<-Filter %>% group_by(Base)  %>% dplyr::summarise(Count=n())
#Base    Count
# A     1062746
# C     1126421
# G     1169719
# T     1026427

fisher.test(matrix(c(11117,	8980,(1169719-11117),(112642-8980)),nrow=2))

#############################################################

Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Median_PV.txt",sep="\t",header = T,stringsAsFactors = F)
head(Mutation)
Mutation$SNP<-factor(Mutation$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Cirseq_Hela_PV.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data = Mutation) +
  geom_point(mapping = aes(x = log10(Hela), y = PV_Nature,color=SNP))+
  scale_color_manual(values = Self1)+
  labs(x='Hela mRNA error(median)',y='PV Nature(median)')#+ geom_abline(intercept=0, slope=1,color="black")
dev.off()

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig3/Cirseq_Hela_PV2.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data = Mutation) +
  geom_point(mapping = aes(x = log10(Hela), y = PV_Nature,color=SNP))+
  scale_color_manual(values = Self1)+
  labs(x='Hela mRNA error(median)',y='PV Nature(median)')+ my_theme2+guides(col=F)
dev.off()


cor.test(Mutation$Hela,Mutation$PV_Nature,method = "s")
