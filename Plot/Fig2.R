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


setwd("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/")
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/SARS_CoV_2_GASID/SCV2_Align.SNP2",sep="\t",header = F,stringsAsFactors = F)

colnames(Count)<-c("Site","SNP","Sample")
head(Count)

#count subtitutions only once in the same position
#and >=2 genome support the subtitution
SNP<-as.data.frame(Count %>% group_by(SNP,Site)  %>% dplyr::summarise(count=n()))
head(SNP)
plot<- as.data.frame(SNP %>% filter(count>1) %>% group_by(SNP)  %>% dplyr::summarise(count=n()))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/SCV2_polymorphism_Hu-1.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=count/100,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self1)+
  labs(x='',y='SNV count(x100)')+my_theme2+guides(fill=F)
dev.off()

fisher.test(matrix(c(1168,174,(5093-1168),(4723-174)),nrow=2))

#########################################################################

#Fig2B correlation between polymorphism and mutation

SCV2<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/VeroSNP_Polymorphism.txt",header = T,sep="\t")
SCV2$SNP<-factor(SCV2$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

head(SCV2)


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/Fig2B.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=SCV2, aes(x=SCV2_count,y=Polymorphism/100,col=SNP,size=3)) + 
  geom_point()+
  scale_color_manual(values = Self1)+
  labs(x='SARS-CoV-2 mutation from Vero',y='SNV count (x100)')+ 
  my_theme2+guides(col=F,size=F)
dev.off()

cor.test(SCV2$SCV2_count,SCV2$Polymorphism,method = "s")


pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/SCV2_VS_JBC1.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=SCV2, aes(x=log10(SCV2_Rate),y=log10(JBC/1000000),col=SNP,size=3)) + 
  geom_point()+
  scale_color_manual(values = Self1)+
  labs(x='log10(SARS-CoV-2 mutation frequency)',y='log10(Poliovirus mutation frequency)')+ 
  scale_x_continuous(limits = c(-6,-4))+
  scale_y_continuous(limits = c(-7,-3))+#my_theme1
  my_theme2+guides(col=F,size=F)
dev.off()

cor.test(log10(SCV2$SCV2_Rate),log10(SCV2$JBC/1000000),method = "s")

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig2/SCV2_VS_Nature.pdf",width = 3,height = 3,useDingbats = F)
ggplot(data=SCV2, aes(x=log10(SCV2_Rate),y=PV_Nature,col=SNP,size=3)) + 
  geom_point()+
  scale_color_manual(values = Self1)+
  labs(x='log10(SARS-CoV-2 mutation frequency)',y='log10(Poliovirus mutation frequency)')+ 
  scale_x_continuous(limits = c(-6,-4))+
  scale_y_continuous(limits = c(-6,-3))+#my_theme1
  my_theme2+guides(col=F,size=F)
dev.off()

cor.test(log10(SCV2$SCV2_Rate),SCV2$PV_Nature,method = "s")
#########################################################################

#############################################################################################################################################################################################################################

#Cirseq
#Hela cell
###############################################################################################################
#mRNA error

temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",pattern = "*.TranscriptionError")
plot_list=list()
Filter<-data.frame()
i=1
head(Mutation)
for (i in 1:7) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  
  colnames(Mutation)<-c("Strand","Chr","Site","Total_reads","SNP","SNP_Number","Fraction","Ref_Reads")
  Mutation$SNP_Ref<-Mutation$Ref_Reads/(Mutation$Total_reads)
  # plot_list[[temp[i]]]<-ggplot(data = Mutation,aes(x = Fraction)) +
  #   geom_density()+#facet_wrap(~SNP)+
  #   geom_vline(xintercept=0.01,col="red",linetype="dashed")+
  #   geom_vline(xintercept=0.05,col="blue",linetype="dashed")+
  #   labs(title = temp[i])
  
  Mid<-as.data.frame(filter(Mutation, !Chr %in% c("MT","Y"),Fraction<=0.01,SNP_Ref>=0.99) %>% #,Fraction<=0.01,SNP_Ref>=0.99
                       group_by(Chr,Site,SNP,SNP_Number,Total_reads)  %>%
                       dplyr::summarise(count=n()) %>%
                       filter(count==1))
  name=sub(".TranscriptionError","",temp[i])
  
  #ATCG<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/",name,".ATGC"),header = F,stringsAsFactors = F)
  #colnames(ATCG)<-c("base","Num")
  Mid$sample<-rep(name,nrow(Mid))
  Filter<-rbind(Filter,Mid)
  
  plot<-as.data.frame(filter(Mid) %>% group_by(SNP)  %>% dplyr::summarise(Count=sum(SNP_Number)))#sum(SNP_Number)
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
  scale_fill_manual(values = Self)+facet_grid(~sample)+
  labs(x='',y='mRNA error(x10^4)')+my_theme
dev.off()

#############################################################
#mutation rate

temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",pattern = "*.ATGC20")
Filter<-data.frame()
i=1
head(Mutation)
for (i in 1:7) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Human/SNP/20210507/Q20/",temp[i]),sep="\t",header = T,stringsAsFactors = F)
  #only exists in 1 gene
  Count<-as.data.frame(Mutation %>% group_by(Chr,Site,Base,count)  %>% dplyr::summarise(Count=n())) %>% filter(Count==1)
  name=sub(".ATGC20","",temp[i])
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
# Base    Count
# 1 A     6259950
# 2 C     6424358
# 3 G     6662586
# 4 T     6035698
fisher.test(matrix(c(22064,	20382,(6662586-22064),(6424358-20382)),nrow=2))

#############################################################

