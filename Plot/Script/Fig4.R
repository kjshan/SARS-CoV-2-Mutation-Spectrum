theme_set(theme_bw()+theme(panel.grid=element_blank(),panel.border=element_rect(size=1,color="black")))
my_theme<-theme(axis.line.x=element_line(size=0,color="black"),axis.line.y=element_line(size=0,color="black"),
                axis.ticks=element_line(size=0.5,color="black"),axis.ticks.length=unit(0.05,"inches"),
                axis.title.x = element_text(size=10),axis.title.y = element_text(size=10),
                axis.text.x = element_text(angle = 0,hjust = 1,size=8,color="black"),
                axis.text.y =  element_text(size=10,color="black"),
                strip.text.x = element_text(size=10,face = "bold"),
                strip.background = element_rect(color = "black",size=1),
                legend.position = "none",
                legend.text = element_text(size=10),legend.title = element_text(size=10))

#########################################################################

setwd("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig4/")

#ACRseq-Cirseq
#Yeast

#############################################################################################################################################################################################################################
#CALL SNP
#Transcript Error(Cirseq)
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/SNP/",pattern = "*.TranscriptionError")
plot_list=list()
Filter<-data.frame()
Cirseq_plot<-data.frame()
i=4
head(Mid)
for (i in 1:length(temp)) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/SNP/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  colnames(Mutation)<-c("Strand","Chr","Site","Total_number","SNP","Number","Fraction","Reads","Ref_Reads")
  Mutation$SNP_Ref<-Mutation$Ref_Reads/(Mutation$Total_number)
  Mid<-filter(Mutation, !SNP %in% c('SNP'),!Chr %in% c("MT","L-BC-La","L-A","L-BC-2","20S_RNA_narnavirus","23S_RNA_narnavirus"), Fraction<=0.01, SNP_Ref>=0.99) %>% #Number==1,
    group_by(Chr,Site,SNP,Number)  %>% 
    dplyr::summarise(count=n()) %>% 
    filter(count==1)# ,!grepl('MT', Chr)
  name=sub(".TranscriptionError","",temp[i])
  Mid$sample<-rep(name,nrow(Mid))
  Filter<-rbind(Filter,as.data.frame(Mid))
  
  plotmid<-as.data.frame(Mid  %>% group_by(Chr,Site,SNP)  %>% dplyr::summarise())
  plot<-as.data.frame(filter(plotmid) %>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
  
  #plot<-as.data.frame(Mid %>% group_by(SNP)  %>% dplyr::summarise(Count=sum(Number)))
  plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
  
  plot_list[[temp[i]]]<-ggplot(data=plot, aes(x=factor(SNP),y=Count/1000,fill=SNP)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_fill_manual(values = Self)+
    labs(x='',y='Transcript error(x 1000)',title = name)+my_theme

  plot$Cirseq<-rep(name,nrow(plot))
  Cirseq_plot<-rbind(Cirseq_plot,plot)
  rm(Mid)
  rm(Mutation)
}
do.call(grid.arrange, plot_list)


P1<-Cirseq_plot
P2<-Cirseq_plot
Lynch_Me<-merge(P1,P2,by=c("SNP","Cirseq"))
head(Lynch_Me)
test<-data.frame()
for (i in unique(Lynch_Me$Cirseq)) {
  mid<-filter(Lynch_Me,Cirseq==i)
  mid<-corfun(mid$Count.x,mid$Count.y)
  rownames(mid)<-i
  test<-rbind(mid,test)
}




ctrl<-filter(Cirseq_plot,Cirseq=="ctrl")
wt<-filter(Cirseq_plot,Cirseq=="wt_yeast")

cor.test(ctrl$Count,
         wt$Count,
         method = "s")
#做检验

Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/wt_ctrl/ctrl.Mismatch_Match.mpileup.content100",sep="\t",header = T,stringsAsFactors = F)
#排除同一位点计算多次的情况，因为可能出现一个位点位于多个基因的情况
plotmid<-Mutation %>% filter(!Chr=="MT") %>% group_by(Chr,Site,count)  %>% dplyr::summarise(Count=n()) %>% filter(Count==1)
merge<-merge(plotmid,Mutation,by = c("Chr","Site","count"))

head(Base_pos_count)

#以基因组中>=100 reads覆盖的位点作为background，如果在p2-p8中有同一位点，只计算一次
Base_pos_count<-merge %>% group_by(Chr,Site,count,Base,Count)  %>% dplyr::summarise(Pos_number=n())  %>% #相同的位点只计算一遍
  group_by(Base) %>% dplyr::summarise(Count=n())

#Base   Count
# A     288263
# C     184247
# G     187704
# T     254766
fisher.test(matrix(c(20873,	10895,(187704-20873),(184247-10895)),nrow=2))

#TranscriptError(Bulk RNAseq)
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Yeast/SNP/",pattern = "*.TranscriptionError.Final")
plot_list=list()
mid<-data.frame()
test<-data.frame()
head(Mutation)
for (i in 1:length(temp)) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Yeast/SNP/",temp[i]),sep="\t",header = F,stringsAsFactors = F)
  colnames(Mutation)<-c("Chr","Site","SNP","SNP_reads_Number","SNP_reads_Frac","Reads","Ref_reads_Frac","Strand")
  name=sub(".TranscriptionError.Final","",temp[i])
  test<-Mutation %>% filter(SNP_reads_Number==1,SNP_reads_Frac<=0.01,Ref_reads_Frac>=0.99,!Chr=="MT") %>% group_by(SNP) %>% dplyr::summarise(count=n())#,Total_Fraction<=0.01 !grepl('MT', Chr),
  test$name<-rep(name,length(test$SNP))
  test$SNP<-factor(test$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
  mid<-rbind(mid,test)
}

bulk_yeast_transcript_error<-mid
rm(mid)


#bulk vs cirseq correlation
bulk_wt<-merge(bulk_yeast_transcript_error,wt,by = c("SNP"))
colnames(bulk_wt)
bulk_ctrl<-merge(bulk_yeast_transcript_error,ctrl,by = c("SNP"))

test<-data.frame()
for (i in unique(bulk_wt$name)) {
  sample<-filter(bulk_wt,name==i)#bulk_ctrl
  mid<-cor.test(sample$count,sample$Count,method = "s",alternative = "t")
  rho_P<-data.frame(sample=i,rho=as.numeric(mid$estimate), P_value=as.numeric(mid$p.value))
  test<-rbind(rho_P,test)
}

ggplot(data = test) +
  geom_point(mapping = aes(x = rho, y = -log10(P_value)))+
  my_theme



#####################################################################

#20S,L-A,L-BC-la mutation spectrum 
#Cirseq

#####################################################################

Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/wt_ctrl/wt_ctrl.Narnavirus.TranscriptionError",sep="\t",header = F,stringsAsFactors = F)
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/20S/wt_H2O2.Narnavirus.TranscriptionError",sep="\t",header = F,stringsAsFactors = F)
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/20S/wt_yeast.Narnavirus.TranscriptionError",sep="\t",header = F,stringsAsFactors = F)
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/SNP/E1103G_yeast.TranscriptionError",sep="\t",header = F,stringsAsFactors = F)
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/SNP/delRPB9_yeast.TranscriptionError",sep="\t",header = F,stringsAsFactors = F)

fisher.test(matrix(c(51, (725-51), 31, (740-31)), nrow=2))

#20S
#G	725
#A	445
#T	604
#C	740

fisher.test(matrix(c(565, (1107-565), 298, (985-298)), nrow=2))

#L-A
#G	1107	
#A	1296	
#T	1191	
#C	985	

fisher.test(matrix(c(52, (1101-52), 24, (860-24)), nrow=2))

#L-BC-La
#G	1101	
#A	1431	
#T	1223	
#C	860	



head(Mutation)
colnames(Mutation)<-c("Chr","Site","SNP","SNP_Number","SNP_Fraction","Reads","Ref_reads")

#E1103G and delRPB9
#colnames(Mutation)<-c("Strand","Chr","Site","Total_reads","SNP","SNP_Number","SNP_Fraction","Reads","Ref_reads") 
plotmid<-as.data.frame(Mutation %>% filter(SNP_Number>=1, Chr %in% c("L-BC-La","L-A","L-BC-2","20S_RNA_narnavirus","23S_RNA_narnavirus")) %>%  group_by(Chr,SNP)  %>% dplyr::summarise(count=n()))
head(plotmid)
#plot<-as.data.frame(filter(plotmid) %>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
plotmid$SNP<-factor(plotmid$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

#pdf("Cirseq.pdf",width = 2,height = 2)
plot_list=list()
mid<-data.frame()
head(Mutation)
for (i in unique(plotmid$Chr)) {
  #if (i %in% c('L-A')){#!i %in% c('L-BC-2')
  plot_list[[i]]<-ggplot(data=filter(plotmid,Chr %in% i), aes(x=factor(SNP),y=count ,fill=SNP)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_fill_manual(values = Self)+#facet_wrap(~Chr)+#facet_grid(SNP_Fraction_bin~Chr)+
    labs(x='',y='Mutation count',title = i)+my_theme
  #}
}
do.call(grid.arrange, plot_list)

levels(plotmid$SNP_Fraction_bin)

Total_Virus<-plotmid %>% group_by(SNP) %>% dplyr::summarise(Count=sum(count))
Cor<-merge(Total_Virus,ctrl,c("SNP"))
colnames(Cor)<-c("SNP","Virus","TE")
Cor$SNP<-factor(Cor$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Cor, aes(x=TE/10000,y=Virus/100,col=SNP)) + 
  geom_point(size=2.5)+
  scale_color_manual(values = Self1)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='SNV count from total RNA virus (x 100)',x='Transcript Error (x 10^4)')+ 
  my_theme
cor.test(Cor$Virus,Cor$TE)

#####################################################################

#PV and Hela correlation

#####################################################################
base<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/0915/Fragment/PV.txt",header = T,sep = "\t")
base$SNP<-factor(base$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
Self1<-c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99" )

pdf("PV_Hela.pdf",width =3,height = 2)

ggplot(data=base, aes(x=log10(Hela),y=log10(PV_JBC*0.000001),col=SNP)) + 
  geom_point()+
  scale_color_manual(values = Self1)+ #scale_x_continuous(limits = c(0,60)) +
  labs(x='log10(Hela transcription error frequency)',y='log10(PV mutation frequency)')+ 
  #scale_y_continuous(limits = c(0,300)) +
  theme_classic()+theme(panel.background=element_rect(colour='black'))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, colour='black'),axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, colour='black'))

dev.off()

cor.test(log10(base$Hela),log10(base$PV_JBC*0.000001),method = "s",alternative = "t")


ggplot(data=base, aes(x=SNP,y=(10^PV_Nature)*1000000,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='PV mutation frequency (x 10^-6)',title = "Poliovirus(+) From Nature \nMeasured by ImageJ" )+ 
  #scale_y_continuous(limits = c(0,300)) +
  my_theme

#####################################################################

#Ebola
Ebola_293T <- as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/Fig1/Ebola/293T.txt",stringsAsFactors = F))
Ebola_EpoNi <-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/Plot/Fig1/Ebola/EpoNi.txt",stringsAsFactors = F))

Ebola_293T$Type <- rep(c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C"),7)
Ebola_293T$Type<-factor(Ebola_293T$Type,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
Ebola_293T$name=paste0(Ebola_293T$Sample,Ebola_293T$Replication)

Ebola_EpoNi$Type <- rep(c("T > G","T > C","T > A","G > T","G > C","G > A","C > T","C > G","C > A","A > T","A > G","A > C"),7)
Ebola_EpoNi$Type<-factor(Ebola_EpoNi$Type,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
Ebola_EpoNi$name=paste0(Ebola_EpoNi$Sample,Ebola_EpoNi$Replication)

ggplot(data=Ebola_293T, aes(x=Type,y=log10(MLE),fill=Type)) + 
  geom_boxplot()+
  scale_fill_manual(values = Self)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='mutation frequency (log10)',title = "Ebola (-) from 293T" )+ 
  facet_grid(~Sample)+
  #scale_y_continuous(limits = c(0,300)) +
  my_theme

ggplot(data=Ebola_293T, aes(x=Type,y=log10(MLE),fill=Type)) + 
  geom_boxplot()+
  scale_fill_manual(values = Self)+ #scale_x_continuous(limits = c(0,60)) +
  labs(y='mutation frequency (log10)',title = "Ebola (-) from 293T" )+ 
  #facet_grid(Replication~Sample)+
  #scale_y_continuous(limits = c(0,300)) +
  my_theme

#############################################################################################################################################################################################################################

#Vero bulk RNA-seq

#############################################################################################################################################################################################################################
head(Mutation)
pdf("Vero_TranscriptionError.pdf",width = 3,height = 3)
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/BulkRNAseq/SNP/Vero/Vero.PE.overlap.V2.TranscriptionError3",sep="\t",header = T,stringsAsFactors = F)
Mutation$SNP_fragment_fraction<-1/Mutation$Fragment_Total_number
Mutation$Ref_fragment_fraction<-Mutation$Fragment_ref_number/Mutation$Fragment_Total_number
name=sub(".PE.overlap.V2.TranscriptionError3","",temp[i])
Mutation$Pos_SNP<-paste(Mutation$Chr,Mutation$Site,sep = ":")
mid<-as.data.frame(table(Mutation$Pos_SNP))
colnames(mid)<-c("Pos_SNP","count")
#head(mid)
Test1<-merge(Mutation,mid,by.x="Pos_SNP", by.y="Pos_SNP")
Vero<-Test1 %>% filter(count==1, !Chr=="MT",SNP_fragment_fraction<=0.01,Ref_fragment_fraction>=0.99) %>% group_by(SNP) %>% dplyr::summarise(count=n())  
Vero$SNP<-factor(Vero$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
ggplot(data=Vero ,aes(x=factor(SNP),y=count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='Mutation count',title = "Vero")+my_theme
dev.off()

Vero<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/0915/Fragment/Vero.txt",sep="\t",header = T,stringsAsFactors = F)
Vero$SNP<-factor(Vero$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Vero, aes(x=Vero_TE,y=Vero_SCV2,col=SNP)) + 
  geom_point(size=2.5)+
  scale_color_manual(values = Self1)+ #scale_x_continuous(limits = c(0,60)) +
  labs(x='Transcript Error in Vero cell',y='SARS-CoV-2 mutation')+ 
  scale_x_continuous(breaks = seq(0,300,100),limits = c(0,300))+
  scale_y_continuous(breaks = seq(0,75,10),limits = c(0,75)) +
  my_theme

cor.test(Vero$Vero_SCV2,Vero$Vero_TE,method = "s")


#做检验

head(Mutation)
  Mutation<-read.table(("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/BulkRNAseq/Vero/Vero.fragment.content100"),sep="\t",header = T,stringsAsFactors = F)
  #排除同一位点计算多次的情况，因为可能出现一个位点位于多个基因的情况
  plotmid<-Mutation %>% filter(!Chr=="MT") %>% group_by(Chr,Site,Count)  %>% dplyr::summarise(count=n()) %>% filter(count==1)
  merge<-merge(plotmid,Mutation,by = c("Chr","Site","Count"))


#以基因组中>=100 reads覆盖的位点作为background
Base_pos_count<-merge %>% group_by(Chr,Site,Base,Count)  %>% dplyr::summarise(Pos_number=n())  %>% #相同的位点只计算一遍
  group_by(Base) %>% dplyr::summarise(Count=n())

fisher.test(matrix(c(89,	25,(8285-89),( 10136-25)),nrow=2))


test<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Fig4.txt",header = T,sep="\t")
test$SNP<-factor(test$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=test, aes(x=log10(Vero_TE),y=log10(Hela_TE),col=SNP)) + 
  geom_point(size=2.5)+
  scale_color_manual(values = Self1)+ #scale_x_continuous(limits = c(0,60)) +
  labs(x='Transcript Error in Vero cell(log10)',y='Transcript Error in Hela cell(log10)')+ 
  my_theme

cor.test(test$Vero_TE,test$Hela_TE)


#TranscriptError(Bulk RNAseq)
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Lucas/",pattern = "*.PE.overlap.V2.TranscriptionError3")
plot_list=list()
mid<-data.frame()
HK293T_plot<-data.frame()
head(Mutation)
for (i in 1:length(temp)) {
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Lucas/",temp[i]),sep="\t",header = T,stringsAsFactors = F)
  Mutation$SNP_fragment_fraction<-1/Mutation$Fragment_Total_number
  Mutation$Ref_fragment_fraction<-Mutation$Fragment_ref_number/Mutation$Fragment_Total_number
  name=sub(".PE.overlap.V2.TranscriptionError3","",temp[i])
  Mutation$Pos_SNP<-paste(Mutation$Chr,Mutation$Site,sep = ":")
  mid<-as.data.frame(table(Mutation$Pos_SNP))
  colnames(mid)<-c("Pos_SNP","count")
  #head(mid)
  Test1<-merge(Mutation,mid,by.x="Pos_SNP", by.y="Pos_SNP")
  HK293T<-Test1 %>% filter(count==1, !Chr=="MT",SNP_fragment_fraction<=0.01,Ref_fragment_fraction>=0.99) %>% group_by(SNP) %>% dplyr::summarise(count=n())  
  HK293T$SNP<-factor(Vero$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
  plot_list[[temp[i]]]<-ggplot(data=HK293T ,aes(x=factor(SNP),y=count,fill=SNP)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_fill_manual(values = Self)+
    labs(x='',y='Mutation count',title = name)+my_theme
  
  HK293T$sample<-rep(name,nrow(HK293T))
  HK293T_plot<-rbind(HK293T_plot,HK293T)
  rm(mid)
  rm(Mutation)
}

do.call(grid.arrange, plot_list)
head(HK293T_plot)

HK293T_plot<-spread(HK293T_plot,sample,count)

test<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Fig4.txt",header = T,sep="\t")
test$SNP<-factor(test$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
test<-merge(test,HK293T_plot,c("SNP"))

test1<-data.frame()
for (i in 2:length(colnames(HK293T_plot))) {
  mid<-cor.test(test[,i],test$Hela_TE,method = "s")
  rho_P<-data.frame(sample=colnames(HK293T_plot)[i],rho=as.numeric(mid$estimate), P_value=as.numeric(mid$p.value))
  test1<-rbind(rho_P,test1)
}

ggplot(data = test1) +
  geom_point(mapping = aes(x = rho, y = -log10(P_value)))+
  my_theme+ geom_hline(yintercept=-log10(0.05),col="red",linetype="dashed")

#############################################################################################################################################################################################################################

#293T bulk RNA-seq

#############################################################################################################################################################################################################################

#Control 1
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Lucas/Control1.PE.overlap.TranscriptionError3",sep="\t",header = F,stringsAsFactors = F)
head(Count)
colnames(Count)<-c("Chr","Site","SNP","SNP_reads_Number","Total_number","Total_Fraction","Reads","Strand")

plotmid<-as.data.frame(Count  %>% filter(Total_Fraction<=0.01,!grepl('MT', Chr))  %>% group_by(Chr,Site,SNP)  %>% dplyr::summarise(count=n()))
plot<-as.data.frame(plotmid%>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("HEK293T_R1.pdf",width = 2,height = 2)
ggplot(data=plot, aes(x=factor(SNP),y=Count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='Mutation count')+my_theme
dev.off()

#  Ref Count
#   A  1019
#   C  1014
#   G   729
#   T   777

fisher.test(matrix(c(56,	28,(729-56),(1014-28)),nrow=2))



#Control 2
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig4/Lucas/Control2.PE.overlap.TranscriptionError3",sep="\t",header = F,stringsAsFactors = F)
head(Count)
colnames(Count)<-c("Chr","Site","SNP","SNP_reads_Number","Total_number","Total_Fraction","Reads","Strand")

plotmid<-as.data.frame(Count  %>% filter(Total_Fraction<=0.01,!grepl('MT', Chr))  %>% group_by(Chr,Site,SNP)  %>% dplyr::summarise(count=n()))
plot<-as.data.frame(plotmid%>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("HEK293T_R2.pdf",width = 2,height = 2)
ggplot(data=plot, aes(x=factor(SNP),y=Count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+
  labs(x='',y='Mutation count')+my_theme
dev.off()

#Ref Count
#   A   886
#   C   866
#   G   661
#   T   664
fisher.test(matrix(c(57,22,(661-57),(866-22)),nrow=2))

#############################################################################################################################################################################################################################

#Fig5 Average
GT_CA<-data.frame(Sample=c("Bat/Human","Monkey/Human"),GT_CA=c(0.341341433,1.073594329))
ggplot(data=GT_CA, aes(x=factor(Sample),y=GT_CA)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.5)+scale_y_continuous(breaks = seq(0,1.25,0.25),limits = c(0,1.25)) +
  #scale_fill_manual(values = Self)+
  labs(x='',y='The ratio of GT/AC odds ratio between two species')+  geom_hline(yintercept=1,col="red",linetype="dashed")+my_theme


#Fig5 point-to-point
GT_CA<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Fig5/test_point.txt",sep="\t",header = T,stringsAsFactors = F)
head(GT_CA)
mid<-c()
for (i in 1:9) {
  for (j in 10:14) {
      mid<-c(GT_CA$Odds_ratio[i]/GT_CA$Odds_ratio[j],mid)
  }
}

Bat<-data.frame(Odss=mid,Sample=rep("Bat/Human",length(mid)))

wilcox.test(GT_CA$Odds_ratio[1:9],GT_CA$Odds_ratio[10:14])

mid<-c()
for (i in 15:21) {
  for (j in 10:14) {
    
    mid<-c(GT_CA$Odds_ratio[i]/GT_CA$Odds_ratio[j],mid)
    
    
  }
}

Vero<-data.frame(Odss=mid,Sample=rep("Vero/Human",length(mid)))
wilcox.test(GT_CA$Odds_ratio[15:21],GT_CA$Odds_ratio[10:14])


GT_CA<-rbind(Bat,Vero)

ggplot(data=GT_CA, aes(x=factor(Sample),y=Odss)) + 
  geom_boxplot() +
  #scale_fill_manual(values = Self)+
  labs(x='',y='The ratio of GT/AC odds ratio between two species')#+  geom_hline(yintercept=1,col="red",linetype="dashed")+my_theme


