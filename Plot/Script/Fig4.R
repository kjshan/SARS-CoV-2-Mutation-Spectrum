#######################################

#ARC-seq

#######################################
Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/WT_Ctrl.noIndel.Match.mRNA.errorQ30",sep="\t",header = F,stringsAsFactors = F)
colnames(Mutation)<-c("Strand","Chr","Site","Total_reads","SNP","Number","Fraction","Ref_Reads")
head(Mutation)
ggplot(data=Mutation, aes(x=Fraction)) + 
  geom_vline(xintercept=0.05,col="black",linetype="dashed")+
  geom_density()+  
  geom_vline(xintercept=0.01,col="red",linetype="dashed")+
  labs(x='SNP fraction',y='')+my_theme

#mRNA
Filter<-filter(Mutation, !Chr %in% c("Mito","20S_RNA_narnavirus","L-A","L-BC-La","23S_RNA_narnavirus"), Fraction<=0.01, Ref_Reads/Total_reads>=0.99)# ,!grepl('MT', Chr)
Filter<-filter(Mutation, Chr %in% c("20S_RNA_narnavirus"), Fraction<=0.01, Ref_Reads/Total_reads>=0.99)# ,!grepl('MT', Chr)

plot<-as.data.frame(Filter %>% group_by(SNP)  %>% dplyr::summarise(Count=sum(Number)))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/WT_Ctrl_mRNA.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=Count*0.0001 ,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+#scale_y_continuous(limits = c(0,8)) +
  labs(x='',y='mRNA Error (x10^4)')+my_theme2+guides(fill=F)

ggplot(data=plot, aes(x=factor(SNP),y=Count ,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+#scale_y_continuous(limits = c(0,8)) +
  labs(x='',y='20S RNA virus')+my_theme2+guides(fill=F)
dev.off()



#20S RNA virus
Filter<-filter(Mutation, Chr %in% c("20S_RNA_narnavirus"), Fraction<=0.01, Ref_Reads/Total_reads>=0.99)# ,!grepl('MT', Chr)

plot<-as.data.frame(Filter %>% group_by(SNP)  %>% dplyr::summarise(Count=sum(Number)))
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/WT_Ctrl_20S.pdf",width = 3,height = 3)
ggplot(data=plot, aes(x=factor(SNP),y=Count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+#scale_y_continuous(limits = c(0,8)) +
  labs(x='',y='20S RNA virus mutation')+my_theme2+guides(fill=F)
dev.off()

# A     302147
# C     492118
# G     478751
# T     417745

fisher.test(matrix(c(46, (478751-46), 22, (492118-22)), nrow=2))
fisher.test(matrix(c(69, (478751-69), 157, (492118-157)), nrow=2))

#mRNA VS RNA virus

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210716/Fig4D_mRNA_20S.pdf",height = 3,width = 3,useDingbats = F)

Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/Ctrl.mRNA.Virus",sep="\t",header = T,stringsAsFactors = F)
Mutation$SNP<-factor(Mutation$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(Mutation, aes(x= log10(X20S_Count), y = log10(mRNA_Count), col = SNP)) +
  geom_point(size=3) +
  scale_color_manual(values = Self1)+ 
  scale_y_continuous(limits = c(3,5))+
  scale_x_continuous(limits = c(0.5,2.5))+
  labs(x='log10(20S RNA virus mutation count)',y='log10(mRNA errors count)')#+#,title = paste0("BY4742")
  #my_theme2+guides(col=F)
dev.off()

cor.test(log10(Mutation$X20S_Count),log10(Mutation$mRNA_Count),method = "p")

#######################################

#CirSeq 

#######################################


Lynch<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Lynch/faster/20160621_4WT_1/20160621_4WT_1.candidates-refined.clean-demult",sep="\t",header = F)
head(Lynch)
SUB<-filter(Lynch, V1 %in% c("BASE_SUB"),!V10 %in% c("NO_GENE","transcript:L-BC-La_mRNA","transcript:L-A_mRNA","transcript:L-BC-2_mRNA"),grepl("mRNA",V10))
head(SUB)
SUB$SNP<-paste0(toupper(SUB$V5)," > ",toupper(SUB$V6))
plot<-as.data.frame(SUB %>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
plot
write.table(plot,"/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/4WT_1.txt",row.names = F,quote = F,sep="\t")
plot<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/4WT_1.txt",sep="\t",header = T)
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

Mutation
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210625/S7A-2.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=plot, aes(x=factor(SNP),y=Count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+ scale_y_continuous(limits = c(0,150))+
  labs(x='',y='mRNA error')#+my_theme2+guides(fill=F)
dev.off()
fisher.test(matrix(c(30,3,(35397808-30),(34239236-3)),nrow=2))
fisher.test(matrix(c(33,140,(35397808-33),(34239236-140)),nrow=2))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210716/S5B-2.pdf",height = 3,width = 3,useDingbats = F)

Mutation<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/STAR/WT_Ctrl/Ctrl.mRNA.Virus",header = T,sep = "\t")
Mutation$SNP<-factor(Mutation$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

ggplot(data=Mutation,aes(x=log10(Count_Lynch ),y=log10(mRNA_Count),col=SNP,size=3) )+ #y=log10(as.numeric(Error_Fraction)+10^(-11))
  geom_point()+#scale_y_continuous(limits = c(-6,-4)) +scale_x_continuous(limits = c(-6,-4.5)) +
  scale_color_manual(values = Self1)+
  labs(x='Lynch',y='ARC-seq') +  
  scale_y_continuous(limits = c(3,5.5))+
  scale_x_continuous(limits = c(0,2.5))#+
  #my_theme2+guides(col=F,size=F)
dev.off()

mid=log10(Mutation$Count_Lynch)
mid[5]<- NA
cor.test(mid,log10(Mutation$mRNA_Count),use="na.or.complete")


