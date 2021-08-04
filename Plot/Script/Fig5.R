#######################################

#CirSeq 

#######################################

Lynch<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Lynch/4WT_1.txt",sep="\t",header = T)

SUB<-filter(lynch, V1 %in% c("BASE_SUB"),!V10 %in% c("NO_GENE","transcript:L-BC-La_mRNA","transcript:L-A_mRNA","transcript:L-BC-2_mRNA"),grepl("mRNA",V10))
head(SUB)
SUB$SNP<-paste0(toupper(SUB$V5)," > ",toupper(SUB$V6))
plot<-as.data.frame(SUB %>% group_by(SNP)  %>% dplyr::summarise(Count=n()))
plot
plot$SNP<-factor(plot$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig5/4WT_1.pdf",height = 3,width = 3,useDingbats = F)
ggplot(data=plot, aes(x=factor(SNP),y=Count,fill=SNP)) + 
  geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
  scale_fill_manual(values = Self)+ scale_y_continuous(limits = c(0,150))+
  labs(x='',y='mRNA error')+my_theme2+guides(fill=F)
dev.off()
fisher.test(matrix(c(30,3,(35397808-30),(34239236-3)),nrow=2))

##################################################################################

#Consensus VS Inconsensus

##################################################################################

Consensus<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/BY4742/Ctrl1/Consensus_mismatch.txt",stringsAsFactors = F))
Inconsensus<-as.data.frame(fread("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/BY4742/Ctrl1/Unconsensus_mismatch.txt",stringsAsFactors = F))

ggplot()+
  geom_point(data=Consensus,aes(x=V1,y=log10(V4)),col="black")+
  geom_point(data=Inconsensus,aes(x=V1,y=log10(V4)),col="red")+
  scale_y_continuous(breaks = seq(-4,0,1),limits = c(-4,0))+
  my_theme2+guides(fill=F)



#######################################

#H2O2 VS Ctrl

#######################################
temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/bin2/",pattern = "*.SNP$")

temp <- list.files(path = "/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/data/SNP/A364A/",pattern = "*.SNP")
mid<-data.frame()
for (i in 1:length(temp)) {
  #Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/bin2/",temp[i]),sep="\t",header = T,stringsAsFactors = F)
  Mutation<-read.table(paste0("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/H2O2_Test/Clean/STAR/data/SNP/A364A/",temp[i]),sep="\t",header = T,stringsAsFactors = F)
  name=sub("\\d.Strand.SNP","",temp[i])
  Mutation$name<-rep(name,nrow(Mutation))
  mid<-rbind(mid,Mutation)
}
head(mid)


test<-data.frame()
for (i in unique(mid$name)) {
  Ctrl <- summarySE(filter(mid,name==i), measurevar="Mutation_Freq", groupvars=c("name","Chr","SNP"))
  colnames(Ctrl)<-c("Name","Chr","SNP",paste0(i,"N"),paste0(i,"Mutation_Freq"),paste0(i,"sd"),paste0(i,"se"),paste0(i,"ci"))
  Ctrl<-tidyr::spread(Ctrl, Name,paste0(i,"Mutation_Freq"))
  
  if (nrow(test)==0) {
    test<-rbind(test,Ctrl)
  }else{
    test<-merge(test,Ctrl,by=c("Chr","SNP"))
  }
}
test[is.na(test)] <- 0

write.table(test,"/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/bin2/BY4742_Freq.txt",row.names = F,quote = F,sep="\t")
write.table(test,"/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/A364A/A364_Freq.txt",row.names = F,quote = F,sep="\t")


#######################################

#CirSeq VS ours

#######################################

Lynch<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Cirseq_Yeast/Lynch/4WT_1.txt",sep="\t",header = T)

test<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/SNP/bin2/BY4742_Freq.txt",header = T,sep="\t")
colnames(test)
Lynch_A364A<-merge(Lynch,filter(test,Chr=="mRNA"),by="SNP")
test<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/A364A/A364_Freq.txt",header = T,sep="\t")
Lynch_A364A<-merge(Lynch_A364A,filter(test,Chr=="mRNA"),by="SNP")
Lynch_A364A$SNP<-factor(Lynch_A364A$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))


lm2 <- lmodel2(log10(Lynch_A364A$Freq_Lynch)~log10(Lynch_A364A$`Ctrl`))
ggplot(Lynch_A364A, aes(x= log10(`Ctrl`), y = log10(Freq_Lynch), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  scale_color_manual(values = Self1)+  
  scale_x_continuous(limits = c(-6,-3.9))+
  scale_y_continuous(limits = c(-8,-5))+
  labs(x='log10(mRNA errors Frequency from ours)',y='log10(mRNA errors Frequency from Lynch)')+#,title = paste0("BY4742")
  my_theme2+guides(col=F)

cor.test(Lynch_A364A$Freq_Lynch,Lynch_A364A$`Ctrl`,method = "s")#,method = "s"
cor.test(Lynch_A364A$Freq_Lynch,Lynch_A364A$`Ctrl`)#,method = "s"


ggplot(Lynch_A364A, aes(x= log10(`A364.Ctrl`), y = log10(Freq_Lynch), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  #geom_abline(intercept=0, slope=1,col="grey",linetype="dashed")+
  #geom_errorbar(aes(xmin = log10(`A364.Ctrl` - `A364.Ctrlse`), xmax = log10(`A364.Ctrl` + `A364.Ctrlse`)), 
  #              position = position_dodge(0.5),width=0,col='black',size=.5)+#facet_grid(~Chr)+
  scale_color_manual(values = Self1)+  
  scale_x_continuous(limits = c(-6,-3.5))+
  scale_y_continuous(limits = c(-8,-5))+
  labs(x='log10(mRNA errors Frequency from ours)',y='log10(mRNA errors Frequency from Lynch)')+
  my_theme2+guides(col=F)#+ 
#geom_abline(intercept=as.numeric(lm2$regression.results[3,2]),
#            slope=as.numeric(lm2$regression.results[3,3]),color='brown')

cor.test(Lynch_A364A$Freq_Lynch,Lynch_A364A$`A364.Ctrl`)#,method = "s"
cor.test(Lynch_A364A$Freq_Lynch,Lynch_A364A$`A364.Ctrl`,method = "s")#

lm2 <- lmodel2(log10(Lynch_A364A$Ctrl)~log10(Lynch_A364A$`A364.Ctrl`))
ggplot(Lynch_A364A, aes(x= log10(`A364.Ctrl`), y = log10(Ctrl), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  #geom_abline(intercept=0, slope=1,col="grey",linetype="dashed")+
  geom_errorbar(aes(ymin = log10(Ctrl - Ctrlse), ymax = log10(Ctrl + Ctrlse)), 
                position = position_dodge(0.5),width=0,col='black',size=.5 )+
  geom_errorbar(aes(xmin = log10(`A364.Ctrl` - `A364.Ctrlse`), xmax = log10(`A364.Ctrl` + `A364.Ctrlse`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5)+#facet_grid(~Chr)+
  scale_color_manual(values = Self1)+
  labs(x='log10(mRNA errors Frequency from our A364A)',y='log10(mRNA errors Frequency from BY4742)',title = paste0("Ours"))+my_theme1

cor.test(Lynch_A364A$Ctrl,Lynch_A364A$`A364.Ctrl`,method = "s")

#######################################

#Ctrl VS  H2O2

#######################################
colnames(Lynch_A364A)
#mRNA
#BY4742
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig5/BY4742_mRNA.pdf",height = 3,width = 3,useDingbats = F)

lm2 <- lmodel2(log10(Lynch_A364A$`H2O2.1h.`)~log10(Lynch_A364A$`Ctrl`))
ggplot(Lynch_A364A, aes(x= log10(`Ctrl`), y = log10(`H2O2.1h.`), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  scale_color_manual(values = Self1)+  
  scale_x_continuous(limits = c(-6,-3.9))+
  labs(x='log10(mRNA errors Frequency from ours)',y='log10(mRNA errors Frequency from Lynch)')+#,title = paste0("BY4742")
  my_theme2+guides(col=F)+
  geom_errorbar(aes(ymin = log10(`H2O2.1h.` - `H2O2.1h.se`), ymax = log10(`H2O2.1h.` + `H2O2.1h.se`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5 )+
  geom_errorbar(aes(xmin = log10(`Ctrl` - `Ctrlse`), xmax = log10(`Ctrl` + `Ctrlse`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5)+
  geom_abline(intercept=as.numeric(lm2$regression.results[3,2]),
              slope=as.numeric(lm2$regression.results[3,3]),color='brown')

dev.off()

#A364A
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig5/A364A_mRNA.pdf",height = 3,width = 3,useDingbats = F)
lm2 <- lmodel2(log10(Lynch_A364A$`A364.H2O2.`)~log10(Lynch_A364A$`A364.Ctrl`))
ggplot(Lynch_A364A, aes(x= log10(`A364.Ctrl`), y = log10(`A364.H2O2.`), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  scale_color_manual(values = Self1)+  
  scale_x_continuous(limits = c(-6,-3.5))+
  labs(x='log10(mRNA errors Frequency from ours)',y='log10(mRNA errors Frequency from Lynch)')+#,title = paste0("BY4742")
  my_theme2+guides(col=F)+
  geom_errorbar(aes(ymin = log10(`A364.H2O2.` - `A364.H2O2.se`), ymax = log10(`A364.H2O2.` + `A364.H2O2.se`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5 )+
  geom_errorbar(aes(xmin = log10(`A364.Ctrl` - `A364.Ctrlse`), xmax = log10(`A364.Ctrl` + `A364.Ctrlse`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5)+
  geom_abline(intercept=as.numeric(lm2$regression.results[3,2]),
             slope=as.numeric(lm2$regression.results[3,3]),color='brown')
dev.off()

#20S RNA virus
test<-read.table("/Dell/Dell13/shankj/projects/Cov/TranscriptionError/Ours/A364A/A364_Freq.txt",header = T,sep="\t")
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210414/Fig5/A364A_20S.pdf",height = 3,width = 3,useDingbats = F)
ssRNA<-filter(test,Chr=="20S")
ssRNA$SNP<-factor(ssRNA$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))
lm2 <- lmodel2(log10(ssRNA$`A364.H2O2.`)~log10(ssRNA$`A364.Ctrl`))
ggplot(ssRNA, aes(x= log10(`A364.Ctrl`), y = log10(`A364.H2O2.`), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  scale_color_manual(values = Self1)+  
  scale_x_continuous(limits = c(-5.5,-3.4))+
  #my_theme2+guides(col=F)+
  geom_errorbar(aes(ymin = log10(`A364.H2O2.` - `A364.H2O2.se`), ymax = log10(`A364.H2O2.` + `A364.H2O2.se`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5 )+
  geom_errorbar(aes(xmin = log10(`A364.Ctrl` - `A364.Ctrlse`), xmax = log10(`A364.Ctrl` + `A364.Ctrlse`)), 
                position = position_dodge(0.5),width=0,col='black',size=.5)+
  geom_abline(intercept=as.numeric(lm2$regression.results[3,2]),
              slope=as.numeric(lm2$regression.results[3,3]),color='brown')
dev.off()

mRNA<-filter(test,Chr=="mRNA")
ssRNA<-filter(test,Chr=="20S")
A364A<-merge(mRNA,ssRNA,by="SNP")
head(A364A)
cor.test(A364A$`A364.Ctrl.y`,A364A$`A364.Ctrl.x`,method = "s")
cor.test(A364A$`A364.H2O2..y`,A364A$`A364.H2O2..x`,method = "s")

A364A$SNP<-factor(A364A$SNP,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

lm2 <- lmodel2(log10(A364A$`A364.Ctrl.y`)~log10(A364A$`A364.Ctrl.x`))
ggplot(A364A, aes(x= log10(`A364.Ctrl.y`), y = log10(`A364.Ctrl.x`), col = SNP)) +
  geom_point(position = position_dodge(0.1),size=3) +
  scale_color_manual(values = Self1)+  
  labs(x='log10(20S RNA virus mutation rate)',y='log10(mRNA errors rate)')+#,title = paste0("BY4742")
  #scale_x_continuous(limits = c(-6,-3.5))+
   my_theme2+guides(col=F)#+
  # geom_errorbar(aes(ymin = log10(`A364.H2O2.` - `A364.H2O2.se`), ymax = log10(`A364.H2O2.` + `A364.H2O2.se`)), 
  #               position = position_dodge(0.5),width=0,col='black',size=.5 )+
  # geom_errorbar(aes(xmin = log10(`A364.Ctrl` - `A364.Ctrlse`), xmax = log10(`A364.Ctrl` + `A364.Ctrlse`)), 
  #               position = position_dodge(0.5),width=0,col='black',size=.5)+
  # geom_abline(intercept=as.numeric(lm2$regression.results[3,2]),
  #             slope=as.numeric(lm2$regression.results[3,3]),color='brown')
