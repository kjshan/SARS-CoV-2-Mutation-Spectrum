#########################################################################

#Fig6

#########################################################################
Count<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210711_Test6/SCV2/Delet/Delete_SNP.countNew",sep="\t",header = T,stringsAsFactors = F)
head(Count)
#colnames(Count)<-c("Branch","SNP","count","G","C","Species")
Count$SNP<-factor(Count$SNP,levels = c("CT","GA","AG","TC","GT","CA","GC","CG","AC","TG","AT","TA"))
Count<-filter(Count,Species=="SCV2",!Branch_Name %in% c("N5:Rc_o319_LC556375.1","N2:BatBM48_31_NC_014470.1","SCV2_Among_Human"))
head(Count)
unique(Count$Branch)


filter(Count,Branch=="N3:Human_Astr")
sum(as.numeric(filter(Count,Branch=="N3:Human_Astr")$count))

plot_list <- list()
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/Fig6S/FigSB-1.pdf",height = 18,width = 18,useDingbats = F)
for (i in unique(Count$Branch)) {
  #i="N5:GX_P5L_EPI_ISL_410540"
  mid<-filter(Count,Branch==i)
  sum1<-sum(mid$Count)
  mid$rate<-mid$Count/sum1
  mid$SNP<-factor(mid$SNP,levels = c("CT","GA","AG","TC","GT","CA","GC","CG","AC","TG","AT","TA"))
  plot_list[[i]] <-ggplot(data=mid, aes(x=factor(SNP),y=rate)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_y_continuous(breaks = seq(0,0.4,0.2),limits = c(0,0.45)) +
    #scale_fill_manual(values = Self)+#facet_grid(~Branch)+
    labs(x='',y='SNV count',title = i)+my_theme#2+guides(fill=F)
}
do.call(grid.arrange, plot_list)


dev.off()

for (i in unique(Count$Branch)) {
  pdf(paste("/Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210711_Test6/",i,".pdf",sep=""),height = 12,width = 12,useDingbats = F)
    print(plot_list[[i]])
  dev.off()
}

#Fig6B
test<-read.table("/Dell/Dell13/shankj/projects/Cov/Genome/Fig6/20210711_Test6/SCV2/Delet/MutationSpectrumSummary.txt",header = T)
head(test)
Fig6B<-filter(test,Branch %in% c(paste0("B",0:9),"pSCV2","SCV2_SNP_Vero","Lung"))
mid<-filter(reshape2::dcast(Fig6B,SNP~Branch,value.var = "Count_Rate"))
head(mid)

Fig6D<-data.frame()
Fig6D_mid
i="B6"
for (i in unique(Fig6B$Branch)) {
  for (j in unique(Fig6B$Branch)) {
    mid1<-filter(Fig6B,Branch==i)
    mid2<-filter(Fig6B,Branch==j)
    mid1$Count_Rate[mid1$Count_Rate==0] <- 0.5
    mid2$Count_Rate[mid2$Count_Rate==0] <- 0.5
    mid1$rate<-mid1$Count_Rate/sum(mid1$Count_Rate)
    mid2$rate<-mid2$Count_Rate/sum(mid2$Count_Rate)
    
    merge<-merge(mid1,mid2,by="SNP")
    test<-cor.test(log10(merge$`rate.x`),log10(merge$`rate.y`))
    Fig6D_mid<-data.frame(Name1=i,Name2=j,Correlation=test$estimate,P=test$p.value)
    Fig6D<-rbind(Fig6D_mid,Fig6D)
  }
}

head(Fig6D)



M <- spread(Fig6D[,1:3], Name2, Correlation)
#nrow(Fig6D)
head(M)
rownames(M)<-M$Name2
#corrplot(as.matrix(M[,-1]),order = "hclust",addrect = 3,cl.lim = c(0.3, 1), col = col2(8), tl.cex = 0.7, is.corr = FALSE )
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/Heatmap.pdf",height = 3,width = 3,useDingbats = F)

pheatmap::pheatmap(as.matrix(M[,-1]),cluster_rows = T,show_rownames = T,legend_breaks=seq(0,1,0.2))
#,legend_breaks=30
dev.off()




pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210716/Fig6B_B0B1.pdf",height = 3,width =3,useDingbats = F)
ggplot(data=mid, aes(x=log10(B0/sumB0),y=log10(B1/sumB1),col=SNP)) + 
    geom_point(size=3)+
    scale_color_manual(values = Self1)+
    #scale_y_continuous(breaks = seq(0,0.4,0.2),limits = c(0,0.45)) +
    scale_x_continuous(limits = c(-2.4,-0.4)) +
    scale_y_continuous(limits = c(-2.4,-0.4)) +
    my_theme2+guides(col=F)+
    labs(x='B0 log10',y='B1 log10')
dev.off()
cor.test(log10(mid$B0/sumB0),log10(mid$B1/sumB1))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210716/Fig6B_B0pSCV2.pdf",height = 3,width =3,useDingbats = F)
ggplot(data=mid, aes(x=log10(B0/sumB0),y=log10(B1/sumB1),col=SNP)) + 
  geom_point(size=3)+
  scale_color_manual(values = Self1)+
  #scale_y_continuous(breaks = seq(0,0.4,0.2),limits = c(0,0.45)) +
  #scale_x_continuous(limits = c(-2.4,-0.4)) +
  #scale_y_continuous(limits = c(-2.4,-0.4)) +
  my_theme2+guides(col=F)+
  labs(x='B0 log10',y='pSCV2 log10')
dev.off()

cor.test(log10(mid$B0/sumB0),log10(mid$pSCV2/sumpSCV2))

###############################################################################################

#Fig7

###############################################################################################

#Fig7A
Fig7<-read.table("/Dell/Dell13/shankj/projects/Cov/Plot/Script/Result/Fig7.txt",header = T,sep="\t")
#Fig7<-filter(Fig7,Branch %in% paste0("B",10:25))

Fig7$SNP<-factor(Fig7$SNP,levels = c("CT","GA","AG","TC","GT","CA","GC","CG","AC","TG","AT","TA"))
#Filter$Branch<-factor(Filter$Branch,levels = paste0("B",0:25))
head(Fig7)

plot_list <- list()
pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210715/FigSB-1.pdf",height = 18,width = 15,useDingbats = F)
for (i in paste0("B",10:25)) {
  #i="N5:GX_P5L_EPI_ISL_410540"
  mid<-filter(Fig7,Branch==i)
  sum1<-sum(mid$Count)
  mid$rate<-mid$Count/sum1
  mid$SNP<-factor(mid$SNP,levels = c("CT","GA","AG","TC","GT","CA","GC","CG","AC","TG","AT","TA"))
  plot_list[[i]] <-ggplot(data=mid, aes(x=factor(SNP),y=rate)) + 
    geom_bar(position="dodge", stat="identity",color="black",width = 0.8)+
    scale_y_continuous(breaks = seq(0,0.4,0.2),limits = c(0,0.45)) +
    #scale_fill_manual(values = Self)+#facet_grid(~Branch)+
    labs(x='',y='SNV count')+my_theme2+guides(fill=F)
}
do.call(grid.arrange, plot_list[paste0("B",10:25)])
dev.off()


################################################################################
#Fig7B
head(Fig7)
mid<-reshape2::dcast(Fig7[,1:3],Branch~SNP,value.var = "Count")

mid[,-1]<- mid[,-1]/rowSums(mid[,-1])
row.names(mid) <- mid$Branch
nrow(mid)
sum(mid[1,-1])

pca <- prcomp(mid[,2:13],scale = T)
head(pca)
summary(pca)
aa <- data.frame(Branch=rownames(mid),pca$x)
head(Fig7)

name<-unique(Fig7[,3:4])
a<-merge(aa,name,by="Branch")
unique(a$Species)  
human #4EA0EA
human+bat #2564BF
camel #9966FF
bat #548235
bat+other #99DD5E
Ne #00553A
Ve #04845A
Er #04BC7F
Mix#08FFAD
Self3<-c("#548235",
         "grey",
         "grey",
         "grey",
         "#9966FF",
         "#4EA0EA",
         "#2564BF",
         "#F90D2F")

nrow(Fig7)/12

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210812/AdditionalFigure/Fig7B.pdf",height = 5,width = 5,useDingbats = F)
ggplot()+
  geom_point(data=filter(a,!Branch %in% c("B21","B22","B23","B24","B25")), mapping=aes(x=PC1,y=PC2,color=Species))+
  #geom_point(data=filter(a,grepl("Human",Species)), mapping=aes(x=PC1,y=PC2),color="blue")+ 
  #geom_point(filter(a,Species=="Bat",!Branch %in% c("B20","B21","B22","B23","B24","B25")), mapping =aes(PC1,PC2),color="green")+  #,color=Species
  #geom_point(data=filter(a,Species=="Human"), 
  #           mapping=aes(x=PC1,y=PC2),color="red")+
  # geom_point(data=filter(a,Species=="Bat"), mapping=aes(x=PC1,y=PC2),color="purple")+ 
  ggrepel::geom_text_repel(data =a, aes(x=PC1,y=PC2,label = Branch),color="#000000")+
  stat_ellipse(filter(a,Species=="Bat",!Branch %in% c("B20","B21","B22","B23","B24","B25")), mapping =aes(PC1,PC2),color="orange")+
  stat_ellipse(filter(a,grepl("Human",Species)), mapping =aes(PC1,PC2),color="blue")+
  stat_ellipse(filter(a,Species=="Bat"), mapping =aes(PC1,PC2),color="green")+
  geom_point(data=filter(a,Species=="Bat",!Branch %in% c("B20","B21","B22","B23","B24","B25")), mapping=aes(x=PC1,y=PC2,color=Species))+
  geom_point(data=filter(a,Branch=="B21"), mapping=aes(x=PC1,y=PC2),color="#00553A",shape=3,size=3)+
  geom_point(data=filter(a,Branch=="B23"), mapping=aes(x=PC1,y=PC2),color="#04845A",shape=4,size=3)+
  geom_point(data=filter(a,Branch=="B25"), mapping=aes(x=PC1,y=PC2),color="grey")+
  geom_point(data=filter(a,Branch%in% c("B22","B24")), mapping=aes(x=PC1,y=PC2),color="#08FFAD",shape=2,size=3)+
  scale_x_continuous(limits = c(-8.5,5)) +
  scale_y_continuous(limits = c(-4,4))+
  my_theme2+guides(colour=F)+
  scale_color_manual(values = Self3)#+my_theme1

dev.off()

################################################################################
#Fig6B
#B0
Fig7<-read.table("/Dell/Dell13/shankj/projects/Cov/SomaticMutations/Fig7.txt",header = T,sep="\t")
SampleB0_1000<-data.frame()
B0_B1<-merge(B0[,1:2],B1[,1:2],by="SNP")
colnames(B0_B1)<-c("SNP","B0","B1")
cor.test(log10(B0_B1$B0/529),log10(B0_B1$B1/628))

for (i in 1:1000) {
B0<-filter(Fig7,Branch=="B0")
B0_SNP<-c(rep(as.character(B0$SNP),B0$Count))
B0_rep<-sum(B0$Count)
B0_indexes <- sample (B0_SNP, B0_rep, replace =T) 
SampleB0<-as.data.frame(t(table(B0_indexes)))[,2:3]
SampleB0$rep<-rep(i,nrow(SampleB0))
SampleB0_1000<-rbind(SampleB0,SampleB0_1000)
}
colnames(SampleB0_1000)<-c("SNP","B0","Sample")

SampleB1_1000<-data.frame()
for (i in 1:1000) {
  B1<-filter(Fig7,Branch=="B1")
  B1_SNP<-c(rep(as.character(B1$SNP),B1$Count))
  B1_rep<-sum(B1$Count)
  B1_indexes <- sample (B1_SNP, B1_rep, replace =T) 
  SampleB1<-as.data.frame(t(table(B1_indexes)))[,2:3]
  SampleB1$rep<-rep(i,nrow(SampleB1))
  SampleB1_1000<-rbind(SampleB1,SampleB1_1000)
  
}
colnames(SampleB1_1000)<-c("SNP","B1","Sample")

SamplepSCV2_1000<-data.frame()
for (i in 1:1000) {
  pSCV2<-filter(Fig7,Branch=="pSCV2")
  pSCV2_SNP<-c(rep(as.character(pSCV2$SNP),pSCV2$Count))
  pSCV2_rep<-sum(pSCV2$Count)
  pSCV2_indexes <- sample (pSCV2_SNP, pSCV2_rep, replace =T) 
  SamplepSCV2<-as.data.frame(t(table(pSCV2_indexes)))[,2:3]
  SamplepSCV2$rep<-rep(i,nrow(SamplepSCV2))
  SamplepSCV2_1000<-rbind(SamplepSCV2,SamplepSCV2_1000)
}
colnames(SamplepSCV2_1000)<-c("SNP","pSCV2","Sample")
setwd("/Dell/Dell13/shankj/projects/Cov/Plot/Response/")
B0_B1<-merge(SampleB0_1000,SampleB1_1000,by=c("SNP","Sample"), all=TRUE)
B0_B1_pSCV2<-merge(B0_B1,SamplepSCV2_1000,by=c("SNP","Sample"), all=TRUE)
B0_B1_pSCV2[is.na(B0_B1_pSCV2)] <- 0
write.csv(B0_B1_pSCV2,"/Dell/Dell13/shankj/projects/Cov/Plot/Response/B0_B1_pSCV2.csv",row.names = F,quote = F)

head(B0_B1_pSCV2,12)

Correlation<-data.frame()
for (i in 1:1000) {
  mid<-filter(B0_B1_pSCV2,Sample==i)
  sumB0<-sum(mid$B0+0.5)
  sumB1<-sum(mid$B1+0.5)
  sumpSCV2<-sum(mid$pSCV2) 
  

  B1<-cor.test(log10((mid$B0+0.5)/sumB0),log10((mid$B1+0.5)/sumB1))
  pSCV2<-cor.test(log10((mid$B0+0.5)/sumB0),log10(mid$pSCV2/sumpSCV2))
  
  Cor_mid<-data.frame(Sample=i,B0VSB1_Cor=B1$estimate,B0VSB1_P=B1$p.value,B0VSpSCV2_Cor=pSCV2$estimate,B0VSpSCV2_P=B1$p.value)
  Correlation<-rbind(Cor_mid,Correlation)
}
Correlation<-as.data.frame(Correlation)
head()

ggplot(mid)+
  labs(x='Correlation')+
  geom_density(aes(x=B0VSB1_Cor),col="red")+
  geom_vline(xintercept = 0.98,color="red",linetype="dotted",size=1)+
  geom_density(aes(x=B0VSpSCV2_Cor),col="blue")+
  geom_vline(xintercept = 0.69,color="blue",linetype="dotted",size=1)+
  my_theme1

summarySE(Correlation,measurevar="B0VSB1_Cor")
summarySE(Correlation,measurevar="B0VSpSCV2_Cor")
mid<-na.omit(Correlation)
Correlation<-data.frame()
for (i in 1:1000) {
  mid<-filter(B0_B1_pSCV2,Sample==i)
  sumB0<-sum(mid$B0)
  sumB1<-sum(mid$B1)
  sumpSCV2<-sum(mid$pSCV2) 
  
  
  B1<-cor.test(log10((mid$B0)/sumB0),log10((mid$B1)/sumB1))
  pSCV2<-cor.test(log10((mid$B0)/sumB0),log10(mid$pSCV2/sumpSCV2))
  
  Cor_mid<-data.frame(Sample=i,B0VSB1_Cor=B1$estimate,B0VSB1_P=B1$p.value,B0VSpSCV2_Cor=pSCV2$estimate,B0VSpSCV2_P=B1$p.value)
  Correlation<-rbind(Cor_mid,Correlation)
}
Correlation<-as.data.frame(Correlation)

