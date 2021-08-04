#####################################################################################################

#We discarded mutations detected in multiple humans to reduce the interference from the potential standing polymorphisms in the population. 
#We also discarded somatic mutations that existed in multiple tissues of the same human to exclude potential RNA editing events.

#####################################################################################################
somatic<-read.table("/Dell/Dell13/shankj/projects/Cov/SomaticMutations/Somatic.txt",header = F,sep = "\t")
mid<-somatic %>% group_by(V2,V3,V5,V6)   %>% dplyr::summarise(count=n()) 
test<-merge(mid,somatic,by=c("V2","V3","V5"))
head(test)
head(mid)
nrow(somatic)
nrow(filter(test,count==1))
count<-filter(test,count==1) %>% group_by(V5,V6)  %>% dplyr::summarise(count=n())
head(count)

#####################################################################################################

#somatic mutation spectrum 

#####################################################################################################
count$V5<-factor(count$V5,levels = c("C > T","G > A","A > G","T > C","G > T","C > A","G > C","C > G","A > C","T > G","A > T","T > A"))

pdf("/Dell/Dell13/shankj/projects/Cov/Plot/20210712/20210715/Organ-1.pdf",height = 18,width = 18,useDingbats = F)

plot_list=list()
for (i in unique(count$V6)) {
  med<-filter(count,V6==i)
  med$rate<-med$count/sum(med$count)
  sum<-sum(med$count)
  plot_list[[i]]<-ggplot(data=med, aes(x=V5,y=rate)) + #geom_abline(intercept=0.5, slope=0,color="red")+
    geom_bar(position="dodge", stat="identity",width = 0.8)+
    scale_y_continuous(breaks = seq(0,0.4,0.2),limits = c(0,0.42))+ 
    #scale_fill_manual(values = Self)+#scale_y_continuous(breaks = seq(0,2000,500),limits = c(0,2000))+
    labs(x='SNP',y='rate',title = paste0(str_replace_all(string=i,pattern="_",replacement = " ")," N=",sum))+#+#facet_grid(~`Sample.Type`)
    my_theme#2+guides(fill=F)
  
}
do.call(grid.arrange, plot_list)
dev.off()

#####################################################################################################

#Fig5

#####################################################################################################
head(count)
odds<-data.frame()
for (i in unique(count$V6)) {
  Fisher1<-fisher.test(matrix(c(filter(count,V6==i,V5=="G > T")$count,
                                filter(count,V6==i,V5=="C > A")$count,
                                (278064958-filter(count,V6==i,V5=="G > T")$count),
                                (277765994-filter(count,V6==i,V5=="C > A")$count)),nrow=2))
  Fisher2<-fisher.test(matrix(c(filter(count,V6==i,V5=="C > T")$count,
                                filter(count,V6==i,V5=="G > A")$count,
                                (277765994-filter(count,V6==i,V5=="C > T")$count),
                                (278064958-filter(count,V6==i,V5=="G > A")$count)),nrow=2))
  
  mid1<-data.frame(sample=i,P1=Fisher1$p.value, odds1=Fisher1$estimate,P2=Fisher2$p.value, odds2=Fisher2$estimate)
  odds<-rbind(odds,mid1)
}

head(odds)
Polymorphism<-data.frame(sample="SARS-CoV-2 polymorphism",P1=3.631382e-188, odds1=6.2237,P2=1.033515e-199, odds2=3.797579)


#fisher.test(matrix(c(51,15,2354397,2185628),nrow = 2))
#fisher.test(matrix(c(6,24,2354397,2185628),nrow = 2))

odds<-rbind(odds,Polymorphism)
tail(odds)


odds$q1<-qvalue(odds$P1,pi0 = 1)$qvalues
odds$q2<-qvalue(odds$P2,pi0 = 1)$qvalues
odds$seq<-1:nrow(odds)
ggplot() + 
  geom_point(data=filter(odds,q2>0.01,q1>0.01), aes(x=log2(odds2),y=log2(odds1)),color="grey")+
  geom_point(data=filter(odds,q1<=0.01), aes(x=log2(odds2),y=log2(odds1)),color="#386CB0")+ 
  geom_point(data=filter(odds,q2<=0.01), aes(x=log2(odds2),y=log2(odds1)),color="#E6AB01")+ 
  geom_point(data=filter(odds,q2<=0.01,q1<=0.01), aes(x=log2(odds2),y=log2(odds1)),color="#F02A7F")+ 
  labs(x='log2(odds ratio C>T/G>A)',y='log2(odds ratio G>T/C>A)')+
  scale_x_continuous(limits = c(-1,4)) +
  scale_y_continuous(limits = c(-1,4)) +
  geom_hline(yintercept=0,col="black",linetype="dashed")+
  geom_vline(xintercept=0,col="black",linetype="dashed")+
  ggrepel::geom_text_repel(data =odds, aes(x=log2(odds2),y=log2(odds1),label = sample),col="black",segment.size = 1,min.segment.length = 1)+
  my_theme2+guides(col=F)



