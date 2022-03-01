##SD Monocytes
library(ggplot2)
library(plyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(precrec)

res.col<-c("RD-R"="#E64B35FF", "RD-nR"="#4DBBD5FF", "pCR"="#00A087FF", "No NAC"="#3C5488FF", "RD"="#E64B35FF", "FALSE" = "grey20", "TRUE"= "#F39B7FFF")

theme_set(theme_classic()) 
theme_update(axis.text=element_text(color="black", size=10), axis.title.y = element_text(size=12), axis.title.x = element_text(size=11))


scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

read.csv("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/Manuscript/TidyDataCode/SD_data.csv")->sd.data

sd.rel.mon<-sd.data[sd.data$type %in% c("MONORE", "Monocy"),] #select just relative monocytes
sd.rel.mon<-na.omit(sd.rel.mon)


##Take the average monocyte value for all measured monocytes in 30 day pre-surgery interval
ids<-unique(sd.rel.mon$ID)
result<-NULL
for (i in ids){
  tmp<-sd.rel.mon[sd.rel.mon$ID %in% i,]
  tmp2<-cbind(i, mean(tmp$value))
  result<-rbind(result, tmp2)
}
colnames(result)<-c("ID", "Avg.Mono")
sd.rel.mon2<-sd.rel.mon[!duplicated(sd.rel.mon$ID),]

avg.data<-merge(result, sd.rel.mon2, by.x="ID")
###
avg.data$Avg.Mono<-as.numeric(avg.data$Avg.Mono)

#Fig 2d

stat<-wilcox_test(avg.data, Avg.Mono~response, alternative = "g")%>% add_significance()
d3<-ggplot(avg.data, aes(response, Avg.Mono)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, aes(color=response), size=.1)+
  labs(y="Relative Monocytes", title="VICC-SD", x=element_blank())+
  scale_color_manual(values=res.col)+
  scale_x_discrete(limits=c("RD", "pCR"))+
  stat_pvalue_manual(stat, y.position = 20, label="p.signif")+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=12))

table(avg.data$response)

her2<-avg.data[avg.data$HER2=="TRUE",]
wilcox_test(her2, Avg.Mono~response, alternative = "g")%>% add_significance()

s3d3<-ggplot(her2, aes(response, Avg.Mono)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, aes(color=response), size=.1)+
  labs(y="Relative Monocytes", title="VICC-SD\nHER2+ Only", x=element_blank())+
  scale_color_manual(values=res.col)+
  scale_x_discrete(limits=c("RD", "pCR"))+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=12))

##ROC for whole SD cohort
precrec_obj <- evalmod(scores = avg.data$Avg.Mono, labels = avg.data$response)
auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.5, y=0.25, label="AUC: 0.61")+
  theme(axis.text = element_text(color="black"), axis.title = element_text(size=14), title = element_blank())


###Follow up data added###
tnbc2<-read.csv("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Monocytes/SD/tnbc_patients.csv")

tnbc2<-tnbc2[!tnbc2$Recur %in% c("InsufficientFU"),]
tnbc2$Recur<-revalue(tnbc2$Recur, c("No"= "RD-nR", "Yes" = "RD-R"))

table(tnbc2$Recur)


#### Take the average monocyte value for all measured monocytes in 30 day pre-surgery interval
## for TNBC patients with follow-up
ids<-unique(tnbc2$ID)
result<-NULL
for (i in ids){
  tmp<-tnbc2[tnbc2$ID %in% i,]
  tmp2<-cbind(i, mean(tmp$value))
  result<-rbind(result, tmp2)
}
colnames(result)<-c("ID", "Avg.Mono")
tnbc3<-tnbc2[!duplicated(tnbc2$ID),]
avg.data<-merge(result, tnbc3, by.x="ID")
avg.data$Avg.Mono<-as.numeric(avg.data$Avg.Mono)


#Supplementary Fig 2H
stat<-wilcox_test(avg.data, Avg.Mono~Recur, alternative = "g") %>% add_significance()
stat<-stat[stat$p.adj<0.05,]

s3d2<-ggplot(avg.data, aes(Recur, Avg.Mono)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, aes(color=Recur), size=.1)+
  labs(y="Relative Monocytes", title="VICC-SD\nTNBC Only", x=element_blank())+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR"))+
  scale_color_manual(values=res.col)+
  stat_pvalue_manual(stat, y.position = c(17), label="p.adj.signif")+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


#Fig 2D
stat<-wilcox_test(avg.data, Avg.Mono~response, alternative = "g") %>% add_significance()

d4<-ggplot(avg.data, aes(response, Avg.Mono)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, aes(color=Recur), size=.1)+
  labs(y="Relative Monocytes", title="VICC-SD\nTNBC Only", x=element_blank())+
  scale_x_discrete(limits=c("RD","pCR"))+
  scale_color_manual(values=res.col)+
  stat_pvalue_manual(stat, y.position = c(16.5), label="p.signif")+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5, size=12))


table(avg.data$Recur)

d1|d2|d3|d4

#ROC ##TNBC only
precrec_obj <- evalmod(scores = avg.data$Avg.Mono, labels = avg.data$response)
auc(precrec_obj)
d5<-autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.7, y=0.05, label="AUC: 0.71")+
  labs(title="VICC-SD\nTNBC only")+
  theme(axis.text = element_text(color="black", size=8), axis.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=12))



avg.data$response2<-avg.data$Recur
avg.data$response2<-revalue(avg.data$Recur, c("pCR"="nR", "RD-nR"="nR", "RD-R"="R"))
precrec_obj <- evalmod(scores = avg.data$Avg.Mono, labels = avg.data$response2)
auc(precrec_obj)
autoplot(precrec_obj, "ROC")




###HR+ Only
hr2<-read.csv("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Monocytes/SD/hr_patients.csv") #read in data with follow-up
hr2<-hr2[!hr2$Recur %in% c("InsufficientFU"),]
hr2$Recur<-revalue(hr2$Recur, c("No"= "RD-nR", "Yes" = "RD-R"))

##averages hr
ids<-unique(hr2$ID)
result<-NULL
for (i in ids){
  tmp<-hr2[hr2$ID %in% i,]
  tmp2<-cbind(i, mean(tmp$value))
  result<-rbind(result, tmp2)
}
colnames(result)<-c("ID", "Avg.Mono")

hr3<-hr2[!duplicated(hr2$ID),]
avg.data<-merge(result, hr3, by.x="ID")
avg.data$Avg.Mono<-as.numeric(avg.data$Avg.Mono)

wilcox_test(avg.data, Avg.Mono~Recur, alternative = "g")
wilcox_test(avg.data, Avg.Mono~response, alternative = "g")

#Supplementary Fig 2H
s3d1<-ggplot(avg.data, aes(response, Avg.Mono)) +geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, aes(color=Recur), size=.1)+
  labs(y="Relative Monocytes", title="VICC-SD\nHR+ Only", x=element_blank())+
  scale_x_discrete(limits=c("RD", "pCR"))+
  scale_color_manual(values=res.col)+
  theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=12))



############Additional Cohorts############

#Supplementary Fig 2G
ivdata<-read.csv("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/Manuscript/TidyDataCode/InstitutoValencianoMono.csv")

stat<-wilcox_test(ivdata, Avg.Mono ~ Metastasis, alternative = "l") %>% add_significance()
stat$pe<-paste0("p=", stat$p)
s3c<-ggplot(ivdata, aes(as.factor(Metastasis), Avg.Mono))+ geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.1, size=1)+labs(x="Metastasis", y= "Relative Monocytes", title="Instituto\nValenciano")+
  scale_x_discrete(labels=c("Yes", "No"))+
  ylim(0,14)+
  stat_pvalue_manual(stat, label="pe", y.position = 13.31)+
  theme(plot.title = element_text(hjust = 0.5, size=12))


precrec_obj <- evalmod(scores = ivdata$Avg.Mono, labels = ivdata$Metastasis)
aucs <- auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.5, y=0.25, label="AUC: 0.73")+
  theme(axis.text = element_text(color="black"), axis.title = element_text(size=14), title = element_blank())



############### Gepar Nuevo data

setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/Monocytes/JITC/")
read.csv("data.csv")->data
data$ID<-gsub("T4", "", data$proben)
data$ID<-gsub("t4", "", data$ID)
data$ID<-as.numeric(data$ID)
data<-data[!duplicated(data$ID),]

read.delim("pCR.txt")->outcome

full<-merge(data, outcome, by="ID")
full$MONO_0<-as.numeric(full$MONO_0)

anyGCSF<-c("228","94", "115", "189", "130", "136", "142", "192", "84", "103", "108") 

noGCSF<-full[!full$ID %in% anyGCSF,]


#no GCSF
stat1<-wilcox_test(noGCSF, formula= MONO_1 ~ ypT0_ypN0, alternative = "l")
stat1$lab<- paste0("p=", stat1$p)

s3a1<-ggplot(noGCSF, aes(ypT0_ypN0, MONO_1))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, size=.1)+
  labs(x=element_blank(), y= "Relative Monocytes" , title="GeparNuevo")+
  scale_x_discrete(limits=c("no pCR", "pCR"), labels=c("RD", "pCR"))+
  ylim(0,31)+
  stat_pvalue_manual(stat1, label= "lab", y.position = max(noGCSF$MONO_1)+1.5)+
  theme(plot.title = element_text(hjust = 0.5, size=12))


stat0<-wilcox_test(noGCSF, formula= MONO_0 ~ ypT0_ypN0, alternative = "l")
stat0$lab<-paste0("p=", stat0$p)

s3a2<-ggplot(noGCSF, aes(ypT0_ypN0, MONO_0))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, size=.1)+
  ylim(0,1.6)+
  labs(x=element_blank(), y= "Absolute Monocytes", title="GeparNuevo")+
  scale_x_discrete(limits=c("no pCR", "pCR"), labels=c("RD", "pCR"))+
  stat_pvalue_manual(stat0, label= "lab", y.position = max(noGCSF$MONO_0)+.1)+
  theme(plot.title = element_text(hjust = 0.5, size=12))

table(noGCSF$ypT0_ypN0)

precrec_obj <- evalmod(scores = noGCSF$MONO_0, labels = noGCSF$ypT0_ypN0)
aucs <- auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.5, y=0.25, label="AUC: 0.64")+
  theme(axis.text = element_text(color="black"), axis.title = element_text(size=14), title = element_blank())

precrec_obj <- evalmod(scores = noGCSF$MONO_1, labels = noGCSF$ypT0_ypN0)
aucs <- auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.5, y=0.25, label="AUC: 0.58")+
  theme(axis.title = element_text(size=14), title = element_blank())



####PIRS nanostring


setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/JITC_RNAseq/NanoString_PIRS/")


read.delim("metadata.txt", row.names=1)->meta
read.delim('normalized.txt', row.names = 1)->data
rownames(meta)<-paste0("X", rownames(meta))
all(rownames(meta)==colnames(data))
log(data+1,2)->data



cty.genes<-c("FGFBP2", "GNLY", "GZMB", "GZMH", "LAG3", "PDCD1", "NKG7", "HLA-G") 
ic.genes<-read.delim("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/IFNCompGenes.txt", row.names=1, check.names = FALSE)


cty<-data[rownames(data) %in% cty.genes,]
ic<-data[rownames(data) %in% ic.genes$Gene,]


t1<-as.data.frame(t(scaleRow(cty)))
t3<-as.data.frame(t(scaleRow(ic)))
all(rownames(t1)==rownames(meta))

t2<-t1[,colnames(t1)%in% c("FGFBP2", "GNLY", "GZMB", "GZMH", "LAG3", "PDCD1", "NKG7")]

meta$cyt.score<-(rowSums(t2)-t1$`HLA-G` )/ ncol(t1)
meta$ic.score<-rowSums(t3)/ ncol(t3)


##PIRS
meta$PIRS<-meta$ic.score-meta$cyt.score

meta<-meta[!rownames(meta) %in% c("X24"),] #exclude for gcsf


t_test(meta, PIRS~Response, alternative = "g")

s3b<-ggplot(meta, aes(Response, PIRS))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(height=0, width=0.15, size=0.75)+
  labs(y="PIRS", title="GeparNeuvo", x=element_blank())+
  scale_x_discrete(labels=c("RD", "pCR"))+
  theme(legend.position="none", axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=11), plot.title = element_text(hjust = 0.5, size=12),
        axis.text = element_text(size=10, color="black"))


precrec_obj <- evalmod(scores = meta$PIRS, labels = meta$Response)
auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.5, y=0.25, label="AUC: 0.41")+
  theme(axis.text = element_text(color="black"), axis.title = element_text(size=14), title = element_blank())
