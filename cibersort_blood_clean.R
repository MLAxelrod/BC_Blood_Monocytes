##CIBERSORT blood
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(plyr)
library(ggfortify)
library(corrplot)
library(colorspace)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(precrec)

theme_set(theme_classic()) 
theme_update(axis.text=element_text(color="black", size=10))


res.col<-c("RD-R"="#E64B35FF", "RD-nR"="#4DBBD5FF", "pCR"="#00A087FF", "No NAC"="#3C5488FF")

scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

### CIBERSORT ####
#Load and format data for metadata:data matching
#Job 12 absolute mode with 500 permutations
#Job 13 relative mode with 500 permutations

setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/Cibersort/")
read.delim("CIBERSORTx_Job12_Results.txt", row.names=1)->abs
read.delim("CIBERSORTx_Job13_Results.txt", row.names=1)->rel
all(rownames(abs)==rownames(rel))
setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/")
read.delim("metadata.txt", row.names=1)->metadata

metadata$RR<-revalue(metadata$RR, c("RD NO" = "RD-nR", "RD YES" = "RD-R"))

#Match data to metadata
order<-rownames(abs)
foo<-sub("(.*?_.*?_.*)_.*", "\\1", order)
meta<-metadata[match(foo, rownames(metadata)),]
all(rownames(meta)==foo)
rownames(abs)<-foo
rownames(rel)<-foo
all(rownames(meta)==rownames(abs))

cib.abs<-cbind(abs, meta)
cib.rel<-cbind(rel, meta)

#### HEATMAP of relative cibersort data ###
#simplify cib cats
simp<-rel
simp$T.cells.CD4<-simp$T.cells.CD4.memory.activated + simp$T.cells.CD4.memory.resting +simp$T.cells.CD4.naive+simp$T.cells.regulatory..Tregs.
simp$B.cells<-simp$B.cells.naive+simp$B.cells.memory+simp$Plasma.cells
simp$NK.cells<-simp$NK.cells.activated+simp$NK.cells.resting
simp$Other.Myeloid<-simp$Macrophages.M0+simp$Macrophages.M1+simp$Macrophages.M2+simp$Dendritic.cells.activated+simp$Dendritic.cells.resting
simp$Granulocytes<-simp$Mast.cells.activated+simp$Mast.cells.resting+simp$Eosinophils+simp$Neutrophils

cols.to.keep<-c("Monocytes", "T.cells.CD4", "T.cells.CD8","B.cells", "NK.cells", "Other.Myeloid", "Granulocytes")
simp2<-simp[,colnames(simp) %in% cols.to.keep]

h<-as.matrix(simp2)
meta$rownam<-rownames(meta)
meta.sort<-meta[order(meta$TNBC),]
RR.order<-c("pCR", "RD-nR", "RD-R", "No NAC")
meta.sort<-meta.sort[order(match(meta.sort$RR, RR.order)),]

h<-h[match(rownames(meta.sort), rownames(h)),]

data.frame(meta.sort[["RR"]], meta.sort[["TNBC"]])->annData
colnames(annData)<-c("Response", "TNBC")
f<-HeatmapAnnotation(df= annData, 
                     col = list(Response=c("No NAC" = "#3C5488FF","RD-R" = "#E64B35FF",
                                           "RD-nR" = "#4DBBD5FF", "pCR" ="#00A087FF"),
                                TNBC= c("NO" = "grey20", "YES"= "#F39B7FFF")))


t<-t(h)
t2<-t[order(-rowSums(t)),] #sort to put most abundant on top
t3<-scaleRow(t2)
rownames(t3)
rownames(t3)<-c("Monocytes", "CD4 T cells", "B cells", "NK cells", "CD8 T cells", "Other Myeloid", "Granulocytes")

##Figure 2a
Heatmap(t3, top_annotation = f, row_names_gp = gpar(fontsize=12),
        cluster_columns = FALSE, cluster_rows = FALSE, show_column_names = FALSE,
        column_split = factor(annData$Response, levels=c("No NAC", "RD-R", "RD-nR", "pCR")),
        cluster_column_slices = FALSE, column_gap= unit(c(3,1,1), "mm"),
        heatmap_legend_param = list(title="Row Z scores", title_position= "lefttop-rot"),
        border=TRUE)


### Plots of cell types by response ##

cibplt<-function(cells, ylab){
  ggplot(cib.rel, aes(Response, cells))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+
    scale_x_discrete(limits=c("No NAC", "RD", "pCR"),
                     labels=c("No\nNAC", "RD", "pCR"))+
    labs(color= "Response", x=element_blank(), y=ylab)+
    scale_color_manual(values=res.col)+
    theme(axis.title.y=element_text(size=12), 
          axis.text.x = element_text(size=10), legend.position = "none")
}

rrplt<-function(cells, ylab){
  ggplot(cib.rel, aes(RR, cells))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+  
    scale_x_discrete(limits=c("No NAC", "RD-R", "RD-nR", "pCR"))+
    labs(y=ylab, x=element_blank())+
    scale_color_manual(values=res.col)+
    theme(legend.position="none", axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=10))
}


##Granulocytes
cib.rel$Granulocytes<-cib.rel$Mast.cells.activated+cib.rel$Mast.cells.resting+cib.rel$Eosinophils+cib.rel$Neutrophils

wilcox_test(data=cib.rel, formula=Granulocytes~Response, p.adjust.method = "fdr")
cibplt(cib.rel$Granulocytes, "Granulocytes")
wilcox_test(data=cib.rel, formula=Granulocytes~RR, p.adjust.method = "fdr")
rrplt(cib.rel$Granulocytes, "Granulocytes (CIBERSORTx)") +   geom_bracket(xmin=c("No NAC"), xmax=c("pCR"), y.position = 0.032, label="*")




#Naive B Cells
stat<-wilcox_test(data=cib.rel, formula=B.cells.naive~Response, p.adjust.method = "fdr")
s1<-cibplt(cib.rel$B.cells.naive, "Naive B Cells (CIBERSORTx)")+
  stat_pvalue_manual(stat, label="p.adj.signif", y.position = c(0.39,0.36,0.31))

#Naive B cells by RR
stat<-wilcox_test(data=cib.rel, formula=B.cells.naive~RR, p.adjust.method = "fdr")
stat<-stat[stat$p.adj<0.05,]
s3<-rrplt(cib.rel$B.cells.naive, "Naive B Cells (CIBERSORTx)")+  
  stat_pvalue_manual(stat, label="p.adj.signif", y.position = c(0.39,0.37,0.35, 0.33))+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

##Figure 2b
#CIB PLOT MONOCYTES
stat<-wilcox_test(data=cib.rel, formula=Monocytes~Response, p.adjust.method = "fdr")
p2b<-cibplt(cib.rel$Monocytes, "Relative Monocytes")+
  stat_pvalue_manual(stat, label="p.adj.signif", y.position = c(0.92,0.8,0.88))+
  labs(title="CIBERSORTx")+ theme(plot.title = element_text(hjust = 0.5, size=12))



#Monocytes by RR
stat<-wilcox_test(data=cib.rel, formula=Monocytes~RR,  p.adjust.method = "fdr")
stat<-stat[stat$p.adj<0.05,]
s2<-rrplt(cib.rel$Monocytes, "Monocytes (CIBERSORTx)") +
  stat_pvalue_manual(stat, label="p.adj.signif", y.position = c(0.87, 0.8))+
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

#supplementary figure 2ab
s1|s2|s3


###do cib plt monocytes on just TNBC or HR or HER2
f1<-ggplot(subset(cib.rel, TNBC=="YES"), aes(Response, Monocytes))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+
  scale_x_discrete(limits=c("No NAC", "RD", "pCR"), labels=c("No\nNAC", "RD", "pCR"))+
  labs(x=element_blank(), y ="Relative Monocytes",
       title="TNBC only")+
  scale_color_manual(values=res.col)+
  theme(axis.title.y=element_text(size=12),legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=12))


f2<-ggplot(subset(cib.rel, ER=="pos"), aes(Response, Monocytes))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+
  scale_x_discrete(limits=c("No NAC", "RD", "pCR"), labels=c("No\nNAC", "RD", "pCR"))+
  labs(x=element_blank(), y ="Relative Monocytes",
       title="ER+ only")+
  scale_color_manual(values=res.col)+
  theme(axis.title.y=element_text(size=12),legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=12))


f3<-ggplot(subset(cib.rel, HER2=="pos"), aes(Response, Monocytes))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+
  scale_x_discrete(limits=c("No NAC", "RD", "pCR"), labels=c("No\nNAC", "RD", "pCR"))+
  labs(x=element_blank(), y ="Relative Monocytes",
       title="HER2+ only")+
  scale_color_manual(values=res.col)+
  theme(axis.title.y=element_text(size=12),legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=12))

##Supplementary figure 2c
f1|f2|f3


#####correlations with clinically measured monocytes
#Supplementary figure 2d
ggscatter(cib.rel, x="Monocytes", y="Monocy", add="reg.line", size=1)+ 
  stat_cor(label.x = .1, label.y = 19)+
  labs(x= "CIBERSORTx Monocytes",
       y= "Clinical Monocytes")



all(rownames(cib.abs)==rownames(cib.rel))
cib.abs$relMon<-cib.rel$Monocytes
cbcdata<-cib.abs[!is.na(cib.abs$MonAbs),]


#Post NAC relative monocytes by RR
stat<-wilcox_test(data=cbcdata, formula=Monocy~RR, p.adjust.method = "fdr")
stat<-stat[stat$p.adj<0.05,]
e2<-ggplot(cbcdata, aes(RR, Monocy))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+  
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR"))+
  labs(title= "Post-NAC", y="Relative Monocytes", x=element_blank())+
  stat_pvalue_manual(stat, label="p.adj.signif", y.position = c(23))+
  scale_color_manual(values=res.col)+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.y=element_text(size=12), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), legend.position = "none")



##Post NAC relative
wilcox_test(data=cbcdata, formula=Monocy~Response, alternative = "g")
stat<-wilcox_test(data=cbcdata, formula=Monocy~Response) %>% add_significance()

d2<-ggplot(cbcdata, aes(Response, Monocy))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=.5)+  
  scale_x_discrete(limits=c("RD", "pCR"))+
  labs(y="Relative Monocytes", x= element_blank(),
       title="Post-NAC")+
  stat_pvalue_manual(stat, label="p.signif", y.position = 24)+
  scale_color_manual(values=res.col)+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.y=element_text(size=12), legend.position = "none")



##Pre tx rel

stat<-wilcox_test(data=cbcdata, formula=PreMonocy~Response)

d1<-ggplot(cbcdata, aes(Response, PreMonocy))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.15, height=0, size=.5)+  
  scale_x_discrete(limits=c("RD", "pCR"))+
  labs(title= "Pre-NAC", y="Relative Monocytes", x=element_blank())+
  scale_color_manual(values=res.col)+
  theme(plot.title = element_text(hjust = 0.5, size=12), axis.title.y=element_text(size=12), legend.position = "none")


cbcdata.complete<-cbcdata[!is.na(cbcdata$PreMonAbs),]
ggpaired(cbcdata.complete, cond1="PreMonocy", cond2 = "Monocy", facet.by = "Response")+
  labs(x="", y="Clinical Monocytes")+
  scale_x_discrete(labels=c("Pre-NAC", "Post-NAC"))+
  facet_grid(~factor(Response, levels=c("RD", "pCR")))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))



ggpaired(cbcdata.complete, cond1="PreMonocy", cond2 = "Monocy")+
  stat_compare_means(paired=TRUE)+labs(x="", y="Relative Monocytes, CBC")+
  scale_x_discrete(labels=c("Pre-NAC", "Post-NAC"))


### ROC Analyses for monocytes ####

sub.cib<-cib.rel[!cib.rel$RR=="No NAC",]
sub.cib$response2<-sub.cib$RR
sub.cib$response2<-revalue(sub.cib$response2, c("pCR"="nR", "RD-nR"="nR", "RD-R"="R"))
table(sub.cib$response2)
#Good vs bad outcome


precrec_obj <- evalmod(scores = sub.cib$Monocytes, labels = sub.cib$Response)
auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.4, y=0.25, label="AUC:0.74")+ 
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(size=14))

precrec_obj <- evalmod(scores = sub.cib$Monocytes, labels = sub.cib$response2)
auc(precrec_obj)
autoplot(precrec_obj, "ROC")


precrec_obj <- evalmod(scores = sub.cib$Monocy, labels = sub.cib$Response)
auc(precrec_obj)
autoplot(precrec_obj, "ROC") +
  annotate("text", x=0.4, y=0.25, label="AUC:0.65")+
  labs(title="Clinical Monocytes, pCR vs RD")+
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(size=14) )

precrec_obj <- evalmod(scores = sub.cib$Monocy, labels = sub.cib$response2)
auc(precrec_obj)
autoplot(precrec_obj, "ROC")



