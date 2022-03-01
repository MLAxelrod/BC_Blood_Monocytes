##### Peripheral blood gene expression predicts outcomes in breast cancer patients ###
##### Differential gene expression, GSEA, signatures ####

library(ggplot2)
library(ggpubr)
library(genefilter)
library(rstatix)
library(DESeq2)
library(apeglm)
library(patchwork)
library(plyr)
library(ggfortify)
library(corrplot)
library(colorspace)
library(ComplexHeatmap)
library(fgsea)
library(tidyverse)
library(limma)
library(systemPipeR)
library(circlize)
library(data.table)
library(ggrepel)

theme_set(theme_classic())
theme_update(axis.text=element_text(color="black", size=10))

#Color vector
res.col<-c("RD-R"="#E64B35FF", "RD-nR"="#4DBBD5FF", "pCR"="#00A087FF", "No NAC"="#3C5488FF")

#zscore function
scaleRow <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
  return(x)
}

filter<-c("RN7SL1", "RN7SL2","HBA1", "HBA2", "HBB", "HBQ1", 
          "HBZ", "HBD", "HBG2",  "HBE1", "HBG1", "HBM",
          "MIR3648-1", "MIR3648-2", "AC104389.6","AC010507.1",
          "SLC25A37",  "SLC4A1", "NRGN", "SNCA", "BNIP3L",  "EPB42", "ALAS2", "BPGM", "OSBP2")
#genes associated with excess RBCs

setwd("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/Manuscript/TidyDataCode/")
read.csv("metadata.csv", row.names=1)->metadata
read.csv("tpm.csv", row.names=1, check.names = FALSE)->tpm
read.csv("featurecounts.csv", row.names=1, check.names = FALSE)->feat
all(rownames(metadata)==colnames(tpm))
all(rownames(metadata)==colnames(feat))



##PCA##
apply(tpm,1,sd)->row.sd
apply(tpm,1,mean)->row.mean
row.sd/(row.mean+1e-6)->row.cv 
tpm[row.cv>0.05,]->tpm.filt 

pcadata<-cbind(t(tpm.filt), metadata)
pca<-prcomp(t(tpm.filt), scale=TRUE)
autoplot(pca, data= pcadata, colour= 'RR') +theme_classic()


#### DeSeq ####
feat<-feat[!rownames(feat) %in% filter, ]

#Take only most variable genes
apply(feat,1,sd)->row.sd
apply(feat,1,mean)->row.mean
row.sd/(row.mean+1e-6)->row.cv 
feat[row.cv>0.05,]->feat.filt 

dds <- DESeqDataSetFromMatrix(countData = feat.filt,
                              colData = metadata,
                              design = ~ Response)
dds$Response <- relevel(dds$Response, ref = "RD") # use this to set control group
dds<-DESeq(dds)
resultsNames(dds)
res<-DESeq2::results(dds, name="Response_pCR_vs_RD")
res<-lfcShrink(dds, coef="Response_pCR_vs_RD", type="apeglm" )
resOrdered <- res[order(res$pvalue),]
summary(res)

############# GSEA #######################
result<-as.data.frame(resOrdered)
sig<-result[result$padj<0.1,]
sig<-na.omit(sig)

#Hallmarks pathways downloaded from database
pathways<- gmtPathways("Y://Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/h.all.v7.1.symbols.gmt")
ranks<-sig$log2FoldChange
names(ranks)<-rownames(sig)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, minSize= 5)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways in pCR vs RD") + 
  theme_minimal()+ theme(axis.text = element_text(color="black"))

###Figure 1b ##
plotEnrichment(pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], ranks)+
  labs(title= "Hallmark Interferon Gamma Response",
       x="Rank", y= "Enrichment Score")+
  theme(axis.title=element_text(size=15), title = element_text(size=14))

plotEnrichment(pathways[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], ranks)+
  labs(title= "Hallmark Interferon Alpha Response",
       x="Rank", y= "Enrichment Score")+
  theme(axis.title=element_text(size=15), title = element_text(size=14))

plotEnrichment(pathways[["HALLMARK_COMPLEMENT"]], ranks)+
  labs(title= "Hallmark Complement",
       x="Rank", y= "Enrichment Score")+
  theme(axis.title=element_text(size=15), title = element_text(size=14))


########GSEA leading edge heatmap ###############

#Venn Diagram of overlap between upregulated pathways
gseagenes<-fgseaResTidy$leadingEdge[1:3]
library(ggvenn)
names(gseagenes)<- c("IFN Gamma", "IFN Alpha", "Complement")

ggvenn(gseagenes, fill_color = c("white", "white", "white"),
       show_elements = FALSE, text_size = 5,
       show_percentage = FALSE)

#foo<-unlist(gseagenes)
#foo2<-unique(foo)
#write.table(foo2, file= "Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/IFNCompGenes.txt", sep="\t")
ic.genes<-read.delim("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/IFNCompGenes.txt", row.names=1, check.names = FALSE)

sig.tpm<-tpm[rownames(tpm) %in% ic.genes$Gene,]
z<-na.omit(scaleRow(sig.tpm))

meta.sort<-metadata[order(metadata$TNBC),]
table(meta.sort$RR)

RR.order<-c("pCR", "RD_nR", "RD_R", "No NAC")
meta.sort<-meta.sort[order(match(meta.sort$RR, RR.order)),]
z.sort<-z[,match(rownames(meta.sort),colnames(z))]
all(rownames(meta.sort)==colnames(z.sort))
#Heatmap annotation
data.frame(meta.sort[["RR"]], meta.sort[["TNBC"]])->annData
colnames(annData)<-c("Response", "TNBC")
annData$Response<-revalue(annData$Response, c("RD_nR" = "RD-nR", "RD_R" = "RD-R"))

f<-HeatmapAnnotation(df= annData, 
                     col = list(Response=c("No NAC" = "#3C5488FF","RD-R" = "#E64B35FF",
                                           "RD-nR" = "#4DBBD5FF", "pCR" ="#00A087FF"),
                                TNBC= c("NO" = "grey20", "YES"= "#F39B7FFF")))

####### Supplementary Figure 1a ###
Heatmap(z.sort, top_annotation = f, row_names_gp = gpar(fontsize=9),
        show_column_names = FALSE, cluster_columns = FALSE, border= TRUE,
        column_split = factor(annData$Response, levels=c("No NAC", "RD-R", "RD-nR", "pCR")),
        cluster_column_slices = FALSE, column_gap= unit(c(3,1,1), "mm"),
        heatmap_legend_param = list(title="Row Z Score", title_position= "lefttop-rot"))



########## Signautres ##########

#Cytotoxic nanoString Signature genes, from CCR publication
nsgenes<-c("FGFBP2", "GNLY", "GZMB", "GZMH", "LAG3", "PDCD1", "NKG7", "HLA-G") 

nsdata<-tpm[rownames(tpm) %in% nsgenes,]
le.data<-tpm[rownames(tpm) %in% ic.genes$Gene,] #leading edge IFN/Comp genes


t1<-as.data.frame(t(scaleRow(nsdata)))
t3<-as.data.frame(t(scaleRow(le.data)))
all(rownames(t1)==rownames(metadata))

metadata$ns.score<-(rowSums(t1[,1:7])-t1$`HLA-G`) / ncol(t1)
metadata$le.score<-rowSums(t3)/ ncol(t3)
metadata$combo<-metadata$le.score-metadata$ns.score #PIRS


metadata$RR<-revalue(metadata$RR, c("RD_nR" = "RD-nR", "RD_R" = "RD-R"))
sub<-metadata[metadata$RR!="No NAC",]


# Scatter plot of scores
## Figure 1c
ggplot(sub, aes(ns.score, le.score))+geom_point(aes(color=RR), size=2)+
  labs(x= "8 Gene Cytotoxic Score",
       y="Interferon/Complement Score", color="Response")+theme_bw()+
  scale_color_manual(values=res.col)+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  geom_vline(xintercept = 0.79, linetype="dashed")+
  geom_hline(yintercept = .86, linetype="dashed")+
  theme(axis.text = element_text(color="black"))


ggscatter(sub, x="ns.score", y="le.score", add="reg.line")+ 
  stat_cor(label.x = .3)



#cytotoxic score score
wilcox_test(data=sub, formula= ns.score ~ RR, p.adjust.method = "fdr")

p1<-ggplot(metadata, aes(RR, ns.score))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y= "Cytotoxic Score", x=element_blank())+
  scale_color_manual(values=res.col)+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=10))




#IFN/Comp leading edge score
wilcox_test(data=sub, formula= le.score ~ RR, p.adjust.method = "fdr")

p2<-ggplot(metadata, aes(RR, le.score))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y= "IFN/Complement Score", x=element_blank())+
  scale_color_manual(values=res.col)+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))


####"Peripheral Immunologic Response Score"

stat<-wilcox_test(data=sub, formula= combo ~ RR, p.adjust.method = "fdr")
plt<-stat[stat$p.adj<0.05,]
p3<-ggplot(metadata, aes(RR, combo))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  stat_pvalue_manual(plt, label="p.adj.signif", y.position = c(4.35,4.8))+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y="PIRS", x=element_blank())+
  scale_color_manual(values=res.col)+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))


##FIGURE 1d
p1|p2|p3

##Supplemental figure 1b

p4<-ggplot(subset(metadata, TNBC=="YES"), aes(RR, combo))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y="PIRS", x=element_blank(),
       title="TNBC only")+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(size=10, color="black", angle = 90, hjust=1,vjust=0.5), title = element_text(size=10))+
  scale_color_manual(values=res.col)


p5<-ggplot(subset(metadata, ER=="pos"), aes(RR, combo))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y="PIRS", x=element_blank(),
       title="ER+ only")+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(size=10, color="black", angle = 90, hjust=1,vjust=0.5), title = element_text(size=10))+
  scale_color_manual(values=res.col)



p6<-ggplot(subset(metadata, HER2=="pos"), aes(RR, combo))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(color=RR),width = 0.05, height=0, size=1)+
  scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
  geom_vline(xintercept=3.5, linetype= "dashed")+
  labs(y="PIRS", x=element_blank(),
       title="HER2+ only")+
  theme(legend.position="none", axis.title = element_text(size=13),
        axis.text.x = element_text(size=10, color="black", angle = 90, hjust=1,vjust=0.5), title = element_text(size=10))+
  scale_color_manual(values=res.col)


p4|p5|p6


##Plot individual genes
foo<-as.data.frame(t(scaleRow(tpm)))
all(rownames(foo)==rownames(metadata))

full.foo<-cbind(metadata, foo)

plt<-function(y, ylab){
  ggplot(full.foo, aes(RR, y))+geom_boxplot(outlier.shape=NA)+
    geom_jitter(aes(color=RR),width = 0.15, height=0, size=0.5)+
    scale_x_discrete(limits=c("RD-R", "RD-nR", "pCR", "No NAC"))+
    geom_vline(xintercept=3.5, linetype= "dashed")+
    labs(x=element_blank(), y=ylab)+
    theme_classic()+ theme(legend.position="none", axis.title.y = element_text(size=13),
                           axis.text = element_text(size=9, color="black"),
                           axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+
    scale_color_manual(values=res.col)
}


##Classical monocyte genes
nac.full<-full.foo[full.foo$RR!="No NAC",]
wilcox_test(data=nac.full, formula= CD14 ~ RR, p.adjust.method = "fdr") #, alternative = "g")
wilcox_test(data=nac.full, formula= CCR2 ~ RR, p.adjust.method = "fdr") #, alternative = "g")

wilcox_test(data=full.foo, formula= CD14 ~ RR, p.adjust.method = "fdr")
p1<-plt(full.foo$CD14, expression(paste(italic(CD14), "  Expression")))+
  geom_bracket(xmin=c("RD-nR"), xmax=c("pCR"), y.position = 2.8, label="*")+ylim(-2,3.1)

plt(full.foo$S100A12, "S100A12")
wilcox_test(data=full.foo, formula= CCR2 ~ RR, p.adjust.method = "fdr")
p2<-plt(full.foo$CCR2, expression(paste(italic(CCR2), "  Expression")))+
  geom_bracket(xmin=c("RD-nR"), xmax=c("pCR"), y.position = 3.3, label="*")+ylim(-1.5, 3.5)

plt(full.foo$CSF3R, "CSF3R")
plt(full.foo$ALOX5AP, "ALOX5AP")

#Non Classical genes
wilcox_test(data=full.foo, formula= FCGR3A ~ RR, p.adjust.method = "fdr")
p3<-plt(full.foo$FCGR3A, expression(paste(italic(FCGR3A), "  Expression")))
wilcox_test(data=full.foo, formula= FCGR3B ~ RR, p.adjust.method = "fdr")
p4<-plt(full.foo$FCGR3B, expression(paste(italic(FCGR3B), "  Expression")))

#(p1|p2)/(p3|p4)
#p2b|p1|p2|p3|p4
p1|p2|p3|p4

plt(full.foo$VMO1, "VMO1")
plt(full.foo$SIGLEC10, "SIGLEC10")


#Intermediate Monocyte genes
plt(full.foo$MARCO, "MARCO")
wilcox_test(data=full.foo, formula= MARCO ~ RR, p.adjust.method = "fdr")

plt(full.foo$GFRA2, "GFRA2")
plt(full.foo$APOBEC3A, "APOBEC3A")
plt(full.foo$TGM2, "TGM2")

#more genes
plt(full.foo$CD274,"CD274")
plt(full.foo$ARG1, "ARG1")
plt(full.foo$IDO1, "IDO1")
plt(full.foo$KLF4, "KLF4")
plt(full.foo$IRF8, "IRF8")
plt(full.foo$CD86, "CD86")
plt(full.foo$SPI1, "SPI1")
plt(full.foo$SELL, "SELL")
plt(full.foo$MRC1, "MRC1")

