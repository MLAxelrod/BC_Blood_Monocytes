########## Whole Blood Single Cell RNA Seq ##############################
#Pt1 = pt 3399; metaplastic TNBC with pCR
#Pt2 = pt 2020-01; TNBC with pCR, newly collected

library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(tidyverse)

combo<-readRDS("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/FilesforNACpaper/scRNAseq/combinedSeurat.rds")

#Simplify cell types
table(combo$CellTypes)
Idents(combo)<-"CellTypes"
combo<-RenameIdents(combo, "Erythrocytes"="Other", "DC"="Other", "Eosinophils"="Other", "HSC"="Other", "Neutrophils"="Other")

#Identities DimPlot
combo$CellTypes<-Idents(combo)
ct<-DimPlot(combo, reduction = "umap", label=TRUE, label.box = TRUE, repel=TRUE)+scale_color_npg()+scale_fill_manual(values=rep("white",6))+NoLegend()
ct

###Cytotoxic 8 Gene Sig###
gene.set<-c("FGFBP2", "GNLY", "GZMB", "GZMH", "NKG7","LAG3", "PDCD1")

# Get mean expression of genes of interest per cell
mean.exp<-rowMeans(FetchData(combo, gene.set), na.rm=TRUE) - rowMeans(FetchData(combo, "HLA-G"), na.rm=TRUE)

# Add mean expression values in 'object@meta.data$cytotoxic.score'
if (all(names(x = mean.exp) == rownames(x = combo@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'combo@meta.data':\n", 
      "adding gene set mean expression values in 'combo@meta.data$cytotoxic.score'")
  combo@meta.data$cytotoxic.score <- mean.exp
}

#cytotoxic score only plots
FeaturePlot(object = combo, features = "cytotoxic.score", cols=c("darkgrey", "grey","orangered","red"))+
  ggtitle("Cytotoxic Score")
  
a1<-VlnPlot(combo, features = c("cytotoxic.score"), group.by="CellTypes",pt.size =0, sort=FALSE, log=FALSE)+
  scale_x_discrete(limits=c("NK cells", "CD8+ T-cells", "CD4+ T-cells", "B-cells", "Monocytes", "Other"))+
  ggtitle("Cytotoxic Score")+xlab("")+ 
  geom_boxplot(width=0.1, outlier.shape = NA)+NoLegend()



####IFN/Comp Gene score ####
IFNComp<-read.delim("Y:/Balko/Data/Maggie Axelrod/BreastCancerPBMC/bloodRNAseq/IFNCompGenes.txt", row.names=1, check.names = FALSE)

mean.exp<-rowMeans(FetchData(combo, IFNComp$Gene), na.rm=TRUE)

# Add mean expression values in 'object@meta.data$cytotoxic.score'
if (all(names(x = mean.exp) == rownames(x = combo@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'combo@meta.data':\n", 
      "adding gene set mean expression values in 'combo@meta.data$ifn.comp.score'")
  combo@meta.data$ifn.comp.score <- mean.exp
}

# IFN/compe score only plots
FeaturePlot(object = combo, features = "ifn.comp.score", cols=c("darkgrey", "grey","orangered","red")) + 
  ggtitle("IFN/Complement Score")

a2<-VlnPlot(combo, features = c("ifn.comp.score"), group.by="CellTypes", pt.size =0, sort=FALSE, log=FALSE)+
  scale_x_discrete(limits=c("NK cells", "CD8+ T-cells", "CD4+ T-cells", "B-cells", "Monocytes", "Other"))+
  ggtitle("IFN/Complement Score")+xlab("")+ 
  geom_boxplot(width=0.2, outlier.shape = NA)+
  NoLegend()


###Combined gene sets figure ######
p<-FeaturePlot(object = combo, features = c("cytotoxic.score", "ifn.comp.score"),
               cols=c("red", "blue"),
               blend=TRUE, combine = TRUE)

p1<-p[[1]]+ggtitle("Cytotoxic Score")
p2<-p[[2]]+ggtitle("IFN/Complement Score")
p3<-p[[3]]+ggtitle("Both Scores")
p4<-p[[4]]+labs(x="Cytotoxic Score", y="IFN/Complement Score", title=element_blank())


(ct|p1|p2|p3|p4)

##PIRS violin (IFN - cyto)

# Get mean expression of genes of interest per cell
mean.exp<-rowMeans(FetchData(combo, "ifn.comp.score"), na.rm=TRUE) - rowMeans(FetchData(combo, "cytotoxic.score"), na.rm=TRUE)


if (all(names(x = mean.exp) == rownames(x = combo@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'combo@meta.data':\n", 
      "adding gene set mean expression values in 'combo@meta.data$PIRS'")
  combo@meta.data$PIRS <- mean.exp
}

a3<-VlnPlot(combo, features = c("PIRS"), group.by="CellTypes", pt.size =0, sort=FALSE, log=FALSE)+
  scale_x_discrete(limits=c("NK cells", "CD8+ T-cells", "CD4+ T-cells", "B-cells", "Monocytes", "Other"))+
  ggtitle("PIRS")+xlab("")+ 
  geom_boxplot(width=0.2, outlier.shape = NA)+
  NoLegend()


a1/a2/a3
