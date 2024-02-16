# Single-cell RNA-seq analysis - Pseudo-bulk DE analysis with MAST; March 2023
# Venn Diagram for P301S/P301S+K18 1.5 ug/ml pseudo bulk scRNA-seq overlap with AD brain
# Otero-Garcia et al. Neuron 2022: Human AD NFT-/+ dataset (Excitatory neurons only)

#To install packages use this function
BiocManager::install("package")

# Load the packages
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(SingleCellExperiment)
library(AnnotationHub)
library(ensembldb)
library(RColorBrewer)
library(scuttle)
library(scater)
library(edgeR)
library(dplyr)
library(magrittr)
library(MAST)
library(ggVennDiagram)
library(GeneOverlap)
library(viridis)

# Bring in Seurat object
seurat <- readRDS("Excitatory.rds")

table(seurat$SORT)

Idents(seurat) <- "SORT"
DefaultAssay(seurat) <- "RNA"

# with logFC threshold
# AT8 (AD NFT+) vs MAP2 (AD NFT-)
DEGs <- FindMarkers(seurat, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, ident.1 = "AT8", ident.2 = "MAP2")
write.csv(DEGs,"Ex_AT8_vs_MAP2_DEGs.csv")

# AT8 (AD NFT+) vs MAP2control (non-AD)
DEGs <- FindMarkers(seurat, logfc.threshold = 0.1, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, ident.1 = "AT8", ident.2 = "MAP2control")
write.csv(DEGs,"Ex_AT8_vs_MAP2control_DEGs.csv")



# no logFC threshold
# AT8 (AD NFT+) vs MAP2 (AD NFT-)
DEGs <- FindMarkers(seurat, logfc.threshold = 0, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, ident.1 = "AT8", ident.2 = "MAP2")
write.csv(DEGs,"Ex_AT8_vs_MAP2_DEGs_2.csv")

# AT8 (AD NFT+) vs MAP2control (non-AD)
DEGs <- FindMarkers(seurat, logfc.threshold = 0, test.use = "MAST", min.pct = 0.25, only.pos = FALSE, ident.1 = "AT8", ident.2 = "MAP2control")
write.csv(DEGs,"Ex_AT8_vs_MAP2control_DEGs_2.csv")
# ^ couldn't run, memory exhausted



# Calculating overlaps

library(viridis)
library(ggplot2)
library(ggrepel)
library(readxl)
library(ggVennDiagram)
library(GeneOverlap)
library(VennDiagram) 

#P301S+K18 1.5/P301S vs AD_AT8+/AD_AT8- up and down
{
  P301S_dn<-read.table(
    file = "ctrlvs1.5K18_SigDown.csv", 
    header = TRUE, sep =",")
  AT8_dn<-read.table(file = "Ex_AT8_vs_MAP2_SigDown.csv", header = TRUE, sep = ",")
  AT8_dn <- subset(AT8_dn, p_val_adj < 0.05)
  P301S_up<-read.table(
    file = "ctrlvs1.5K18_SigUp.csv", 
    header = TRUE, sep =",")
  AT8_up<-read.table(file = "Ex_AT8_vs_MAP2_SigUp.csv", header = TRUE, sep = ",")
  AT8_up <- subset(AT8_up, p_val_adj < 0.05)

  vd<-list(P301S_dn$Gene, AT8_dn$Gene, P301S_up$Gene, AT8_up$Gene)
  names(vd)<-c("P301S+K18 vs P301S","AT8+ vs AT8-", "P301S+K18 vs P301S","AT8+ vs AT8-")
  
  
  #for list of overlapping genes
  overlap_dn<-intersect(P301S_dn$Gene, AT8_dn$Gene)
  overlap_dn<-data.frame(overlap_dn)
  write.csv(overlap_dn, "AD_1.5_overlap_dn.txt")
  overlap_up<-intersect(P301S_up$Gene, AT8_up$Gene)
  overlap_up<-data.frame(overlap_up)
  write.csv(overlap_up, "AD_1.5_overlap_up.txt")
  
  #global (multi-set plots)
  global<-ggVennDiagram(vd[1:4])
  
  #pairwise ggplots
  P301S_K18vsP301S__AT8pos_AT8neg_Down<-(ggVennDiagram(c(vd[1],vd[2]))+ scale_color_brewer(palette = "Paired")) 
  P301S_K18vsP301S__AT8pos_AT8neg_Up<-(ggVennDiagram(c(vd[3],vd[4]))+ scale_color_brewer(palette = "Paired")) 
  
  #For Figure - up
  venn.diagram(list(A = 1:243, B = 227:326), fill = c("skyblue", "darkblue"), 
               alpha = c(0.5, 0.5), lwd =0, ext.text = FALSE, "venn_diagram_1.5_AD_up.tiff")

  GSEA_up <- read.csv("AD_1.5_overlap_up_GSEA.csv",header=T)
  colnames(GSEA_up) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
  GSEA_up$FDR <- as.numeric(GSEA_up$FDR)
  GSEA_up$logP <- -log10(GSEA_up$FDR)
  GSEA_up$enrich <- paste(GSEA_up$genes_in_data, "/", GSEA_up$genes_in_gsea, sep=" ")
  GSEA_up <- GSEA_up[order(-GSEA_up$logP),]
  GSEA_up$name <-  factor(GSEA_up$name, levels=rev(GSEA_up$name))
  
  ggplot(data=GSEA_up, aes(x=name  , y=logP)) +
    theme_classic() +
    ylab("-Log(FDR)") + xlab(NULL) +
    geom_bar(stat="Identity",  width=0.7, fill="red", alpha=0.8) +
    geom_text(aes(label=enrich), vjust=0.4, 
              hjust=0, size=4, color="black", 
              stat="Identity", y=0.05*max(GSEA_up$logP)) +
    coord_flip() + 
    theme(aspect.ratio = 1.5)+
    ggtitle("P301S+K18 vs P301S/AT8+ vs. AT8- - Up")
  
  
  
  AT8_dn_list <- AT8_dn[,-2:-6]
  writeLines(AT8_dn_list, "AT8_dn.txt")
  
  AT8_up_list <- AT8_up[,-2:-6]
  writeLines(AT8_up_list, "AT8_up.txt")
  
  #GeneOverlap statistics
  
  #Down
  P301S_dn.list<-read.table("ctrlvs1.5K18_SigDown.txt")
  AT8_dn.list<-read.table("AT8_dn.txt")
  
  sapply(P301S_dn.list, length)
  sapply(AT8_dn.list, length)
  
  data(GeneOverlap)
  go.obj <- newGeneOverlap(P301S_dn.list$V1, 
                           AT8_dn.list$V1, 
                           gs.RNASeq)
  go.obj <- testGeneOverlap(go.obj)
  go.obj  # show.
  print(go.obj)  # more details.
  getContbl(go.obj)  # contingency table.
  
  
  
  
  #Up
  P301S_up.list<-read.table("ctrlvs1.5K18_SigUp.txt")
  AT8_up.list<-read.table("AT8_up.txt")
  
  sapply(P301S_up.list, length)
  sapply(AT8_up.list, length)
  
  data(GeneOverlap)
  go.obj <- newGeneOverlap(P301S_up.list$V1, 
                           AT8_up.list$V1, 
                           gs.RNASeq)
  go.obj <- testGeneOverlap(go.obj)
  go.obj  # show.
  print(go.obj)  # more details.
  getContbl(go.obj)  # contingency table.
}

#P301S+K18 1.5/P301S vs AD_AT8+/non-AD up and down
{
  P301S_dn<-read.table(
    file = "ctrlvs1.5K18_SigDown.csv", 
    header = TRUE, sep =",")
  nonAD_dn<-read.table(file = "Ex_AT8_vs_MAP2control_SigDown.csv", header = TRUE, sep = ",")
  nonAD_dn <- subset(nonAD_dn, p_val_adj < 0.05)
  P301S_up<-read.table(
    file = "ctrlvs1.5K18_SigUp.csv", 
    header = TRUE, sep =",")
  nonAD_up<-read.table(file = "Ex_AT8_vs_MAP2control_SigUp.csv", header = TRUE, sep = ",")
  nonAD_up <- subset(nonAD_up, p_val_adj < 0.05)
  
  vd<-list(P301S_dn$Gene, nonAD_dn$Gene, P301S_up$Gene, nonAD_up$Gene)
  names(vd)<-c("P301S+K18 vs P301S","AT8+ vs AT8-", "P301S+K18 vs P301S","AT8+ vs AT8-")
  
  
  #for list of overlapping genes
  overlap_dn<-intersect(P301S_dn$Gene, nonAD_dn$Gene)
  overlap_dn<-data.frame(overlap_dn)
  overlap_up<-intersect(P301S_up$Gene, nonAD_up$Gene)
  overlap_up<-data.frame(overlap_up)
  
  #global (multi-set plots)
  global<-ggVennDiagram(vd[1:4])
  
  #pairwise ggplots
  P301S_K18vsP301S__AT8pos_nonAD_Down<-(ggVennDiagram(c(vd[1],vd[2]))+ scale_color_brewer(palette = "Paired")) 
  P301S_K18vsP301S__AT8pos_nonAD_Up<-(ggVennDiagram(c(vd[3],vd[4]))+ scale_color_brewer(palette = "Paired")) 
  
  
  #For Figure - up
  venn.diagram(list(A = 1:512, B = 482:581), fill = c("skyblue", "darkblue"), 
               alpha = c(0.5, 0.5), lwd =0, ext.text = FALSE, "venn_diagram_1.5_nonAD_up.tiff")
  
  GSEA_up <- read.csv("nonAD_1.5_overlap_up_GSEA.csv",header=T)
  colnames(GSEA_up) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
  GSEA_up$FDR <- as.numeric(GSEA_up$FDR)
  GSEA_up$logP <- -log10(GSEA_up$FDR)
  GSEA_up$enrich <- paste(GSEA_up$genes_in_data, "/", GSEA_up$genes_in_gsea, sep=" ")
  GSEA_up <- GSEA_up[order(-GSEA_up$logP),]
  GSEA_up$name <-  factor(GSEA_up$name, levels=rev(GSEA_up$name))
  
  ggplot(data=GSEA_up, aes(x=name  , y=logP)) +
    theme_classic() +
    ylab("-Log(FDR)") + xlab(NULL) +
    geom_bar(stat="Identity",  width=0.7, fill="red", alpha=0.8) +
    geom_text(aes(label=enrich), vjust=0.4, 
              hjust=0, size=4, color="black", 
              stat="Identity", y=0.05*max(GSEA_up$logP)) +
    coord_flip() + 
    theme(aspect.ratio = 1.5)+
    ggtitle("P301S+K18 vs P301S/AT8+ vs. non-AD - Up")
  
  
  
  
  nonAD_dn_list <- nonAD_dn[,-2:-6]
  writeLines(nonAD_dn_list, "nonAD_dn.txt")
  
  nonAD_up_list <- nonAD_up[,-2:-6]
  writeLines(nonAD_up_list, "nonAD_up.txt")
  
  
  
  
  #GeneOverlap statistics
  
  #Down
  P301S_dn.list<-read.table("ctrlvs1.5K18_SigDown.txt")
  nonAD_dn.list<-read.table("nonAD_dn.txt")
  
  sapply(P301S_dn.list, length)
  sapply(nonAD_dn.list, length)
  
  data(GeneOverlap)
  go.obj <- newGeneOverlap(P301S_dn.list$V1, 
                           nonAD_dn.list$V1, 
                           gs.RNASeq)
  go.obj <- testGeneOverlap(go.obj)
  go.obj  # show.
  print(go.obj)  # more details.
  getContbl(go.obj)  # contingency table.
  
  #Up
  P301S_up.list<-read.table("ctrlvs1.5K18_SigUp.txt")
  nonAD_up.list<-read.table("nonAD_up.txt")
  
  sapply(P301S_up.list, length)
  sapply(nonAD_up.list, length)
  
  data(GeneOverlap)
  go.obj <- newGeneOverlap(P301S_up.list$V1, 
                           nonAD_up.list$V1, 
                           gs.RNASeq)
  go.obj <- testGeneOverlap(go.obj)
  go.obj  # show.
  print(go.obj)  # more details.
  getContbl(go.obj)  # contingency table.
  
}







