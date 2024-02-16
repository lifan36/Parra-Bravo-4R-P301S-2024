#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/Users/lifan/Desktop/data_analysis/Alice_sc/integration_2023")
MC1_integrated <- readRDS(file = "MC1_integrated_PCA_0.1.rds")

Idents(MC1_integrated) <- "seurat_clusters"
pdf("MC1_integrated_umap.pdf", width=4, height=3)
DimPlot(MC1_integrated, reduction = 'umap', label = T)
dev.off()

pdf("MC1_integrated_umap_split_individual.pdf", width=5.5, height=3)
DimPlot(MC1_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 2)
dev.off()

VlnPlot(MC1_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)

Idents(MC1_integrated) <- "orig.ident"
sc_bulk_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(MC1_integrated, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(sc_bulk_TauK18_1.5_vs_TauControl_DEGs, "sc_bulk_TauK18_1.5_vs_TauControl_DEGs.csv")

MC1 <- MC1_integrated
Idents(MC1) <- "seurat_clusters"
Cluster_0 <- subset(MC1, idents = "0")
Cluster_1 <- subset(MC1, idents = "1")
Cluster_2 <- subset(MC1, idents = "2")
Cluster_3 <- subset(MC1, idents = "3")
Cluster_4 <- subset(MC1, idents = "4")
Cluster_5 <- subset(MC1, idents = "5")

Idents(Cluster_0) <- "orig.ident"
Idents(Cluster_1) <- "orig.ident"
Idents(Cluster_2) <- "orig.ident"
Idents(Cluster_3) <- "orig.ident"
Idents(Cluster_4) <- "orig.ident"
Idents(Cluster_5) <- "orig.ident"


setwd("/Users/lifan/Desktop/data_analysis/Alice_sc/integration_2023/DEGs")

Cluster_0_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_0, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_0_TauK18_1.5_vs_TauControl_DEGs, "Cluster_0_TauK18_1.5_vs_TauControl_DEGs.csv")
Cluster_1_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_1, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_1_TauK18_1.5_vs_TauControl_DEGs, "Cluster_1_TauK18_1.5_vs_TauControl_DEGs.csv")
Cluster_2_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_2, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_2_TauK18_1.5_vs_TauControl_DEGs, "Cluster_2_TauK18_1.5_vs_TauControl_DEGs.csv")
Cluster_3_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_3, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_3_TauK18_1.5_vs_TauControl_DEGs, "Cluster_3_TauK18_1.5_vs_TauControl_DEGs.csv")
Cluster_4_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_4, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_4_TauK18_1.5_vs_TauControl_DEGs, "Cluster_4_TauK18_1.5_vs_TauControl_DEGs.csv")
Cluster_5_TauK18_1.5_vs_TauControl_DEGs <- FindMarkers(Cluster_5, ident.1 = "Alice_2", ident.2 = "Alice_1", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(Cluster_5_TauK18_1.5_vs_TauControl_DEGs, "Cluster_5_TauK18_1.5_vs_TauControl_DEGs.csv")









