
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
setwd("/athena/ganlab/scratch/lif4001/Alice_sc/DF_2ndRound_2023")
Alice_1 <- readRDS(file = "Alice_1_singlets_PCA.rds")
Alice_2 <- readRDS(file = "Alice_2_singlets_PCA.rds")
#Alice_3 <- readRDS(file = "Alice_3_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Alice_sc/integration_2023")

MC1 <- c(Alice_1, Alice_2)
anchors_MC1 <- FindIntegrationAnchors(object.list = MC1, dims = 1:30)
MC1_integrated <- IntegrateData(anchorset = anchors_MC1, dims = 1:30)
rm(Alice_1, Alice_2, MC1)

#saveRDS(MC1_integrated, file = "MC1_integrated.rds")

#MC1_integrated <- readRDS("MC1_integrated.rds")
DefaultAssay(MC1_integrated) <- 'integrated'

#MC1_integrated <- NormalizeData(MC1_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
#MC1_integrated <- FindVariableFeatures(MC1_integrated, selection.method = "vst", nfeatures = 3000)

MC1_integrated <- ScaleData(MC1_integrated, verbose = FALSE)
MC1_integrated <- RunPCA(MC1_integrated, features = VariableFeatures(object = MC1_integrated), verbose = FALSE)

MC1_integrated <- FindNeighbors(MC1_integrated, dims = 1:20)
MC1_integrated <- FindClusters(MC1_integrated, resolution = 0.075)
MC1_integrated <- RunUMAP(MC1_integrated, dims = 1: 20)

str(MC1_integrated)

DefaultAssay(MC1_integrated) <- 'RNA'
MC1_integrated <- NormalizeData(MC1_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
MC1_integrated <- ScaleData(MC1_integrated, features = rownames(MC1_integrated))


pdf("MC1_integrated_umap.pdf", width=4, height=3)
DimPlot(MC1_integrated, reduction = 'umap', label = T)
dev.off()
pdf("MC1_integrated_umap_split_individual.pdf", width=6, height=3)
DimPlot(MC1_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 2)
dev.off()

pdf("MC1_QC.pdf", width=9, height=4)
Idents(MC1_integrated) <- "orig.ident"
VlnPlot(object = MC1_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(MC1_integrated, file = 'MC1_integrated_PCA_0.1.rds')

#MC1_integrated <- readRDS("MC1_integrated_PCA_0.1.rds")

DefaultAssay(MC1_integrated) <- 'RNA'
pdf("MC1_integrated_umap_test.pdf", width=4, height=3)
DimPlot(MC1_integrated, reduction = 'umap', label = T)
dev.off()


#Add marker genes

pdf("MC1_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("SYT1","SNAP25","GRIN1","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2","PLP1", "MBP", "MOBP","AQP4","GFAP", 
           "CD74","CSF1R","C3","PDGFRA","VCAN","EBF1","IGFBP7","FLT1","CLDN5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = MC1_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

DefaultAssay(MC1_integrated) <- 'RNA'

MC1_markers <- FindAllMarkers(MC1_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(MC1_markers, "MC1_markers.csv")

write.csv(table(MC1_integrated$seurat_clusters, MC1_integrated$orig.ident), "cell_counts_cluster_sample.csv")



