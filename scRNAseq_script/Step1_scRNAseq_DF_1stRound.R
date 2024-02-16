#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Alice_sc/DF_1stRound_2023")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Alice_sc/cellranger/Alice_1/outs')
sc = autoEstCont(sc)
Alice_1.counts = adjustCounts(sc)
Alice_1 <- CreateSeuratObject(counts = Alice_1.counts, project = "Alice_1", min.cells = 3, min.features = 200)
rm(Alice_1.counts)
#vizualize QC metrics and filtering====
Alice_1[["percent.mt"]] <- PercentageFeatureSet(object = Alice_1, pattern = "^MT-") #recognize mitochondrial transcripts
all <- Alice_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Alice_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf("Alice_1_QC.pdf", width=12, height=4)
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 30%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 30)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Alice_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Alice_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"Alice_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Alice_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Alice_sc/cellranger/Alice_2/outs')
sc = autoEstCont(sc)
Alice_2.counts = adjustCounts(sc)
Alice_2 <- CreateSeuratObject(counts = Alice_2.counts, project = "Alice_2", min.cells = 3, min.features = 200)
rm(Alice_2.counts)
#vizualize QC metrics and filtering====
Alice_2[["percent.mt"]] <- PercentageFeatureSet(object = Alice_2, pattern = "^MT-") #recognize mitochondrial transcripts
all <- Alice_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Alice_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf("Alice_2_QC.pdf", width=12, height=4)
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 30%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 30)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Alice_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Alice_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"Alice_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Alice_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Alice_sc/cellranger/Alice_3/outs')
sc = autoEstCont(sc)
Alice_3.counts = adjustCounts(sc)
Alice_3 <- CreateSeuratObject(counts = Alice_3.counts, project = "Alice_3", min.cells = 3, min.features = 200)
rm(Alice_3.counts)
#vizualize QC metrics and filtering====
Alice_3[["percent.mt"]] <- PercentageFeatureSet(object = Alice_3, pattern = "^MT-") #recognize mitochondrial transcripts
all <- Alice_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Alice_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
pdf("Alice_3_QC.pdf", width=12, height=4)
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 30%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 30)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Alice_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Alice_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"Alice_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Alice_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################

