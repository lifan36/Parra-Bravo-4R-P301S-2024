# 4R/4R-P301S -/+ Fibril Treatment Bulk RNA-Seq; September 2021

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(data.table)
library(biomaRt)
library(reshape2)
library(circlize)
library(stringr)


# Gene names and remove duplicates ####
ensembl <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", useEnsembl("ENSEMBL_MART_ENSEMBL", mirror = 'useast'))

df <- read.csv("raw_counts_genes.clean.csv")
genes <- df$ensembl_gene_id
G_list <- getBM(attributes= c("ensembl_gene_id", "external_gene_name",filters= "ensembl_gene_id"),values=genes,mart=mart)
total <- merge(df,G_list,by="ensembl_gene_id")
total <- total[,-1]
total <- total[,-21]
total <- subset(total, !duplicated(total$external_gene_name))
row.names(total) <- total$external_gene_name
total <- total[,-20]
write.csv(total,"readcounts_genenames.csv")

# DESeq #### use commented code to plot only genotype or treatment
data <- read.csv("readcounts_genenames.csv", header=T, row.names=1)
#data <- data[c(1,2,3,4,5,11,12,13,14,15)]
#data <- data[c(6,7,8,9,10,16,17,18,19)]
filtered = data[rowSums(data)>15,]
sample <- colnames(data)
genotype <- c(rep("Tau4R",10),rep("P301S",9))
treatment <- c(rep("CTL",5),rep("K18",5),rep("CTL",5),rep("K18",4))
#genotype <- c(rep("Tau4R",5),rep("P301S",5))
#treatment <- c(rep("CTL",5),rep("CTL",5))
#genotype <- c(rep("Tau4R",5),rep("P301S",4))
#treatment <- c(rep("K18",5),rep("K18",4))
genotype_treatment <- paste(genotype,treatment,sep="_")
meta <- data.frame(sample=sample, genotype=genotype, treatment=treatment, genotype_treatment=genotype_treatment)

all(colnames(filtered) %in% meta$sample)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~genotype_treatment)
rld <- rlog(dds, blind = T)

plotPCA(rld, intgroup = "genotype_treatment", ntop = 500)+
  geom_text_repel(aes(label = genotype_treatment))+
  theme_classic()

rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
pheatmap(rld_cor)

dds <- DESeq(dds)

# Pairwise Comparisons ####
# 1. P301S vs. 4R
contrast_P301Svs4R <- c("genotype_treatment","P301S_CTL","Tau4R_CTL")
res_P301Svs4R_unshrunken <- results(dds,contrast=contrast_P301Svs4R,alpha=0.05)
res_P301Svs4R <- lfcShrink(dds,contrast=contrast_P301Svs4R,res=res_P301Svs4R_unshrunken, type="normal")
write.csv(res_P301Svs4R, "DE_P301Svs4R.csv")

# 2. P301S K18 vs. 4R K18
contrast_P301SK18vs4RK18 <- c("genotype_treatment","P301S_K18","Tau4R_K18")
res_P301SK18vs4RK18_unshrunken <- results(dds,contrast=contrast_P301SK18vs4RK18,alpha=0.05)
res_P301SK18vs4RK18 <- lfcShrink(dds,contrast=contrast_P301SK18vs4RK18,res=res_P301SK18vs4RK18_unshrunken, type="normal")
write.csv(res_P301SK18vs4RK18, "DE_P301SK18vs4RK18.csv")

# 3. P301S K18 vs.  P301S CTL
contrast_P301SK18vsP301S <- c("genotype_treatment","P301S_K18","P301S_CTL")
res_P301SK18vsP301S_unshrunken <- results(dds,contrast=contrast_P301SK18vsP301S,alpha=0.05)
res_P301SK18vsP301S <- lfcShrink(dds,contrast=contrast_P301SK18vsP301S,res=res_P301SK18vsP301S_unshrunken, type="normal")
write.csv(res_P301SK18vsP301S, "DE_P301SK18vsP301S.csv")

# 4. 4R K18 vs.  4R CTL
contrast_4RK18vs4R <- c("genotype_treatment","Tau4R_K18","Tau4R_CTL")
res_4RK18vs4R_unshrunken <- results(dds,contrast=contrast_4RK18vs4R,alpha=0.05)
res_4RK18vs4R <- lfcShrink(dds,contrast=contrast_4RK18vs4R,res=res_4RK18vs4R_unshrunken, type="normal")
write.csv(res_4RK18vs4R, "DE_4RK18vs4R.csv")

# Volcano Plots ####
# 1. P301S vs. 4R
P301Svs4R <- read.csv("DE_P301Svs4R.csv",header=T,row.names=1)
P301Svs4R$color[P301Svs4R$log2FoldChange > 0 & P301Svs4R$padj < 0.05] <- "red"
P301Svs4R$color[P301Svs4R$log2FoldChange < 0 & P301Svs4R$padj < 0.05] <- "blue"
P301Svs4R$color[P301Svs4R$log2FoldChange < 0 & P301Svs4R$padj > 0.05] <- "grey"
P301Svs4R$color[P301Svs4R$log2FoldChange < 0 & P301Svs4R$padj > 0.05] <- "grey"
labels <- subset(P301Svs4R, log2FoldChange < -1 | log2FoldChange > 1 & padj < 0.05)

ggplot(P301Svs4R, aes(x = log2FoldChange, y = -log(padj), color = color))+
  geom_point()+
  theme_classic()+
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("P301S vs. 4R")+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels), max.overlaps = getOption("ggrepel.max.overlaps", default = 50)))

# 2. P301S K18 vs. 4R K18
P301SK18vs4RK18 <- read.csv("DE_P301SK18vs4RK18.csv",header=T,row.names=1)
P301SK18vs4RK18$color[P301SK18vs4RK18$log2FoldChange > 0 & P301SK18vs4RK18$padj < 0.05] <- "red"
P301SK18vs4RK18$color[P301SK18vs4RK18$log2FoldChange < 0 & P301SK18vs4RK18$padj < 0.05] <- "blue"
P301SK18vs4RK18$color[P301SK18vs4RK18$log2FoldChange < 0 & P301SK18vs4RK18$padj > 0.05] <- "grey"
P301SK18vs4RK18$color[P301SK18vs4RK18$log2FoldChange < 0 & P301SK18vs4RK18$padj > 0.05] <- "grey"
labels <- subset(P301SK18vs4RK18, log2FoldChange < -1 | log2FoldChange > 1 & padj < 0.05)

ggplot(P301SK18vs4RK18, aes(x = log2FoldChange, y = -log(padj), color = color))+
  geom_point()+
  theme_classic()+
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("P301S K18 vs. 4R K18")+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels), max.overlaps = getOption("ggrepel.max.overlaps", default = 50)))

# 3. P301S K18 vs.  P301S CTL
P301SK18vsP301S <- read.csv("DE_P301SK18vsP301S.csv",header=T,row.names=1)
P301SK18vsP301S$color[P301SK18vsP301S$log2FoldChange > 0 & P301SK18vsP301S$padj < 0.05] <- "red"
P301SK18vsP301S$color[P301SK18vsP301S$log2FoldChange < 0 & P301SK18vsP301S$padj < 0.05] <- "blue"
P301SK18vsP301S$color[P301SK18vsP301S$log2FoldChange < 0 & P301SK18vsP301S$padj > 0.05] <- "grey"
P301SK18vsP301S$color[P301SK18vsP301S$log2FoldChange < 0 & P301SK18vsP301S$padj > 0.05] <- "grey"
labels <- subset(P301SK18vsP301S, log2FoldChange < -1 | log2FoldChange > 1 & padj < 0.05)

ggplot(P301SK18vsP301S, aes(x = log2FoldChange, y = -log(padj), color = color))+
  geom_point()+
  theme_classic()+
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("P301S K18 vs. P301S CTL")+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels)))

# 4. 4R K18 vs.  4R CTL
RK18vs4R <- read.csv("DE_4RK18vs4R.csv",header=T,row.names=1)
RK18vs4R$color[RK18vs4R$log2FoldChange > 0 & RK18vs4R$padj < 0.05] <- "red"
RK18vs4R$color[RK18vs4R$log2FoldChange < 0 & RK18vs4R$padj < 0.05] <- "blue"
RK18vs4R$color[RK18vs4R$log2FoldChange < 0 & RK18vs4R$padj > 0.05] <- "grey"
RK18vs4R$color[RK18vs4R$log2FoldChange < 0 & RK18vs4R$padj > 0.05] <- "grey"
labels <- subset(RK18vs4R, log2FoldChange < -1 | log2FoldChange > 1 & padj < 0.05)

ggplot(RK18vs4R, aes(x = log2FoldChange, y = -log(padj), color = color))+
  geom_point()+
  theme_classic()+
  xlab("log2FoldChange")+
  ylab("-log(padj)")+
  ggtitle("4R K18 vs. 4R CTL")+
  geom_vline(xintercept=0, linetype="dashed", col="gray") + geom_hline(yintercept=-log10(0.05), linetype="dashed", col="gray")+
  scale_color_manual(values = c("red" = "red", "blue" = "blue","grey"="grey"),
                     name = "DEG",
                     breaks = c("blue","red","grey"),
                     labels = c("Downregulated","Upregulated","No Change"))+
  geom_text_repel(labels, mapping = aes(label = row.names(labels)))



# GSEA #### Enter list into https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp to perform GO
# Significant Gene Lists
P301Svs4R <- read.csv("DE_P301Svs4R.csv",header=T,row.names=1)
P301Svs4R.up <- subset(P301Svs4R, log2FoldChange > 0 & padj < 0.05)
P301Svs4R.dn <- subset(P301Svs4R, log2FoldChange < 0 & padj < 0.05)
write.csv(P301Svs4R.up, "P301Svs4R_SigUp_fromGillian.csv")
write.csv(P301Svs4R.dn, "P301Svs4R_SigDown_fromGillian.csv")

P301S_K18vs4R_K18 <- read.csv("DE_P301SK18vs4RK18.csv",header=T,row.names=1)
P301S_K18vs4R_K18.up <- subset(P301S_K18vs4R_K18, log2FoldChange > 0 & padj < 0.05)
P301S_K18vs4R_K18.dn <- subset(P301S_K18vs4R_K18, log2FoldChange < 0 & padj < 0.05)
write.csv(P301S_K18vs4R_K18.up, "P301S_K18vs4R_K18_SigUp.csv")
write.csv(P301S_K18vs4R_K18.dn, "P301S_K18vs4R_K18_SigDown.csv")

P301S_K18vsP301S <- read.csv("DE_P301SK18vsP301S.csv",header=T,row.names=1)
P301S_K18vsP301S.up <- subset(P301S_K18vsP301S, log2FoldChange > 0 & padj < 0.05)
P301S_K18vsP301S.dn <- subset(P301S_K18vsP301S, log2FoldChange < 0 & padj < 0.05)
# no significantly changed genes

RK18vs4R <- read.csv("DE_4RK18vs4R.csv",header=T,row.names=1)
RK18vs4R.up <- subset(RK18vs4R, log2FoldChange > 0 & padj < 0.05)
RK18vs4R.dn <- subset(RK18vs4R, log2FoldChange < 0 & padj < 0.05)
# no significantly changed genes

# P301S vs. 4R (CTL) Down
GSEA_dn <- read.csv("P301Svs4R_Down_GSEA.csv",header=T)
colnames(GSEA_dn) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
GSEA_dn$FDR <- as.numeric(GSEA_dn$FDR)
GSEA_dn$logP <- -log10(GSEA_dn$FDR)
GSEA_dn$enrich <- paste(GSEA_dn$genes_in_data, "/", GSEA_dn$genes_in_gsea, sep=" ")
GSEA_dn <- GSEA_dn[order(-GSEA_dn$logP),]
GSEA_dn$name <-  factor(GSEA_dn$name, levels=rev(GSEA_dn$name))

ggplot(data=GSEA_dn, aes(x=name  , y=logP)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="blue", alpha=0.8) +
  geom_text(aes(label=enrich), vjust=0.4, 
            hjust=0, size=4, color="black", 
            stat="Identity", y=0.05*max(GSEA_dn$logP)) +
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  ggtitle("P301S vs. 4R - Down")



# P301S vs. 4R (CTL) Up
GSEA_up <- read.csv("P301Svs4R_Up_GSEA.csv",header=T)
colnames(GSEA_up) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
GSEA_up$FDR <- as.numeric(GSEA_up$FDR)
GSEA_up$logP <- -log10(GSEA_up$FDR)
GSEA_up$enrich <- paste(GSEA_up$genes_in_data, "/", GSEA_up$genes_in_gsea, sep=" ")
GSEA_up <- GSEA_up[order(GSEA_up$FDR),]
GSEA_up$Name <-  gsub("GO_", "", GSEA_up$name)
GSEA_up$Name <- gsub("*_", " ", GSEA_up$Name)
# for some reason, you need this section twice
GSEA_up$Name <-  factor(GSEA_up$Name, levels=(GSEA_up$Name))
GSEA_up$Name <-  factor(GSEA_up$Name, levels=rev(GSEA_up$Name))
GSEA_up <- head(GSEA_up, 10)
ggplot(data=GSEA_up, aes(x=Name  , y=logP)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="red", alpha=0.8) +
  geom_text(aes(label=enrich), vjust=0.4,
            hjust=0, size=4, color="black",
            stat="Identity", y=0.05*max(GSEA_up$logP)) +
  coord_flip() +
  theme(aspect.ratio = 1.5)+ ggtitle("P301S vs. 4R - Up")


# P301S vs. 4R (K18) Down
GSEA_dn <- read.csv("P301Svs4R_K18_Down_GSEA.csv",header=T)
colnames(GSEA_dn) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
GSEA_dn$FDR <- as.numeric(GSEA_dn$FDR)
GSEA_dn$logP <- -log10(GSEA_dn$FDR)
GSEA_dn$enrich <- paste(GSEA_dn$genes_in_data, "/", GSEA_dn$genes_in_gsea, sep=" ")
GSEA_dn <- GSEA_dn[order(-GSEA_dn$logP),]
GSEA_dn$name <-  factor(GSEA_dn$name, levels=rev(GSEA_dn$name))

ggplot(data=GSEA_dn, aes(x=name  , y=logP)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="blue", alpha=0.8) +
  geom_text(aes(label=enrich), vjust=0.4, 
            hjust=0, size=4, color="black", 
            stat="Identity", y=0.05*max(GSEA_dn$logP)) +
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  ggtitle("P301S_K18 vs. 4R_K18 - Down")

# P301S vs. 4R (K18) Up
GSEA_up <- read.csv("P301Svs4R_K18_Up_GSEA.csv",header=T)
colnames(GSEA_up) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
GSEA_up$FDR <- as.numeric(GSEA_up$FDR)
GSEA_up$logP <- -log10(GSEA_up$FDR)
GSEA_up$enrich <- paste(GSEA_up$genes_in_data, "/", GSEA_up$genes_in_gsea, sep=" ")
GSEA_up <- GSEA_up[order(-GSEA_up$logP),]
GSEA_up$Name <-  factor(GSEA_up$name, levels=rev(GSEA_up$name))

ggplot(data=GSEA_up, aes(x=name  , y=logP)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="red", alpha=0.8) +
  geom_text(aes(label=enrich), vjust=0.4, 
            hjust=0, size=4, color="black", 
            stat="Identity", y=0.05*max(GSEA_up$logP)) +
  coord_flip() + 
  theme(aspect.ratio = 1.5)+
  ggtitle("P301S_K18 vs. 4R_K18 - Up")




# Heatmap for P301Sv4R Trafficking DEGs from GSEA
colors <- rev(brewer.pal(10,"RdBu"))
annotations <- as.data.frame(colData(dds)["genotype_treatment"])
counts <- counts(dds,normalized = T)

traffic <- read.csv("P301Sv4R_DEG_SigDown_trafficking_genelist.csv", header=T)
row.names(traffic) <- traffic$Gene
overlap.traffic <- subset(counts, row.names(counts) %in% row.names(traffic))

pheatmap(overlap.traffic[,c(1:5,11:15)],
         show_rownames = T,
         color = colors,
         main = '4R-P301S vs 4R',
         scale = "row")
