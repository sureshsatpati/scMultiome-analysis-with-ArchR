library(dplyr)
library(Seurat)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(sctransform)
library(SingleCellExperiment)
library(cluster)
library(factoextra)
library(intrinsicDimension)
library(DoubletFinder)

rm()
gc()

D1_Tex_1 <- Read10X(data.dir = "/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/1.Raw/D1_Tex_1/D1_Tex_1/outs/filtered_feature_bc_matrix")
D1_Tex_2 <- Read10X(data.dir = "/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/1.Raw/D1_Tex_2/D1_Tex_2/outs/filtered_feature_bc_matrix")
D1_Tex_3 <- Read10X(data.dir = "/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/1.Raw/D1_Tex_3/D1_Tex_3/outs/filtered_feature_bc_matrix")
D1_Tex_4 <- Read10X(data.dir = "/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/1.Raw/D1_Tex_4/D1_Tex_4/outs/filtered_feature_bc_matrix")
D1_Tex_5 <- Read10X(data.dir = "/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/1.Raw/D1_Tex_5/D1_Tex_5/outs/filtered_feature_bc_matrix")

sdata.D1_Tex_1 <- CreateSeuratObject(counts = D1_Tex_1$'Gene Expression', project = "D1_Tex_1")
sdata.D1_Tex_2 <- CreateSeuratObject(counts = D1_Tex_2$'Gene Expression', project = "D1_Tex_2")
sdata.D1_Tex_3 <- CreateSeuratObject(counts = D1_Tex_3$'Gene Expression', project = "D1_Tex_3")
sdata.D1_Tex_4 <- CreateSeuratObject(counts = D1_Tex_4$'Gene Expression', project = "D1_Tex_4")
sdata.D1_Tex_5 <- CreateSeuratObject(counts = D1_Tex_5$'Gene Expression', project = "D1_Tex_5")

sdata.D1_Tex_1$type = "D1_Tex_1"
sdata.D1_Tex_2$type = "D1_Tex_2"
sdata.D1_Tex_3$type = "D1_Tex_3"
sdata.D1_Tex_4$type = "D1_Tex_4"
sdata.D1_Tex_5$type = "D1_Tex_5"


sdata.D1_Tex_1$condition = "Naive"
sdata.D1_Tex_2$condition = "Activated"
sdata.D1_Tex_3$condition = "Exhausted1"
sdata.D1_Tex_4$condition = "Exhausted2"
sdata.D1_Tex_5$condition = "Exhausted3"


alldata <- merge(sdata.D1_Tex_1, c(sdata.D1_Tex_2,sdata.D1_Tex_3,sdata.D1_Tex_4,sdata.D1_Tex_5), add.cell.ids = c("sdata.D1_Tex_1", "sdata.D1_Tex_2", "sdata.D1_Tex_3", "sdata.D1_Tex_4", "sdata.D1_Tex_5"))
saveRDS(alldata, "RNA-raw-all.rds")

alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

pdf("RNA_Raw_plot.pdf")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + geom_boxplot(width=0.1,fill="white") + NoLegend()
dev.off()

names(alldata@meta.data)

table(alldata$type)

selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 5]

selected_c <- WhichCells(alldata, expression = nFeature_RNA > 400)

first.filt <- subset(alldata, features = selected_f, cells = selected_c)

table(first.filt$orig.ident)

table(alldata$type)

saveRDS(first.filt, "filtered_RNAseq_data_after_first_filt.rds")













