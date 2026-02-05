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
set.seed(1)
rm()
gc()

library(ArchR)
addArchRThreads(threads = 1)

library(Seurat)
library(SingleCellExperiment)
library("BSgenome.Hsapiens.UCSC.hg38")

addArchRGenome("hg38")
#change1
#_BEGIN_ multiple lines comment
#if (FALSE) {
#inputFiles <- c("/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_1/atac_fragments.tsv.gz","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_2/atac_fragments.tsv.gz","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_3/atac_fragments.tsv.gz","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_4/atac_fragments.tsv.gz","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_5/atac_fragments.tsv.gz")
#names(inputFiles) <- c("D1_Tex_1","D1_Tex_2","D1_Tex_3","D1_Tex_4","D1_Tex_5")

#inputFiles




#ArrowFiles <- createArrowFiles(
#inputFiles = inputFiles,
#sampleNames = names(inputFiles),
#force = TRUE
#)


#change2
#seRNA_new <- import10xFeatureMatrix(input = c("/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_1/filtered_feature_bc_matrix.h5","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_2/filtered_feature_bc_matrix.h5","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_3/filtered_feature_bc_matrix.h5","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_4/filtered_feature_bc_matrix.h5","/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/input/D1_Tex_5/filtered_feature_bc_matrix.h5"), names = c('D1_Tex_1','D1_Tex_2','D1_Tex_3','D1_Tex_4','D1_Tex_5'))

#seRNA_new

#proj <- ArchRProject(ArrowFiles)

#below command we have to repeat because we have to subset so not using here, but u can use it will only increase your run time.
#proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA_new, force = TRUE, strictMatch = FALSE)

#archrCells <- ArchR::getCellNames(proj)
#scRNACells <- colnames(seRNA_new)
#commonCells <- intersect(archrCells, scRNACells)
#archrproj <- subsetArchRProject(
#proj,
#commonCells,
#force = TRUE
#)

#archrproj <- addGeneExpressionMatrix(input = archrproj, seRNA = seRNA_new, force = TRUE, strictMatch = TRUE)

#archrproj <- archrproj[archrproj$TSSEnrichment > 5 & archrproj$nFrags > 2500 & !is.na(archrproj$Gex_nUMI)]

#archrproj <- addDoubletScores(archrproj)

#archrproj <- filterDoublets(archrproj)

#archrproj <- addIterativeLSI(
#ArchRProj = archrproj,
#clusterParams = list(
#resolution = 0.2,
#sampleCells = 10000,
#n.start = 10
#),
#saveIterations = FALSE,
#useMatrix = "TileMatrix",
#depthCol = "nFrags",
#name = "LSI_ATAC"
#)

#archrproj <- addIterativeLSI(
#ArchRProj = archrproj,
#clusterParams = list(
#resolution = 0.2,
#sampleCells = 10000,
#n.start = 10
#),
#saveIterations = FALSE,
#useMatrix = "GeneExpressionMatrix",
#depthCol = "Gex_nUMI",
#varFeatures = 2500,
#firstSelection = "variable",
#binarize = FALSE,
#name = "LSI_RNA"
#)

#Combined Dims
#archrproj <- addCombinedDims(archrproj,reducedDims=c("LSI_ATAC","LSI_RNA"),name="LSI_Combined")

#colnames(archrproj@reducedDims$LSI_Combined$matRD) <- paste0('LSI',1:length(colnames(getReducedDims(archrproj,'LSI_Combined'))))

#UMAPs
#archrproj <- addUMAP(archrproj,reducedDims="LSI_ATAC",name="UMAP_ATAC",minDist=0.8,force=TRUE)
#p1 <- plotEmbedding(archrproj,name="Clusters",embedding="UMAP_ATAC",size=1.5,labelAsFactors=F,labelMeans=F)

#archrproj <- addUMAP(archrproj,reducedDims="LSI_RNA",name="UMAP_RNA",minDist=0.8,force=TRUE)
#p2 <- plotEmbedding(archrproj,name="Clusters",embedding="UMAP_RNA",size = 1.5,labelAsFactors=F,labelMeans=F)

#archrproj <- addUMAP(archrproj,reducedDims="LSI_Combined",name="UMAP_Combined",minDist=0.8,force=TRUE)


#p3 <- plotEmbedding(archrproj, name="Clusters",embedding="UMAP_Combined",size=1.5,labelAsFactors=F,labelMeans=F)

#archrproj <- addClusters(archrproj,reducedDims="LSI_Combined",name="Clusters",resolution=0.4,force = TRUE)



#p1 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)

#p2 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)

#p3 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)

#plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", ArchRProj = archrproj, addDOC = FALSE)


#markersGS <- getMarkerFeatures(
#ArchRProj = archrproj, 
#useMatrix = "GeneScoreMatrix", 
#groupBy = "Clusters",
#bias = c("TSSEnrichment", "log10(nFrags)"),
#testMethod = "wilcoxon"
#)

#markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
#heatmapGS <- plotMarkerHeatmap(
#seMarker = markersGS,
#cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
#transpose = TRUE
#)
#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#write.table(markerList, file="AllvsAllmarkerList.Differentialmarkers_FDR_0.01_and_Log2FC_1.25.xls", quote=FALSE, sep="\t", row.names=FALSE)

#plotPDF(heatmapGS, name = "Marker-Heatmap-FDR-0.01-FC-1.25", width = 30, height = 16, ArchRProj = archrproj, addDOC = FALSE)


#saveArchRProject(ArchRProj = archrproj, outputDirectory = "Save-Proj-after-differential_marker", load = FALSE)

#archrproj <- addGroupCoverages(ArchRProj = archrproj, groupBy = "Clusters")

#archrproj <- addReproduciblePeakSet(
#ArchRProj = archrproj,
#groupBy = "Clusters",
#pathToMacs2 = "/risapps/rhel7/python/2.7.13/bin/macs2"
#)


#getPeakSet(archrproj)
#PeakSet <- getPeakSet(archrproj)
#write.table(PeakSet, file="PeakSet.xls", quote=FALSE, sep="\t", row.names=FALSE)

#archrproj <- addPeakMatrix(archrproj)

#getAvailableMatrices(archrproj)


#saveArchRProject(ArchRProj = archrproj, outputDirectory = "ArchR-Save-afteradded-PeakMatrix", load = FALSE)
#}
#Multiple lines comment _END_ here
###Change path #change3
archrprojPM <- loadArchRProject(path = "ArchR-Save-afteradded-PeakMatrix", force = FALSE, showLogo = TRUE)


markersPeaks <- getMarkerFeatures(ArchRProj = archrprojPM, useMatrix = "PeakMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), maxCells = 200000, testMethod = "wilcoxon")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

write.table(markerList, file="markerpeak.markers.final.all.xls", quote=FALSE, sep="\t", row.names=FALSE)

heatmapPeaks <- plotMarkerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", transpose = TRUE)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "new-2M-Peak-Marker-Heatmap-FDR-0.01-FC-0.5.pdf", width = 28, height = 16, ArchRProj = archrprojPM, addDOC = FALSE)

#change4 if required
##Motif analysis #file downloaded:wget https://jeffgranja.s3.amazonaws.com/ArchR/Annotations/ArchR-Hg38-v1.Anno
archrprojPM <- addArchRAnnotations(ArchRProj = archrprojPM, db = "/rsrch3/home/genomic_med/ssatpati/Amit_Moran/Common_cells/new_module/ArchR-Hg38-v1.Anno", collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(
seMarker = markersPeaks,
ArchRProj = archrprojPM,
peakAnnotation = "EncodeTFBS",
cutOff = "FDR <= 0.01 & Log2FC >= 1"
)



heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 1000, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap.pdf", width = 20, height = 20, ArchRProj = archrprojPM, addDOC = FALSE)

##homer
archrprojPM <- addMotifAnnotations(ArchRProj = archrprojPM, motifSet = "homer", name = "Motif", force = TRUE)
motifs <- peakAnnoEnrichment(
seMarker = markersPeaks,
ArchRProj = archrprojPM,
peakAnnotation = "Motif",
cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

heatmapMotif <-  plotEnrichHeatmap(motifs, n = 1000, transpose = TRUE)
ComplexHeatmap::draw(heatmapMotif, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapMotif, name = "MotifTFBS-Enriched-Marker-Heatmap.pdf", width = 20, height = 20, ArchRProj = archrprojPM, addDOC = FALSE)


#Confusion Matrix
table(archrprojPM$Clusters)
cM <- confusionMatrix(paste0(archrprojPM$Clusters), paste0(archrprojPM$Sample))
cM


#Different colors: https://rdrr.io/github/GreenleafLab/ArchR/src/R/ColorPalettes.R
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whitePurple"),
    border_color = "black"
)

plotPDF(p, name = "Confusion-Matrix", addDOC = FALSE)

saveArchRProject(ArchRProj = archrprojPM, outputDirectory = "ArchR-Save-after-Motif", load = FALSE)

#Use below section for MarkerTest that for differential marker from differential cluster peaks. #didn’t worked for markerPeaks
