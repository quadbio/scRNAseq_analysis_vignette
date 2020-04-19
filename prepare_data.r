library(Matrix)
counts <- readMM("/links/groups/treutlein/USERS/zhisong_he/Data/cerebral_organoids_10X_Nature/human_cell_counts_GRCh38.mtx.gz")
meta <- read.csv("/links/groups/treutlein/USERS/zhisong_he/Data/cerebral_organoids_10X_Nature/metadata_human_cells.tsv.gz", stringsAsFactors=F, sep="\t")
genes <- read.csv("/links/groups/treutlein/USERS/zhisong_he/Data/cerebral_organoids_10X_Nature/genes_GRCh38.txt", stringsAsFactors=F, sep="\t", header=F)
genes[,3] <- "Gene Expression"
nfeatures <- colSums(counts>0)
mito <- colSums(counts[grep("^MT-", genes[,2]),]) / colSums(counts)

# DS1: organoid h409B2_60d_org1
idx_1 <- which(meta$Sample %in% c("h409B2_60d_org1") & (meta$in_LineComp | mito > 0.05 | nfeatures < 500 | nfeatures > 4000))
counts_1 <- counts[,idx_1]
barcodes <- meta$Barcode[idx_1]

## prepare data
writeMM(counts_1, file="/links/groups/treutlein/USERS/zhisong_he/Data/scRNAseq_analysis_vignette/data/DS1/matrix.mtx")
write.table(genes, file="/links/groups/treutlein/USERS/zhisong_he/Data/scRNAseq_analysis_vignette/data/DS1/features.tsv", row.names=F, col.names=F, quote=F, sep="\t")
write.table(barcodes, file="/links/groups/treutlein/USERS/zhisong_he/Data/scRNAseq_analysis_vignette/data/DS1/barcodes.tsv", row.names=F, col.names=F, quote=F)

## run Seurat and prepare figures
rownames(counts_1) <- make.names(genes[,2], unique=T)
colnames(counts_1) <- barcodes
seurat <- CreateSeuratObject(counts_1, project="DS1")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
png("images/vlnplot_QC_nopt.png", height=300, width=500)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
dev.off()
png("images/vlnplot_QC.png", height=300, width=500)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("images/scatterplot_QCmetrics.png", height=300, width=700)
plot1 + plot2
dev.off()

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)
seurat <- NormalizeData(seurat) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 50) %>% RunUMAP(dims = 1:20) %>% RunTSNE(dims = 1:20)
seurat <- FindNeighbors(seurat, dims = 1:20) %>% FindClusters(resolution = 1)

png("images/tsne_umap_featureplots.png", height=1200, width=800)
plot1 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"), ncol=3, reduction = "umap")
plot1 / plot2
dev.off()

png("images/variablefeatures.png", height=300, width=1000)
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2
dev.off()

png("images/elbowplot.png", height=300, width=500)
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))
dev.off()

png("images/pcheatmap.png", height=30, width=30, unit="cm", res=300)
PCHeatmap(seurat, dims = 1:20, cells = 500, nfeatures = 30, balanced = TRUE, ncol = 3, raster=F)
dev.off()

png("images/tsne_umap_nogroup.png", height=300, width=800)
plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2
dev.off()

png("images/tsne_umap_cluster.png", height=450, width=1000)
plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2
dev.off()

