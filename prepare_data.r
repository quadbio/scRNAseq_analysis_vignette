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

ct_markers <- c("MKI67","NES","DCX","FOXG1","DLX2","DLX5","ISL1","SOX6","NKX2.1","EMX1","PAX6","GLI3","EOMES","NEUROD6","OTX2","LHX9","TFAP2A","HOXB2","HOXB5")
png("images/heatmap_ctmarkers.png", height=400, width=800)
DoHeatmap(seurat, features = ct_markers) + NoLegend()
dev.off()

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) %>% print(n = 50)

library(presto)
cl_markers_presto <- wilcoxauc(seurat)
cl_markers_presto %>% group_by(group) %>% top_n(n = 2, wt = logFC) %>% print(n = 40, width=Inf)

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
png("images/heatmap_clmarkers.png", height=1000, width=1200)
DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()
dev.off()

png("images/featureplot_vlnplot_examples.png", height=400, width=800)
plot1 <- FeaturePlot(seurat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(seurat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))
dev.off()

library(voxhunt)
load_aba_data('/links/groups/treutlein/PUBLIC_DATA/tools/voxhunt_rds')
genes_use <- variable_genes('E13', 300)$gene
vox_map <- voxel_map(seurat, genes_use=genes_use)
png("images/voxhunt.png", height=500, width=500)
plot_map(vox_map)
dev.off()

new_ident <- setNames(c("Dorsal telen. NPC","Midbrain-hindbrain boundary neuron","Dorsal telen. neuron","Dien. and midbrain excitatory neuron",
                        "Ventral telen. neuron","G2M dorsal telen. NPC","Dorsal telen. IP","Dien. and midbrain NPC",
                        "Dien. and midbrain IP and excitatory early neuron","G2M Dien. and midbrain NPC","G2M dorsal telen. NPC","Dien. and midbrain inhibitory neuron",
                        "Dien. and midbrain IP and early inhibitory neuron","Ventral telen. neuron","Unknown 1","Unknown 2"),
                      levels(seurat))
seurat <- RenameIdents(seurat, new_ident)
png("images/umap_annot.png", height=350, width=400)
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

seurat_dorsal <- subset(seurat, subset = RNA_snn_res.1 %in% c(0,2,5,6,10))
seurat_dorsal <- FindVariableFeatures(seurat_dorsal, nfeatures = 2000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), unlist(cc.genes))
seurat_dorsal <- CellCycleScoring(seurat_dorsal, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
seurat_dorsal <- ScaleData(seurat_dorsal)
seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
png("images/featureplot_examples_dorsal.png", height=300, width=1500)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)
dev.off()
seurat_dorsal <- ScaleData(seurat_dorsal, vars.to.regress = c("S.Score", "G2M.Score"))
seurat_dorsal <- RunPCA(seurat_dorsal) %>% RunUMAP(dims = 1:20)
png("images/featureplot_examples_dorsal_cc.png", height=300, width=1500)
FeaturePlot(seurat_dorsal, c("MKI67","GLI3","EOMES","NEUROD6"), ncol = 4)
dev.off()

library(destiny)
dm <- DiffusionMap(Embeddings(seurat_dorsal, "pca")[,1:20])
dpt <- DPT(dm)
seurat_dorsal$dpt <- rank(dpt$dpt)
png("images/dpt_featureplots_dorsal.png", height=300, width=1500)
FeaturePlot(seurat_dorsal, c("dpt","GLI3","EOMES","NEUROD6"), ncol=4)
dev.off()

png("images/dpt_scatterplots_dorsal.png", height=300, width=900)
plot1 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["GLI3",]), xlab="Dpt", ylab="Expression", main="GLI3") + geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot2 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["EOMES",]), xlab="Dpt", ylab="Expression", main="EOMES") + geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot3 <- qplot(seurat_dorsal$dpt, as.numeric(seurat_dorsal@assays$RNA@data["NEUROD6",]), xlab="Dpt", ylab="Expression", main="NEUROD6") + geom_smooth(se = FALSE, method = "loess") + theme_bw()
plot1 + plot2 + plot3
dev.off()

save(seurat, seurat_dorsal, file="../scRNAseq_analysis_vignette_seurat/DS1_seurat_objs.rdata")
saveRDS(seurat, file="../scRNAseq_analysis_vignette_seurat/DS1_seurat_obj.rds")
saveRDS(seurat_dorsal, file="../scRNAseq_analysis_vignette_seurat/DS1_seurat_obj_dorsal.rds")


