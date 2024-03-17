# data download -----------------------------------------------------------
# https://www.brainimmuneatlas.org/download.php for GSE128855
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118948 for GSE118948

version

# load packages -----------------------------------------------------------
library(tidyverse)
library(Seurat)
library(parallel)
library(future)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


# save Rimage --------------------------------------------------------------
save.image(file = "Rproject_image.RData")

# GSE128855 Nat. Neurosci. data ------------------------

# read annotation matrix
annot_fullAggr <- read_csv("annot_fullAggr.csv")

# read 10X matrices
NN_agg_mtx <- ReadMtx(
  mtx="matrix.mtx",
  cells = "barcodes.tsv",
  features = "genes.tsv",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE)

# create seurat object
NN_agg_seurat <- CreateSeuratObject(counts = NN_agg_mtx)

# enable parallel
plan("multisession", workers = 4)
plan()

# add cell annotations provided by the original study
cell_type <- annot_fullAggr$cluster
names(cell_type) <- annot_fullAggr$cell
sample_location <- annot_fullAggr$sample
names(sample_location) <- annot_fullAggr$cell
NN_agg_seurat$cell_type <- cell_type
NN_agg_seurat$sample_location <- sample_location

# filter out low quality or dying cells: pct_mt >5, pct_ribo > 20, pct_Kcnq1ot1 > 0.01
NN_agg_seurat$percent_mt <- PercentageFeatureSet(NN_agg_seurat, pattern = "^mt-")
NN_agg_seurat$percent_ribo <- PercentageFeatureSet(NN_agg_seurat, pattern = "^Rp[ls]")
NN_agg_seurat$percent_Kcnq1ot1 <- PercentageFeatureSet(NN_agg_seurat, pattern = "Kcnq1ot1")
VlnPlot(NN_agg_seurat, c("percent_mt","percent_ribo", "percent_Kcnq1ot1"))
summary(NN_agg_seurat$percent_mt)
summary(NN_agg_seurat$percent_ribo)
summary(NN_agg_seurat$percent_Kcnq1ot1)
NN_agg_seurat<- subset(x = NN_agg_seurat,
                       subset = percent_mt < 5 & percent_ribo < 20 & percent_Kcnq1ot1 < 0.01)

# standard seurat pipeline
NN_agg_seurat <- NormalizeData(object = NN_agg_seurat, 
                               verbose = FALSE, 
                               scale.factor = 1e6) #ln(1+CPM)
NN_agg_seurat <- FindVariableFeatures(object = NN_agg_seurat,
                                      selection.method = "vst", 
                                      nfeatures = 2000, 
                                      verbose = FALSE)
NN_agg_seurat <- ScaleData(NN_agg_seurat)
NN_agg_seurat <- RunPCA(NN_agg_seurat)
# determine PCA
ElbowPlot(NN_agg_seurat,50)
# UMAP
NN_agg_seurat <- RunUMAP(NN_agg_seurat, reduction = "pca", dims = 1:5)

# condense cell identity names
Idents(NN_agg_seurat) <- 'cell_type'
levels(Idents(NN_agg_seurat))
NN_agg_seurat <- RenameIdents(object = NN_agg_seurat, 
                              "pDC" = "DCs",
                              "cDC2" = "DCs",
                              "migDC" = "DCs",
                              "cDC1" = "DCs",
                              "B cells" = "Lymphocytes",
                              "T/NKT cells" = "Lymphocytes",
                              "yd T cells" = "Lymphocytes",
                              "Neutorphils" = "Neutrophils",
                              "Monocytes/Mdc" = "Monocytes",
                              "BAM" = "BAMs",
                              "NA" = "NA"
)
NN_agg_seurat$cell_type_condense <- Idents(NN_agg_seurat)

# exclude NA cells
seurat_noNA <- subset(x = NN_agg_seurat,
                      subset = cell_type != c("NA"))
Idents(seurat_noNA) <- 'cell_type_condense'
seurat_noNA$cell_type_condense <- factor(seurat_noNA$cell_type_condense, levels = c("Microglia", "BAMs", "Monocytes", "Neutrophils", "Lymphocytes", "DCs", "NK cells"))
Idents(seurat_noNA) <- 'cell_type_condense'
# UMAP plot
DimPlot(seurat_noNA, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"))

# find CAMs markers over all populations
NN_DEGs_CAMs_all <- FindMarkers(
  subset(NN_agg_seurat,subset = cell_type!= c("NA")), 
  ident.1 = "CAMs", 
  min.pct = 0.4, 
  logfc.threshold = 4,
  only.pos = T) %>% arrange(desc(avg_log2FC))

# violin plot of CAMs specific markers
VlnPlot(subset(seurat_noNA, downsample = 200), c("Ms4a7"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Pf4"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Mrc1"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Igf1"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Dab2"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Cbr2"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Folr2"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

VlnPlot(subset(seurat_noNA, downsample = 200), c("Mgl2"), y.max = 10, cols = c("#43AA8B", "#f9844a","#277da1","#f94144","#f9c74f", "#90be6d","#f15bb5"), pt.size = 0.01, ncol = 1,adjust = 1,)+FontSize(x.text = 0, y.text = 6, x.title = 0, y.title = 0, main = 6)+NoLegend()

# calculate correlation matrix between monocytes, microglia, and CAMs
av.exp <- AverageExpression(
  subset(seurat_noNA, 
         idents = c("BAMs",  "Microglia",  "Monocytes")))$RNA
cor.exp <- as.data.frame(cor(av.exp))
cor.exp$x <- rownames(cor.exp)
write.csv(cor.exp, file = "3_cell_corr.csv")

# calculate correlation matrix between monocytes and CAMs at different locations
seurat_noNA$tissue_cell <- str_c(seurat_noNA$sample_location, "_", seurat_noNA$cell_type_condense)
unique(seurat_noNA$tissue_cell)
Idents(seurat_noNA) <- 'tissue_cell'
av.exp1 <- AverageExpression(subset(seurat_noNA, 
                                    idents = c("Dura mater_BAMs",  
                                               "Dura mater_Monocytes",
                                               "Enr. SDM_BAMs",
                                               "Enr. SDM_Monocytes",
                                               "Choroid plexus_BAMs",
                                               "Choroid plexus_Monocytes")))$RNA

cor.exp1 <- as.data.frame(cor(av.exp1))
cor.exp1$x <- rownames(cor.exp1)
write.csv(cor.exp1, file = "cell_location_corr.csv")


# GSE118948 Science data ------------------------------------------
# prepare raw data
Sci_agg_mtx <- list.files(pattern = "*_Naive_") %>% 
  lapply(.,read.table, sep = "\t", header = T, row.names = 1) %>% 
  lapply(.,CreateSeuratObject)

Sci_agg_seurat <- merge(Sci_agg_mtx[[1]], y = c(Sci_agg_mtx[[2]],Sci_agg_mtx[[3]],Sci_agg_mtx[[4]],Sci_agg_mtx[[5]],Sci_agg_mtx[[6]],Sci_agg_mtx[[7]],Sci_agg_mtx[[8]]), add.cell.ids = c("blood", "CP1", "CP2", "Men1", "Men2", "PV1", "PV2", "PV3"))

# qc
Sci_agg_seurat$percent_mt <- PercentageFeatureSet(Sci_agg_seurat, pattern = "^mt-")
Sci_agg_seurat$percent_ribo <- PercentageFeatureSet(Sci_agg_seurat, pattern = "^Rp[ls]")
Sci_agg_seurat$percent_Kcnq1ot1 <- PercentageFeatureSet(Sci_agg_seurat, pattern = "Kcnq1ot1")
VlnPlot(Sci_agg_seurat, features = c("percent_mt", "percent_ribo", "percent_Kcnq1ot1"))
Sci_agg_seurat<- subset(x = Sci_agg_seurat, subset = percent_mt < 5 & percent_ribo < 10 & percent_Kcnq1ot1 < 0.15)

# standard pipeline
Sci_agg_seurat <- NormalizeData(object = Sci_agg_seurat, verbose = FALSE, scale.factor = 1e6)
Sci_agg_seurat <- FindVariableFeatures(object = Sci_agg_seurat,
                                       selection.method = "vst", 
                                       nfeatures = 2000, 
                                       verbose = FALSE)
Sci_agg_seurat <- ScaleData(Sci_agg_seurat)
Sci_agg_seurat <- RunPCA(Sci_agg_seurat)
ElbowPlot(Sci_agg_seurat,20)
Sci_agg_seurat <- RunTSNE(Sci_agg_seurat, reduction = "pca", dims = 1:9, check_duplicates = FALSE)
DimPlot(Sci_agg_seurat, reduction = "tsne")
Sci_agg_seurat <- FindNeighbors(object = Sci_agg_seurat, dims = 1:4)
Sci_agg_seurat <- FindClusters(object = Sci_agg_seurat, resolution = 1) 
Idents(Sci_agg_seurat) <- "seurat_clusters"

DimPlot(Sci_agg_seurat, label = T, reduction = "tsne")

#plot microglia markers
FeaturePlot(Sci_agg_seurat, features = c("Tmem119", "Hexb", "P2ry12")) 
#plot monocyte markers
FeaturePlot(Sci_agg_seurat, features = c("Ccr2", "Cfp", "Cxcr4","Ly6c2","Cd300e", "Ace", "Treml4", "Spn"))
#plot BAMs markers
FeaturePlot(Sci_agg_seurat, features = c("Cd163", "Mrc1", "Ms4a7","Lyve1", "Pf4", "Mertk"))
#plot T lymphocyte markers
FeaturePlot(Sci_agg_seurat, features = c("Cd3d"))
#plot B lymphocyte markers
FeaturePlot(Sci_agg_seurat, features = c("Cd20"))
#plot DC cell markers
FeaturePlot(Sci_agg_seurat, features = c("Flt3", "Zbtb46", "Batf3", "Itgae", "Clec9a"))
#plot neutrophils markers
FeaturePlot(Sci_agg_seurat, features = "Ly6g")

# Name cell identity
Sci_agg_seurat <- RenameIdents(object = Sci_agg_seurat, 
                               `0` = "BAMs", 
                               `1` = "Microglia",
                               `2` = "Monocyte-derived",
                               `3` = "Microglia",
                               `4` = "Microglia",
                               `5` = "BAMs",
                               `6` = "Monocyte-derived",
                               `7` = "Granulocytes")

Sci_agg_seurat$cell_type <- Idents(Sci_agg_seurat)
Idents(Sci_agg_seurat) <- "cell_type"
DimPlot(Sci_agg_seurat, label = T, reduction = "tsne")

# Subset Microglia, CAMs, monocytes from Sci_agg_seurat
Sci_seurat_3 <- subset(Sci_agg_seurat, idents = c("BAMs","Microglia", "Monocyte") )
DimPlot(Sci_seurat_3)
Sci_seurat_3 <- NormalizeData(object = Sci_seurat_3, 
                              verbose = FALSE, 
                              scale.factor = 1e6)
Sci_seurat_3 <- FindVariableFeatures(object = Sci_seurat_3,
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = FALSE)
Sci_seurat_3 <- ScaleData(Sci_seurat_3)
Sci_seurat_3 <- RunPCA(Sci_seurat_3)
ElbowPlot(Sci_seurat_3,20)
Sci_seurat_3 <- RunTSNE(Sci_seurat_3, reduction = "pca", dims = 1:9)
Idents(Sci_seurat_3) <- "cell_type"
DimPlot(Sci_seurat_3, reduction = "tsne", 
        cols = c("#993399","#2d862d","#333399"),
        label.size = 7,
        pt.size = 1)

# subset monocytes from Sci_seurat_3
Idents(Sci_seurat_3) <- "cell_type"
Sci_seurat_mono <- subset(Sci_seurat_3, idents = c("Monocyte-derived"))
Sci_seurat_mono <- FindNeighbors(Sci_seurat_mono, dims = 1:11)
Sci_seurat_mono <- FindClusters(Sci_seurat_mono, resolution = 0.1)
Idents(Sci_seurat_mono) <- "seurat_clusters"
DimPlot(Sci_seurat_mono, 
        reduction = "tsne", 
        shape.by = "seurat_clusters",
        cols = c("#333399","#333399"),
        pt.size = 2,
        label.size = 7)
Idents(Sci_seurat_mono) <- "tissue_origin"
DimPlot(Sci_seurat_mono, 
        reduction = "tsne", 
        shape.by = "seurat_clusters",
        cols = c("#993333", "#336699", "#ffff00", "#339966"),
        pt.size = 2,
        label.size = 7)

# plot 2 types of monocytes
FeaturePlot(Sci_seurat_mono, "Ccr2", cols = c("grey","red"), reduction = "tsne")
FeaturePlot(Sci_seurat_mono, "Ly6c2", cols = c("grey","red"), reduction = "tsne")
FeaturePlot(Sci_seurat_mono, "Cx3cr1", cols = c("grey","red"), reduction = "tsne")
FeaturePlot(Sci_seurat_mono, "Sell", cols = c("grey","red"), reduction = "tsne")
FeaturePlot(Sci_seurat_mono, "Spn", cols = c("grey","red"), reduction = "tsne")
FeaturePlot(Sci_seurat_mono, "Treml4", cols = c("grey","red"), reduction = "tsne")

#Get gene counts in 2 monocyte subclusters
genelist <- c("Ccr2", "Ly6c2", "Cx3cr1", "Sell", "Spn", "Treml4")
Idents(Sci_seurat_mono) <- "seurat_clusters"
monocytes_2cluster_counts <- as.data.frame(Idents(Sci_seurat_mono))
for (i in genelist) {
  print(i)
  monocytes_2cluster_counts <- 
    cbind(monocytes_2cluster_counts,
          as.data.frame(Sci_seurat_mono@assays$RNA@data[i,]))
  print("done")
}
head(monocytes_2cluster_counts)
colnames(monocytes_2cluster_counts) <- c("Identity",genelist)
head(monocytes_2cluster_counts)
write.csv(monocytes_2cluster_counts, "monocytes_2cluster_normalized_count.csv")
