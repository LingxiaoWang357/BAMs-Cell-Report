# random seed -------------------------------------------------------------
set.seed(123)

# load library ------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony)
library(SCORPIUS)

# save image --------------------------------------------------------------
save.image(file = "BAMs_project_image.RData")


# load data ---------------------------------------------------------------

# Ctrl sample
Ctrl1 <- ReadMtx(
  mtx="../GSM4525522_Ctrl-1_matrix.mtx.gz",
  cells = "../GSM4525522_Ctrl-1_barcodes.tsv.gz",
  features = "../GSM4525522_Ctrl-1_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl1 <- CreateSeuratObject(counts = Ctrl1, min.cells = 10, min.features = 200)
Ctrl1[["Condition"]] <-  c('Ctrl')

Ctrl1[["sample_name"]] <-  c('Ctrl-1')
Ctrl1[["percent_mt"]] <- PercentageFeatureSet(Ctrl1, pattern = "^mt-")

Ctrl2 <- ReadMtx(
  mtx="../GSM4525523_Ctrl-2_matrix.mtx.gz",
  cells = "../GSM4525523_Ctrl-2_barcodes.tsv.gz",
  features = "../GSM4525523_Ctrl-2_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl2 <- CreateSeuratObject(counts = Ctrl2, min.cells = 10, min.features = 200)
Ctrl2[["Condition"]] <-  c('Ctrl')

Ctrl2[["sample_name"]] <-  c('Ctrl-2')
Ctrl2[["percent_mt"]] <- PercentageFeatureSet(Ctrl2, pattern = "^mt-")


Ctrl3 <- ReadMtx(
  mtx="../GSM4525524_Ctrl-3_matrix.mtx.gz",
  cells = "../GSM4525524_Ctrl-3_barcodes.tsv.gz",
  features = "../GSM4525524_Ctrl-3_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
Ctrl3 <- CreateSeuratObject(counts = Ctrl3, min.cells = 10, min.features = 200)
Ctrl3[["Condition"]] <-  c('Ctrl')

Ctrl3[["sample_name"]] <-  c('Ctrl-3')
Ctrl3[["percent_mt"]] <- PercentageFeatureSet(Ctrl3, pattern = "^mt-")

#Day0 sample
D0_1 <- ReadMtx(
  mtx="../GSM4525525_D0-1_matrix.mtx.gz",
  cells = "../GSM4525525_D0-1_barcodes.tsv.gz",
  features = "../GSM4525525_D0-1_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_1 <- CreateSeuratObject(counts = D0_1, min.cells = 10, min.features = 200)
D0_1[["Condition"]] <-  c('D0')

D0_1[["sample_name"]] <-  c('D0-1')
D0_1[["percent_mt"]] <- PercentageFeatureSet(D0_1, pattern = "^mt-")


D0_2 <- ReadMtx(
  mtx="../GSM4525526_D0-2_matrix.mtx.gz",
  cells = "../GSM4525526_D0-2_barcodes.tsv.gz",
  features = "../GSM4525526_D0-2_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_2 <- CreateSeuratObject(counts = D0_2, min.cells = 10, min.features = 200)
D0_2[["Condition"]] <-  c('D0')

D0_2[["sample_name"]] <-  c('D0-2')
D0_2[["percent_mt"]] <- PercentageFeatureSet(D0_2, pattern = "^mt-")


D0_3 <- ReadMtx(
  mtx="../GSM4525527_D0-3_matrix.mtx.gz",
  cells = "../GSM4525527_D0-3_barcodes.tsv.gz",
  features = "../GSM4525527_D0-3_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D0_3 <- CreateSeuratObject(counts = D0_3, min.cells = 10, min.features = 200)
D0_3[["Condition"]] <-  c('D0')

D0_3[["sample_name"]] <-  c('D0-3')
D0_3[["percent_mt"]] <- PercentageFeatureSet(D0_3, pattern = "^mt-")

# Day2 sample
D2_1 <- ReadMtx(
  mtx="../GSM4525528_D2-1_matrix.mtx.gz",
  cells = "../GSM4525528_D2-1_barcodes.tsv.gz",
  features = "../GSM4525528_D2-1_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D2_1 <- CreateSeuratObject(counts = D2_1, min.cells = 10, min.features = 200)
D2_1[["Condition"]] <-  c('Repop_D2')

D2_1[["sample_name"]] <-  c('D2-1')
D2_1[["percent_mt"]] <- PercentageFeatureSet(D2_1, pattern = "^mt-")


D2_2 <- ReadMtx(
  mtx="../GSM4525529_D2-2_matrix.mtx.gz",
  cells = "../GSM4525529_D2-2_barcodes.tsv.gz",
  features = "../GSM4525529_D2-2_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D2_2 <- CreateSeuratObject(counts = D2_2, min.cells = 10, min.features = 200)
D2_2[["Condition"]] <-  c('Repop_D2')

D2_2[["sample_name"]] <-  c('D2-2')
D2_2[["percent_mt"]] <- PercentageFeatureSet(D2_2, pattern = "^mt-")


D2_3 <- ReadMtx(
  mtx="../GSM4525530_D2-3_matrix.mtx.gz",
  cells = "../GSM4525530_D2-3_barcodes.tsv.gz",
  features = "../GSM4525530_D2-3_features.tsv.gz",
  cell.column = 1,
  feature.column = 2,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
D2_3 <- CreateSeuratObject(counts = D2_3, min.cells = 10, min.features = 200)
D2_3[["Condition"]] <-  c('Repop_D2')

D2_3[["sample_name"]] <-  c('D2-3')
D2_3[["percent_mt"]] <- PercentageFeatureSet(D2_3, pattern = "^mt-")

# data integration --------------------------------------------------------
seurat_all <- merge(Ctrl1, y = c(Ctrl2, Ctrl3, D2_1, D2_2, D2_3, D0_1, D0_2, D0_3))
VlnPlot(seurat_all, "percent_mt")
seurat_all <- subset(seurat_all, percent_mt < 10)
VlnPlot(seurat_all, "percent_ribo")
seurat_all <- subset(seurat_all, percent_ribo < 30)
VlnPlot(seurat_all, "nFeature_RNA")
seurat_all <- subset(seurat_all, nFeature_RNA > 1000)
VlnPlot(seurat_all, "nCount_RNA")
seurat_all <- subset(seurat_all, subset = nCount_RNA > 200 & nCount_RNA < 50000)
seurat_all <- seurat_all %>%
  NormalizeData(.,scale.factor = 1e6) %>% 
  ScaleData(.) %>%
  FindVariableFeatures(.) %>% 
  RunPCA(.)

# harmony
seurat_all <- seurat_all %>% 
  RunHarmony(c("sample_name"), 
             plot_convergence = F, 
             dims.use = 1:20)

seurat_all <- RunUMAP(seurat_all, 
                      reduction = "harmony", 
                      dims = 1:50)

Idents(seurat_all) <- "sample_name"
DimPlot(seurat_all, split.by = "Condition")
DimPlot(seurat_all)

seurat_all <- FindNeighbors(seurat_all,reduction = "harmony", dims = 1:50)
seurat_all <- FindClusters(seurat_all, resolution = c(0.3,0.4,0.5))

Idents(seurat_all) <- "RNA_snn_res.0.3"
DimPlot(seurat_all, label = T)
DimPlot(seurat_all, label = T, split.by = "Condition")

seurat_all <- RenameIdents(seurat_all, 
                           "0"="Microglia",
                           "1"="Microglia",
                           "2"="Microglia",
                           "3"="Microglia",
                           "4"="Microglia",
                           "5"="BAMs/MdCs/DCs",
                           "6"="Neutrophils",
                           "7"="BAMs/MdCs/DCs",
                           "8"="Microglia",
                           "9"="Endothelia",
                           "10"="Endothelia",
                           "11"="Astrocytes",
                           "12"="NK/NKT",
                           "13"= "Choroid Plexus")
FeaturePlot(seurat_all, c("Mrc1","Ttr"))

seurat_all$cell_type <- Idents(seurat_all)
DimPlot(seurat_all, label = T, cols = c('#f3722c', '#F94144','#43AA8B', '#90be6d','#f9c74f','#577590', '#f8961e'))
DimPlot(seurat_all, cols = c('#f3722c', '#F94144','#43AA8B', '#90be6d','#f9c74f','#577590', '#f8961e'))
DimPlot(seurat_all, label = T, split.by = "Condition")

DoHeatmap(subset(seurat_all, downsample = 50),
          group.colors = c('#f3722c', '#F94144','#43AA8B', '#90be6d','#f9c74f','#577590', '#f8961e'),
          features = as.vector(c(
            "Cx3cr1","Hexb","Tmem119","P2ry12","Sall1", #microglia
            "Apoe","Ms4a6c","Lyz2","Dab2","Ms4a7","Mrc1","Pf4","Mki67","Lyve1",  # BAMs/cMdC/DCs
            
            "Flt3","H2-Aa", "H2-Ab1", "H2-Eb1", #DC
            "Ccr2","Ly6c2","Clec4e", #monocytes
            "Ly6g", #Neutrophils
            "Cd79a", "Ms4a1", # B cell
            "Cd3e", "Trac", "Trbc2", #T cell
            "Cldn5", "Itm2a", "Adgrl4", "Vtn", "Acta2", #endothelium,
            "Aldh1l1", "Slc1a2","Aqp4", #Astrocytes
            "Klrb1c", "Ncr1", "Gzma", # NK cell
            "Ttr"#Choroid Plexus
            
          )),
          slot = "data",
          assay = "RNA",
          raster = T,lines.width = 2) + 
  scale_fill_gradientn(colors = c("lightgrey","#e60000"),na.value = "white") +
  theme(plot.margin = unit(c(3,2,1,1), "cm"),
        legend.position = "bottom", 
        text = element_text(size=8))

# Analyze monocytes n BAMs populations -------------------------------
seurat_subset <- subset(seurat_all, idents = c("BAMs/MdCs/DCs"))

seurat_subset <- seurat_subset %>% 
  RunHarmony(c("sample_name"), plot_convergence = F)

seurat_subset <- RunUMAP(seurat_subset, reduction = "harmony", dims = 1:50)
seurat_subset <- FindNeighbors(seurat_subset,reduction = "harmony", dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, 
                              resolution = c(0.1, 0.15, 0.2, 0.25, 0.3))

Idents(seurat_subset) <- "RNA_snn_res.0.15"
DimPlot(seurat_subset)

seurat_subset <- subset(seurat_subset, subset = RNA_snn_res.0.15 %in% c("5","7","8"), invert = T)

seurat_subset <- RunUMAP(seurat_subset, reduction = "harmony", dims = 1:50)
seurat_subset <- FindNeighbors(seurat_subset,reduction = "harmony", dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, 
                              resolution = c(0.1, 0.15, 0.2, 0.25, 0.3))

Idents(seurat_subset) <- "RNA_snn_res.0.15"
DimPlot(seurat_subset)


seurat_subset <- RenameIdents(seurat_subset, 
                              "0"="ncMonocytes",
                              "1"="DCs",
                              "2"="cMdC/MHCII-hi BAMs",
                              "3"="cMonocytes",
                              "4" = "MHCII-lo BAMs")
seurat_subset$cell_type <- Idents(seurat_subset)

levels(seurat_subset) <- c("ncMonocytes", "cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs", "DCs")

DimPlot(seurat_subset, cols = c('#0f4c5c', '#fb8b24', '#cc5a12','#9a031e', '#330822'))

DimPlot(subset(seurat_subset, 
               subset=cell_type %in% 
                 c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs")), 
        cols = c('#fb8b24', '#cc5a12','#9a031e'))

DimPlot(subset(seurat_subset, 
               subset=cell_type %in% 
                 c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs")), 
        cols = c('#fb8b24', '#cc5a12','#9a031e'), split.by = "Condition")


DimPlot(seurat_subset, split.by = "Condition", cols = c('#0f4c5c', '#fb8b24', '#cc5a12','#9a031e', '#330822'))

DoHeatmap(subset(seurat_subset, downsample = 50),
          group.colors = c('#0f4c5c', '#490b31','#fb8b24', '#cc5a12','#9a031e'),
          features = as.vector(c(
            "Cx3cr1", "Spn", "Treml4", #ncMonocytes
            "Flt3","Ccr7", "Itgax", "Siglech", "H2-Aa", "H2-Ab1", "H2-Eb1",
            "Ccr2","Ly6c2","Sell", #monocytes#DC
            "Dab2","Ms4a7","Msr1","Pf4","Lyve1", "Folr2", "Clec4n", "Nrp1", "Mrc1"#BAMs
          )),
          slot = "data",
          assay = "RNA",
          raster = T,lines.width = 2) + 
  scale_fill_gradientn(colors = c("lightgrey","#e60000"),na.value = "white") +
  theme(plot.margin = unit(c(3,2,1,1), "cm"),
        legend.position = "bottom", 
        text = element_text(size=8))

FeaturePlot(seurat_subset, 
            c("Csf1r"), split.by = "Condition", 
            min.cutoff = "q20",
            max.cutoff = "q98", cols = c("grey", "red"))
FeaturePlot(seurat_subset, c("Ccr2"), split.by = "Condition", 
            #cells = c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs"),
            min.cutoff = "q10",
            max.cutoff = "q90",
            cols = c("grey", "red"))
FeaturePlot(subset(seurat_subset, 
                   subset=cell_type %in% 
                     c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs")),
                   c("Ccr2"), 
                   split.by = "Condition", 
            min.cutoff = "q1",
            max.cutoff = "q85",
            cols = c("grey", "red"))
FeaturePlot(seurat_subset, c("H2-Ab1"), split.by = "Condition",
            min.cutoff = "q20",
            max.cutoff = "q90",
            cols = c("grey", "red"))
FeaturePlot(subset(seurat_subset, 
                   subset=cell_type %in% 
                     c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs")),
            c("H2-Ab1"), 
            split.by = "Condition", 
            min.cutoff = "q20",
            max.cutoff = "q90",
            cols = c("grey", "red"))
FeaturePlot(subset(seurat_subset, 
                   subset=cell_type %in% 
                     c("cMonocytes","cMdC/MHCII-hi BAMs","MHCII-lo BAMs")),
            c("Mrc1"), 
            split.by = "Condition", 
            min.cutoff = "q20",
            max.cutoff = "q90",
            cols = c("grey", "red"))



VlnPlot(seurat_subset, c("H2-Ab1"), split.by = "Condition",
        idents = c("cMdC/MHCII-hi BAMs"),cols = c('#cc5a12','#cc5a12','#cc5a12'))

VlnPlot(seurat_subset, c("Ccr2"), split.by = "Condition",
        idents = c("cMdC/MHCII-hi BAMs"),cols = c('#cc5a12','#cc5a12','#cc5a12'))

FeaturePlot(seurat_subset, c("H2-Aa"), split.by = "Condition", min.cutoff = "q20")

DEGs_cMdC <- FindMarkers(seurat_subset,  
                         subset.ident = "2", 
                         group.by = "Condition",
                         ident.1 = "Repop_D2",
                         ident.2 = "D0")

# SCOPIUS -----------------------------------------------------------------
expression <- t(as.matrix(subset(seurat_subset, 
                                 subset = cell_type %in% c("cMonocytes", 
                                                           "cMdC/MHCII-hi BAMs",
                                                           "MHCII-lo BAMs"))@assays$RNA@data))

group_name <- subset(seurat_subset, 
                     subset = cell_type %in% c("cMonocytes", 
                                               "cMdC/MHCII-hi BAMs",
                                               "MHCII-lo BAMs"))@meta.data$cell_type 
space <- reduce_dimensionality(expression, dist = "pearson", 
                               ndim = 5)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space,k = 4)
traj_color <- c('#fb8b24', '#cc5a12','#9a031e')
names(traj_color) <- c("cMonocytes", 
                       "cMdC/MHCII-hi BAMs",
                       "MHCII-lo BAMs")
draw_trajectory_plot(
  space, 
  progression_group = group_name,
  path = traj$path,
  contour = T,
  progression_group_palette = traj_color)


expression <- as.matrix(expression)
gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:250,]
expr_sel <- expression[,gene_sel$gene]
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)

draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules,
                        show_labels_row = T,
                        fontsize_row = 5,
                        progression_group_palette = traj_color)

write.csv(modules, file = "SCORPIUS_modules.csv")


# get expression data -----------------------------------------------------
Ccr2_expression <- as.data.frame(subset(seurat_subset, 
                                        subset=cell_type %in% c("cMdC/MHCII-hi BAMs"))@assays$RNA@data[c("Ccr2"),])
H2_Ab1_expression <- as.data.frame(subset(seurat_subset, 
                                          subset=cell_type %in% c("cMdC/MHCII-hi BAMs"))@assays$RNA@data[c("H2-Ab1"),])

Ccr2_expression <- cbind(rownames(Ccr2_expression), Ccr2_expression, H2_Ab1_expression)
colnames(Ccr2_expression) <- c("CellID","Ccr2", "H2-Ab1")

Ccr2_expression_name <- as.data.frame(cbind(subset(seurat_subset, 
                                                   subset=cell_type %in% 
                                                     c("cMdC/MHCII-hi BAMs"))@meta.data$Condition, 
                                            subset(seurat_subset, 
                                                   subset=cell_type %in% 
                                                     c("cMdC/MHCII-hi BAMs"))@meta.data$sample_name,
                                            subset(seurat_subset,
                                                   subset=cell_type %in% 
                                                     c("cMdC/MHCII-hi BAMs"))@assays$RNA@data@Dimnames[[2]]))

colnames(Ccr2_expression_name) <- c("Condition","SampleName", "CellID")

write.csv(left_join(Ccr2_expression, Ccr2_expression_name,"CellID"), "cMdCs_Ccr2_H2Ab1_Expression.csv")


