
# Set working directory
setwd("~/Steven Lijiajun/Steven/SLE")

# Load libraries
library(Seurat)
library(plyr)
library(scater)
library(stringr)
library(future)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggplot2)
library(data.table)
library(tidyverse)
library(future.apply)
library(harmony)
library(RColorBrewer)
library(SingleCellExperiment)

# Common gene file
genes_file <- "GSE135779_genes.tsv.gz"

# All barcode and matrix files
all_files <- list.files("GSE135779_RAW", full.names = TRUE)
barcode_files <- all_files[grep("_barcodes.tsv.gz$", all_files)]
matrix_files  <- all_files[grep("_matrix.mtx.gz$", all_files)]

# Extract sample names (remove suffix)
sample_names <- gsub("_barcodes.tsv.gz", "", basename(barcode_files))

# Mapping table: sample -> sample1 -> group
sample_mapping <- data.frame(
  sample = c(
    "GSM4029896_JB17001", "GSM4029897_JB17002", "GSM4029898_JB17003", "GSM4029899_JB17004",
    ...
  ),
  sample1 = c(
    "cSLE1", "cSLE2", "cSLE3", "cSLE4", ...
  )
)
sample_mapping$group <- gsub("[0-9]+$", "", sample_mapping$sample1)
rownames(sample_mapping) <- sample_mapping$sample

# Initialize four lists
cSLE_list <- list()
aSLE_list <- list()
cHD_list  <- list()
aHD_list  <- list()

# Read and classify samples
for (i in seq_along(sample_names)) {
  sample <- sample_names[i]
  group <- sample_mapping[sample, "group"]
  
  mat <- ReadMtx(
    mtx = matrix_files[i],
    cells = barcode_files[i],
    features = genes_file,
    feature.column = 2
  )
  obj <- CreateSeuratObject(counts = mat, project = sample)
  obj$sample <- sample
  obj$group <- group
  
  # Add to corresponding list
  if (group == "cSLE") {
    cSLE_list[[sample]] <- obj
  } else if (group == "aSLE") {
    aSLE_list[[sample]] <- obj
  } else if (group == "cHD") {
    cHD_list[[sample]] <- obj
  } else if (group == "aHD") {
    aHD_list[[sample]] <- obj
  }
}

# Merge each group into a single Seurat object
cSLE <- merge(cSLE_list[[1]], y = cSLE_list[-1], add.cell.ids = names(cSLE_list), project = "cSLE")
aSLE <- merge(aSLE_list[[1]], y = aSLE_list[-1], add.cell.ids = names(aSLE_list), project = "aSLE")
cHD  <- merge(cHD_list[[1]],  y = cHD_list[-1],  add.cell.ids = names(cHD_list),  project = "cHD")
aHD  <- merge(aHD_list[[1]],  y = aHD_list[-1],  add.cell.ids = names(aHD_list),  project = "aHD")

# Merge adult/child samples
aCombined <- merge(aSLE, y = aHD, add.cell.ids = c("aSLE", "aHD"), project = "Adult")
cCombined <- merge(cSLE, y = cHD, add.cell.ids = c("cSLE", "cHD"), project = "Child")
Combined <- merge(cSLE, y = list(cHD, aHD, aSLE), add.cell.ids = c("cSLE", "cHD", "aHD", "aSLE"), project = "SLE_HC")

# Add disease information
aCombined$cell_origin <- ifelse(grepl("^aSLE_", colnames(aCombined)), "SLE", "HC")
aCombined@meta.data$disease <- aCombined$cell_origin

cCombined$cell_origin <- ifelse(grepl("^cSLE_", colnames(cCombined)), "SLE", "HC")
cCombined@meta.data$disease <- cCombined$cell_origin

Combined$cell_origin <- ifelse(grepl("^[ac]SLE_", colnames(Combined)), "SLE", "HC")
Combined$disease <- Combined$cell_origin

table(aCombined$disease)
table(cCombined$disease)
table(Combined$disease)

## Seurat pipeline

nfeatures = 2000
ndim = 15
neigh = 75
dist = 0.5
res = 0.2

Combined <- NormalizeData(Combined, scale.factor = 10000,
                          normalization.method = "LogNormalize")

Combined <- FindVariableFeatures(Combined, nfeatures = nfeatures, 
                                 selection.method = "vst")

Combined <- ScaleData(Combined, features = VariableFeatures(Combined))

Combined <- RunPCA(Combined, assay = 'RNA', slot = 'scale.data')

Combined <- RunHarmony(Combined, group.by.vars = "sample",dims.use = 1:50,
                       assay.use = "RNA")

Combined <- FindNeighbors(Combined, k.param = neigh,
                          dims = 1:ndim, reduction = "harmony")

Combined <- FindClusters(Combined, resolution = 0.4, n.iter = 50)

Combined <- RunUMAP(Combined, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "harmony", reduction.name = "umap_harmony")

Combined <- RunTSNE(Combined, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "harmony", reduction.name = "tsne_harmony")

dimplot_new(data = Combined,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("seurat_clusters"))

dimplot_new(data = Combined,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("group"))

i = "CLEC9A"

my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')

FeaturePlot(Combined, features = i, reduction = "umap_harmony", pt.size = 1, raster = TRUE) + 
  theme_bw() +  
  theme(
    panel.grid = element_blank(),           
    axis.ticks = element_blank(),           
    axis.text = element_blank(),            
    axis.title = element_text(colour = "black", size = 15),  
    plot.title = element_text(size = 17, hjust = 0.5),  
    panel.border = element_blank(),         
    legend.position = "right"                
  ) + 
  scale_color_gradientn(colors = my_colors) +  
  labs(x = ' ', y = ' ', title = i)

# Annotation

Combined_cleaned <- subset(Combined, idents = c(0:9,11:13))

Combined_cleaned@meta.data <- Combined_cleaned@meta.data %>%
  mutate(celltype_major = recode(seurat_clusters,
                                 "0" = "c01: CD4_T",
                                 "1" = "c07: CD14_Mono",
                                 "2" = "c03: acCD8_T",
                                 "3" = "c02: nCD8_T",
                                 "4" = "c01: CD4_T",
                                 "5" = "c05: B",
                                 "6" = "c04: NK",
                                 "7" = "c08: CD16_Mono",
                                 "8" = "c11: Eryth",
                                 "9" = "c01: CD4_T",
                                 "10" = "c12: Platelet",
                                 "11" = "c09: cDC",
                                 "12" = "c10: pDC",
                                 "13" = "c06: Plasma B"))

Combined_cleaned$celltype_major <- factor(Combined_cleaned$celltype_major,
                                          levels = c(
                                            "c01: CD4_T",
                                            "c02: nCD8_T",
                                            "c03: acCD8_T",
                                            "c04: NK",
                                            "c05: B",
                                            "c06: Plasma B",
                                            "c07: CD14_Mono",
                                            "c08: CD16_Mono",
                                            "c09: cDC",
                                            "c10: pDC",
                                            "c11: Eryth"))

table(Combined_cleaned$celltype_major)

dimplot_new(data = Combined_cleaned,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("celltype_major"))
