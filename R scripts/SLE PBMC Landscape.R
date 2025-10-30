setwd("~/Steven Lijiajun/Steven/SLE")

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

# 公用基因文件
genes_file <- "GSE135779_genes.tsv.gz"

# 所有 barcode 和 matrix 文件
all_files <- list.files("GSE135779_RAW", full.names = TRUE)
barcode_files <- all_files[grep("_barcodes.tsv.gz$", all_files)]
matrix_files  <- all_files[grep("_matrix.mtx.gz$", all_files)]

# 提取样本名（去掉后缀）
sample_names <- gsub("_barcodes.tsv.gz", "", basename(barcode_files))

# 映射表：样本名 -> 分组（sample1） -> 四类标签（group）
sample_mapping <- data.frame(
  sample = c(
    "GSM4029896_JB17001", "GSM4029897_JB17002", "GSM4029898_JB17003", "GSM4029899_JB17004",
    "GSM4029900_JB17005", "GSM4029901_JB17006", "GSM4029902_JB17007", "GSM4029903_JB17008",
    "GSM4029904_JB17015", "GSM4029905_JB17016", "GSM4029906_JB17014", "GSM4029907_JB17010",
    "GSM4029908_JB17019", "GSM4029909_JB17020", "GSM4029910_JB17021", "GSM4029911_JB17022",
    "GSM4029912_JB17023", "GSM4029913_JB17024", "GSM4029914_JB17017", "GSM4029915_JB17018",
    "GSM4029916_JB18063", "GSM4029917_JB18064", "GSM4029918_JB18065", "GSM4029919_JB18066",
    "GSM4029920_JB18067", "GSM4029921_JB18068", "GSM4029922_JB18069", "GSM4029923_JB18070",
    "GSM4029924_JB18071", "GSM4029925_JB18072", "GSM4029926_JB18073", "GSM4029927_JB18074",
    "GSM4029928_JB18075", "GSM4029929_JB18076", "GSM4029930_JB18077", "GSM4029931_JB18078",
    "GSM4029932_JB18079", "GSM4029933_JB18080", "GSM4029934_JB18081", "GSM4029935_JB18082",
    "GSM4029936_JB18085", "GSM4029937_JB18086", "GSM4029938_JB18083", "GSM4029939_JB18084",
    "GSM4029940_JB19001", "GSM4029942_JB19003", "GSM4029943_JB19004", "GSM4029944_JB19006",
    "GSM4029945_JB19007", "GSM4029946_JB19008", "GSM4029947_JB19009", "GSM4029948_JB19010",
    "GSM4029949_JB19011", "GSM4029950_JB19013", "GSM4029951_JB19014", "GSM4029952_JB19015"
  ),
  sample1 = c(
    "cSLE1", "cSLE2", "cSLE3", "cSLE4", "cSLE5", "cSLE6", "cSLE7", "cSLE8",
    "cSLE10", "cSLE11", "cSLE9", "cHD1", "cSLE12", "cSLE13", "cSLE14", "cSLE15",
    "cSLE16", "cSLE17", "cHD2", "cHD3", "cSLE18", "cSLE19", "cSLE20", "cSLE21",
    "cSLE22", "cSLE23", "cHD4", "cHD5", "cSLE24", "cSLE25", "cSLE26", "cSLE27",
    "cSLE28", "cSLE29", "cHD6", "cHD7", "cSLE30", "cSLE31", "cSLE32", "cSLE33",
    "cHD10", "cHD11", "cHD8", "cHD9", "aHD1", "aHD3", "aSLE1", "aSLE2",
    "aSLE3", "aSLE4", "aHD4", "aHD5", "aHD6", "aSLE5", "aSLE6", "aSLE7"
  )
)
sample_mapping$group <- gsub("[0-9]+$", "", sample_mapping$sample1)
rownames(sample_mapping) <- sample_mapping$sample

# 初始化四个列表
cSLE_list <- list()
aSLE_list <- list()
cHD_list  <- list()
aHD_list  <- list()

# 读取并分类
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
  
  # 加入对应列表
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

# 合并每组为单个 Seurat 对象
cSLE <- merge(cSLE_list[[1]], y = cSLE_list[-1], add.cell.ids = names(cSLE_list), project = "cSLE")
aSLE <- merge(aSLE_list[[1]], y = aSLE_list[-1], add.cell.ids = names(aSLE_list), project = "aSLE")
cHD  <- merge(cHD_list[[1]],  y = cHD_list[-1],  add.cell.ids = names(cHD_list),  project = "cHD")
aHD  <- merge(aHD_list[[1]],  y = aHD_list[-1],  add.cell.ids = names(aHD_list),  project = "aHD")

# 成人/儿童分类
aCombined <- merge(aSLE, y = aHD, add.cell.ids = c("aSLE", "aHD"), project = "Adult")
cCombined <- merge(cSLE, y = cHD, add.cell.ids = c("cSLE", "cHD"), project = "Child")
Combined <- merge(cSLE, y = list(cHD, aHD, aSLE), add.cell.ids = c("cSLE", "cHD", "aHD", "aSLE"), project = "SLE_HC")

# 添加disease信息
aCombined$cell_origin <- ifelse(grepl("^aSLE_", colnames(aCombined)), "SLE", "HC")
aCombined@meta.data$disease <- aCombined$cell_origin

cCombined$cell_origin <- ifelse(grepl("^cSLE_", colnames(cCombined)), "SLE", "HC")
cCombined@meta.data$disease <- cCombined$cell_origin

Combined$cell_origin <- ifelse(grepl("^[ac]SLE_", colnames(Combined)), "SLE", "HC")
Combined$disease <- Combined$cell_origin

table(aCombined$disease)
table(cCombined$disease)
table(Combined$disease)

## Seurat Pipeline

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

#col <- colorRampPalette(c("#F0F0F0","#E5E0ED","#BFC6DD","#8CADCC","#4E92BA","#1871A8","#085889","#003758"))(100)

FeaturePlot(Combined, features = i, reduction = "umap_harmony", pt.size = 1, raster = TRUE) + 
  theme_bw() +  # Use a clean, white background
  theme(
    panel.grid = element_blank(),            # Remove grid lines
    axis.ticks = element_blank(),            # Remove axis ticks
    axis.text = element_blank(),             # Remove axis text
    axis.title = element_text(colour = "black", size = 15),  # Customize axis title
    plot.title = element_text(size = 17, hjust = 0.5),  # Customize plot title, center-align it
    panel.border = element_blank(),          # Remove the panel border
    legend.position = "right"                # Keep the color bar on the right
  ) + 
  scale_color_gradientn(colors = my_colors) +  # Apply the custom color palette for continuous data
  labs(x = ' ', y = ' ', title = i)

# 注释

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

## 组间描述性分析

# 背对背柱状图 SLE/HC
prop_back2back(datafilt = Combined_cleaned,
               group = "disease",
               cluster = "celltype_major",
               order = TRUE)

# 背对背柱状图 Pre/Post
prop_back2back(datafilt = Combined_cleaned,
               group = "sample_time",
               cluster = "celltype_major",
               order = TRUE)

# 棒棒糖图 R/NR
prop_back2back_lollipop(datafilt = Combined_cleaned,
                        group = "disease",
                        group1 = "HC",
                        group2 = "SLE",
                        cluster = "celltype_major")

# 棒棒糖图 Pre/Post
prop_back2back_lollipop(datafilt = Combined_cleaned,
                        group = "sample_time",
                        group1 = "post",
                        group2 = "pre",
                        cluster = "celltype_minor")

#### 出图代码

Fig1b <- dimplot_new(data = Combined_cleaned,
                     reduction = "umap_harmony",
                     pt.size = 1, label = F,
                     group.by = c("group"))

Fig1a <- dimplot_new(data = Combined_cleaned,
                     reduction = "umap_harmony",
                     pt.size = 1, label = F,
                     group.by = c("seurat_clusters"))

Fig1c <- dimplot_new(data = Combined_cleaned,
                     reduction = "umap_harmony",
                     pt.size = 1, label = F,
                     group.by = c("celltype_major"))

source('~/Steven Lijiajun/Steven/预后模型研究/BRCA_端粒+干性/scRNA-seq部分/图谱部分/ks_scAverExp.R')
source('~/Steven Lijiajun/Steven/预后模型研究/BRCA_端粒+干性/scRNA-seq部分/图谱部分/ks_scDotplot.R')

unique(Combined_cleaned@meta.data$celltype_major)

# dotplot
genes <- list("c01: CD4_T" = c("CD3D","CD3E","CCR7","IL7R"),
              "c02: nCD8_T" = c("CD8A","CD8B"),
              "c03: acCD8_T" = c("GZMA","GZMB","GZNH","GZMK"),
              "c04: NK" = c("KLRB1","NKG7"),
              'c05: B' = c("CD79A","MS4A1"),
              "c06: Plasma B" = c("XBP1","IGJ","MZB1"),
              "c07: CD14_Mono" = c("S100A8", "CD14","CST3","LGALS2","LYZ"),
              'c08: CD16_Mono' = c("FCER1G","FCGR3A"),
              "c09: cDC" = c("CLEC9A","ITGAX","ISG15"),
              "c10: pDC" = c("IFI6","IRF7","BCL11A"),
              "c11: Eryth" = c("HBB"))

featureSets <-genes 

Combined_cleaned@meta.data$celltype <- Combined_cleaned$celltype_major

Idents(Combined_cleaned) <- Combined_cleaned$celltype_major

Exp_scRNA <- ks_scAverExp(obj = Combined_cleaned,
                          features = featureSets,
                          CellTypes = unique(Combined_cleaned$celltype))


#rm(BC10,BC11,BC16,BC17,BC2,BC20,BC21,BC22,BC5,BC3,BC6)
#rm(scRNA.anchors,scRNAa,scRNA1)
#rm(scRNAlist)

#BiocManager::install("dittoSeq")
library(dittoSeq)

cols<- c("c01: CD4_T" = "#A6D719",
         "c02: nCD8_T" = "#176EBF",
         "c03: acCD8_T" = "#00A8DE",
         "c04: NK" = "#AEE0E8",
         'c05: B' = "#00A9A3",
         "c06: Plasma B" = "#FBD324",
         "c07: CD14_Mono" = "#F28A24",
         'c08: CD16_Mono' = "#A52828",
         "c09: cDC" = "#A37CB7",
         "c10: pDC" = "#F2D7EE",
         "c11: Eryth" = "#CD6981")

fig1d = scRNA_dotplot = ks_scDotplot(obj = Combined_cleaned,
                                  features = featureSets,
                                  CellTypes = unique(Combined_cleaned$celltype),
                                  avgPctMat= Exp_scRNA,
                                  pal = cols)
fig1d

# 背对背柱状图 SLE/HC
prop_back2back(datafilt = Combined_cleaned,
               group = "disease",
               cluster = "celltype_major",
               order = TRUE)

Combined_cleaned@meta.data$sample1 <- sample_mapping$sample1[match(
  Combined_cleaned@meta.data$sample,
  sample_mapping$sample
)]

table(Combined_cleaned$sample1)

dimplot_new(data = Combined_cleaned,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("sample1"))

select <- c("CD3D","CD4","CD8A","KLRB1","MS4A1","IGJ","CD14","FCER1G","ITGAX")
select2 <- c("GZMB","ISG15","HBB")

for (i in select2) {
  # 创建绘图
  setwd("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/Fig.1")
  
  my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')
  
  fig2 = FeaturePlot(Combined_cleaned, features = i, reduction = "umap_harmony",
                     pt.size = 1, raster=TRUE) + 
    theme_bw() +  # Use a clean, white background
    theme(
      panel.grid = element_blank(),            # Remove grid lines
      axis.ticks = element_blank(),            # Remove axis ticks
      axis.text = element_blank(),             # Remove axis text
      axis.title = element_text(colour = "black", size = 15),  # Customize axis title
      plot.title = element_text(size = 17, hjust = 0.5),  # Customize plot title, center-align it
      panel.border = element_blank(),          # Remove the panel border
      legend.position = "right"                # Keep the color bar on the right
    ) + 
    scale_color_gradientn(colors = my_colors) +  # Apply the custom color palette for continuous data
    labs(x = ' ', y = ' ', title = i)  # Place the gene label "CD3D" as the plot title
  
  # 构建文件名
  file_name <- paste0(i, "umap.png")
  
  # 保存图形到 PDF 文件
  ggsave(file_name, fig2, height = 3.5, width = 4.3)
}

seurat_circ_prop(datafilt = Combined_cleaned,
                 group = "disease",
                 cluster = "celltype",
                 xlab_order = NULL,
                 ncol = 2)
