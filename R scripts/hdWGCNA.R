
## 删除杂群

table(CD14Mono_mRNA$GSH_state)

dimplot_new(data = CD14Mono_mRNA,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("seurat_clusters"))

CD14Mono_mRNA <- subset(CD14Mono_mRNA, idents = c(0:2))

nfeatures = 2000
ndim = 15
neigh = 75
dist = 0.5
res = 0.2

CD14Mono_mRNA <- FindNeighbors(CD14Mono_mRNA, k.param = neigh,
                               dims = 1:ndim, reduction = "harmony")

CD14Mono_mRNA <- FindClusters(CD14Mono_mRNA, resolution = 0.4, n.iter = 50)

CD14Mono_mRNA <- RunUMAP(CD14Mono_mRNA, dims = 1:ndim,
                         n.neighbors = neigh, min.dist = dist, 
                         reduction = "harmony", reduction.name = "umap_harmony")

dimplot_new(data = CD14Mono_mRNA,
            reduction = "umap_harmony",
            pt.size = 1, label = F,
            group.by = c("seurat_clusters"))

## 重新计算 GSH state

CD14Mono_mRNA@meta.data <- CD14Mono_mRNA@meta.data[,1:11]

path = "~/Steven Lijiajun/交大仁济肿瘤所/细胞系富集/all_genesets_gem.rds"
geneset1 = "metabolite"
geneset2 = "subsystem"

subsystem_score <- seurat_score(data = CD14Mono_mRNA,
                                source = path,
                                geneset = geneset2,
                                min.sz = 5)



CD14Mono_mRNA <- AddMetaData(CD14Mono_mRNA, subsystem_score)

CD14Mono_mRNA$Glutathione.metabolism

plot <- featureplot_new(data = CD14Mono_mRNA,
                        reduction = "umap_harmony",
                        pt.size = 0.000000000000000001, 
                        color = "blue2red",
                        features = "Glutathione.metabolism",
                        raster = NULL,
                        outlier.rm = F)
plot

# 经典红灰绿

summary(CD14Mono_mRNA$Glutathione.metabolism)

CD14Mono_mRNA$GSH_state <- ifelse(
  CD14Mono_mRNA$Glutathione.metabolism< -0.01961, "c01: LGSH",
  ifelse(CD14Mono_mRNA$Glutathione.metabolism > 0.08705, "c03: HGSH", "c02: DTGSH")
)

table(CD14Mono_mRNA$GSH_state)

plot <- dimplot_new_cqw_DTYMK_nolegend_LDTH(data = CD14Mono_mRNA,
                    reduction = "umap_harmony",
                    pt.size = 0.000000000000000001, label = F,
                    group.by = c("GSH_state"))
plot

## hdWGCNA pipeline

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)

Idents(CD14Mono_mRNA) <- CD14Mono_mRNA$GSH_state

seurat_obj <- SetupForWGCNA(
  CD14Mono_mRNA,
  gene_select = "fraction",     # 选择在一定比例细胞中表达的基因
  fraction = 0.05,              # 至少在 5% 细胞中表达
  wgcna_name = "tutorial"       # 设置分析名称，结果保存在 @misc$tutorial 下
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("GSH_state"),   # 分组条件
  reduction = "umap_harmony",                # 降维空间
  k = 25,                               # 每个细胞邻居数
  max_shared = 5,                      # 每对 metacell 最多共享细胞数
  ident.group = "GSH_state"              # 设置 metacell Seurat 对象的 Idents
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "c03: HGSH",   # 你感兴趣的细胞群（必须与 celltype 中一致）
  group.by = 'GSH_state',        # meta.data 中用于分组的列
  assay = 'RNA',                # 使用哪个 assay（通常为 RNA）
  slot = 'data'                 # 使用哪个数据层，"data" 是 log-normalized counts
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

# 可视化结果
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol = 2)

seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = 5,                   # 你已测试出的最佳软阈值
  setDatExpr = FALSE,              # 如果之前已经 SetDatExpr，则设为 FALSE
  corType = "pearson",             
  networkType = "signed",          
  TOMType = "signed",              
  detectCutHeight = 0.995,         # 用于初步切断树枝
  minModuleSize = 50,              # 模块最小基因数
  mergeCutHeight = 0.2,            # 模块合并阈值
  tom_outdir = "TOM",              # TOM 矩阵输出路径
  tom_name = "HGSH",                # TOM 文件名前缀
)

PlotDendrogram(seurat_obj, main = "HGSH hdWGCNA Dendrogram")

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use = "linear",  # 默认线性模型
  assay = NULL,                # 默认当前 default assay
  pc_dim = 1                   # 只取第一个主成分作为模块代表表达
)

hMEs <- GetMEs(seurat_obj)  # harmonized = TRUE（默认），整合 metacell 表达
MEs  <- GetMEs(seurat_obj, harmonized = FALSE)  # 原始表达，不做 metacell 合并

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'GSH_state',
  group_name = 'c03: HGSH',
  corFnc = "bicor",
  corOptions = "use='p'",
  harmonized = TRUE,
  assay = NULL,
  slot = "data"
)

seurat_obj <- ResetModuleNames(seurat_obj, new_name = "GSH_state-M")

modules <- GetModules(seurat_obj)
head(modules[, 1:6])

print(levels(modules$module))

p <- PlotKMEs(
  seurat_obj,
  ncol = 5,
  n_hubs = 10,
  text_size = 2,
  plot_widths = c(3, 2)
)
p

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method = "Seurat"
)

library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method = "UCell"
)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = "umap_harmony",
  features = "hMEs",
  order_points = TRUE,
  restrict_range = TRUE,
  point_size = 0.5,
  alpha = 1,
  label_legend = FALSE,
  raster_dpi = 500,
  raster_scale = 1,
  plot_ratio = 1,
  title = TRUE
)

wrap_plots(plot_list, ncol = 3)

ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE,
                  features = "hMEs")

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']  # 去除未分类模块

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(seurat_obj, features = mods, group.by = 'GSH_state')

p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high = "#b11e24", mid = "grey95", low = "#077d7b")
p

p <- DotPlot(seurat_obj, features = mods, group.by = 'GSH_state') +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high = "#b11e24", mid = "grey95", low = "#077d7b") +
  scale_size(range = c(1, 8))  # 调整点大小范围（最小值，最大值）
p

p1 <- ModuleNetworkPlot(
  seurat_obj,
  mods = "all",                # 绘制所有模块的网络图
  outdir = "ModuleNetworks",   # 网络图输出目录，将自动创建 PDF/PNG 文件
  plot_size = c(6, 6),         # 图的尺寸（宽，高，单位：英寸）
  label_center = FALSE,        # 不在图中心显示模块名
  edge.alpha = 0.25,           # 边透明度（越小越淡）
  vertex.label.cex = 1,        # 节点标签字体大小
  vertex.size = 6              # 节点圆圈大小
)
p1

HubGeneNetworkPlot(
  seurat_obj,
  mods = "all",        # 绘制所有模块的网络图
  n_hubs = 3,          # 每个模块选择 top 3 hub genes
  n_other = 6,         # 再随机选取 6 个非 hub gene 作为连接补充
  edge_prop = 0.75     # 显示边权值排名前 75% 的连接
)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,        # 每个模块取 top 10 hub genes
  n_neighbors = 15,   # UMAP 近邻参数
  min_dist = 0.1      # 控制 UMAP 点间距离稠密程度
)

umap_df <- GetModuleUMAP(seurat_obj)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    color = umap_df$color,
    size = umap_df$kME * 2
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.25,         # 边的透明度
  sample_edges = TRUE,       # 是否对边做子抽样（减少图形拥挤）
  edge_prop = 0.1,           # 随机保留前 10% 强度的边
  label_hubs = 0,            # 每个模块在图中标注前 2 个 hub genes
  keep_grey_edges = FALSE    # 是否保留 grey 模块的边，通常建议不保留
)

pdf("ModuleUMAPPlot.pdf", width = 10, height = 8)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha = 0.5,
  sample_edges = TRUE,
  edge_prop = 0.5,
  label_hubs = 0,
  keep_grey_edges = TRUE
)
dev.off()

## 提取模块基因
modules <- GetModules(seurat_obj)

table(modules$module)

genes_M1_M3 <- modules %>%
  filter(module %in% c("GSH_state-M1", "GSH_state-M3")) %>%
  pull(gene_name)


## HGSH vs LGSH

diff_HL <- seurat_diff2(datafilt = CD14Mono_mRNA,
                        group.by = "GSH_state",
                        group1 = "c03: HGSH",
                        group2 = "c01: LGSH",
                        assay = "RNA",
                        min.pct = 0.1,
                        thres.fc = 0.1)

filtered_genes <- diff_HL %>%
  filter(logfc > 0, pvalue < 0.05) %>%
  pull(gene)

## 取交集

common <- intersect(filtered_genes, genes_M1_M3)

common

write.table(common, file =  "common_diffHL_genesM1M3.txt", sep = "/t")

## 和bulk里面上调的基因取交集

deg_bulk <- diff_GSE65391 %>%
  filter(logfc > 0, pvalue < 0.05) %>%
  pull(id)

common_final <- intersect(common, deg_bulk)

common_final

