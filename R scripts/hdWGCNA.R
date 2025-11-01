
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

## Recalculate GSH state

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
                        pt.size = 1e-18, 
                        color = "blue2red",
                        features = "Glutathione.metabolism",
                        raster = NULL,
                        outlier.rm = F)
plot

# Classic red-gray-green coloring

summary(CD14Mono_mRNA$Glutathione.metabolism)

CD14Mono_mRNA$GSH_state <- ifelse(
  CD14Mono_mRNA$Glutathione.metabolism < -0.01961, "c01: LGSH",
  ifelse(CD14Mono_mRNA$Glutathione.metabolism > 0.08705, "c03: HGSH", "c02: DTGSH")
)

table(CD14Mono_mRNA$GSH_state)

plot <- dimplot_new_cqw_DTYMK_nolegend_LDTH(data = CD14Mono_mRNA,
                    reduction = "umap_harmony",
                    pt.size = 1e-18, label = F,
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

# Setup for WGCNA
seurat_obj <- SetupForWGCNA(
  CD14Mono_mRNA,
  gene_select = "fraction",     # Select genes expressed in a fraction of cells
  fraction = 0.05,              # Expressed in at least 5% of cells
  wgcna_name = "tutorial"       # Set analysis name, results saved in @misc$tutorial
)

# Construct metacells by groups
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("GSH_state"),   # Grouping condition
  reduction = "umap_harmony",  # Dimensionality reduction space
  k = 25,                       # Number of neighbors per cell
  max_shared = 5,               # Maximum shared cells between metacells
  ident.group = "GSH_state"     # Set Idents for metacell Seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)

# Set expression matrix for a specific group
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "c03: HGSH",   # Target cell group
  group.by = 'GSH_state',      # Column used for grouping in meta.data
  assay = 'RNA',               # Assay to use
  slot = 'data'                # Data slot to use ("data" is log-normalized counts)
)

# Test soft-thresholding powers
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

# Visualize results
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol = 2)

# Construct network
seurat_obj <- ConstructNetwork(
  seurat_obj, 
  soft_power = 5,                 # Best soft power tested
  setDatExpr = FALSE,             
  corType = "pearson",            
  networkType = "signed",          
  TOMType = "signed",              
  detectCutHeight = 0.995,         # Initial dendrogram cut height
  minModuleSize = 50,              # Minimum module size
  mergeCutHeight = 0.2,            # Module merging threshold
  tom_outdir = "TOM",              # TOM matrix output path
  tom_name = "HGSH"                # TOM file prefix
)

PlotDendrogram(seurat_obj, main = "HGSH hdWGCNA Dendrogram")

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Calculate module eigengenes
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  scale.model.use = "linear",  
  assay = NULL,                
  pc_dim = 1                   
)

hMEs <- GetMEs(seurat_obj)  # Harmonized metacell expression
MEs  <- GetMEs(seurat_obj, harmonized = FALSE)  # Raw expression

# Compute module connectivity
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

# Plot module eigengenes
p <- PlotKMEs(
  seurat_obj,
  ncol = 5,
  n_hubs = 10,
  text_size = 2,
  plot_widths = c(3, 2)
)
p

# Module expression scores
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

# Module feature plots
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

# Module correlation plot
ModuleCorrelogram(seurat_obj,
                  exclude_grey = TRUE,
                  features = "hMEs")

MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']  # Exclude unassigned modules

seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# Dotplot of modules
p <- DotPlot(seurat_obj, features = mods, group.by = 'GSH_state') +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high = "#b11e24", mid = "grey95", low = "#077d7b") +
  scale_size(range = c(1, 8))  # Adjust point size range
p

# Module network plots
p1 <- ModuleNetworkPlot(
  seurat_obj,
  mods = "all",                
  outdir = "ModuleNetworks",   
  plot_size = c(6, 6),         
  label_center = FALSE,        
  edge.alpha = 0.25,           
  vertex.label.cex = 1,        
  vertex.size = 6              
)
p1

# Hub gene network
HubGeneNetworkPlot(
  seurat_obj,
  mods = "all",        
  n_hubs = 3,          
  n_other = 6,         
  edge_prop = 0.75     
)

# Run UMAP on modules
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,        
  n_neighbors = 15,   
  min_dist = 0.1      
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
  edge.alpha = 0.25,         
  sample_edges = TRUE,       
  edge_prop = 0.1,           
  label_hubs = 0,            
  keep_grey_edges = FALSE    
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

## Extract module genes
modules <- GetModules(seurat_obj)

table(modules$module)

genes_M1_M3 <- modules %>%
  filter(module %in% c("GSH_state-M1", "GSH_state-M3")) %>%
  pull(gene_name)

## HGSH vs LGSH differential expression

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

## Take intersection with module genes

common <- intersect(filtered_genes, genes_M1_M3)
common

write.table(common, file =  "common_diffHL_genesM1M3.txt", sep = "/t")

## Take intersection with upregulated bulk genes

deg_bulk <- diff_GSE65391 %>%
  filter(logfc > 0, pvalue < 0.05) %>%
  pull(id)

common_final <- intersect(common, deg_bulk)
common_final
