
diff <- seurat_diff2(datafilt = CD14Mono,
                     group.by = "disease",
                     group1 = "SLE",
                     group2 = "HC",
                     assay = "RNA",
                     min.pct = 0.1,
                     thres.fc = 0)

CD14Mono_exp <- as.data.frame(t(CD14Mono_exp))

fea_sel

input <- data.frame(type = CD14Mono@meta.data$disease,
                    value = CD14Mono_exp$IFIT2)

ecdf_plot(input = input,
          ylim = c(0,1),
          title = "IFIT2")

i = "IFIT2"

my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')

#col <- colorRampPalette(c("#F0F0F0","#E5E0ED","#BFC6DD","#8CADCC","#4E92BA","#1871A8","#085889","#003758"))(100)

FeaturePlot(GSE135779_combined, features = i, reduction = "umap_harmony", pt.size = 1, raster = TRUE) + 
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

select <- fea_sel

my_colors <- c("#F0F0F0",'#EDD1D8', '#f4a3a8', '#e38191', '#cc607d', '#ad466c', '#8b3058', '#672044')

for (i in select) {
  # 创建绘图
  setwd("~/Steven Lijiajun/Steven/SLE/1. SLE CD14_Mono GSHstate/机器学习_113机器学习方案/113机器学习方案/10个umap")
  
  fig2 = FeaturePlot(GSE135779_combined, features = i, reduction = "umap_harmony",
                     pt.size = 1, raster = TRUE) + 
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
  file_name <- paste0(i, "umap_harmony.png")
  
  # 保存图形到 PDF 文件
  ggsave(file_name, fig2, height = 3.5, width = 4.7)
}

diff <- diff[diff$gene %in% select,]

## 分析开始

## Protein-encoding mRNA only

human_mrna <- read.table("/home/ug1268/tools/human_sc/mRNA_list.txt",sep="\t",header=T,check.names=F,row.names = 1)

common <- intersect(rownames(CD14Mono), rownames(human_mrna))
CD14Mono_mRNA <- CD14Mono[common,]

Idents(CD14Mono_mRNA) <- CD14Mono_mRNA$seurat_clusters


