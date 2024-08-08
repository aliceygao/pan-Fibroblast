library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(ggrastr)
library(ggpubr)
library(ggpubr)
library(ggrepel)
library(grDevices)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(data.table)
library(ggthemes)
library(RColorBrewer)
library(pheatmap)

# ------------- Figure 1 --------------
# UMAP of fibroblast clusters
scPalette <- function(n) {
  colorSpace <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F29403", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#222F75", "#1B9E77", "#B2DF8A",
    "#E3BE00", "#FB9A99", "#E7298A", "#910241", "#00CDD1", "#A6CEE3", "#CE1261", "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", "#B3DE69", "#8DD3C7", "#999999"
  )
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

seu <- readRDS("./data/all_fib_cells.rds")
metaData <- seu@meta.data
plotData <- metaData %>%
  dplyr::select(umap_1, umap_2, cluster) %>%
  dplyr::mutate(cellClusters = cluster)

IdentLabel <- plotData %>%
  group_by(cellClusters) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

IdentLabel$Label <- IdentLabel$cellClusters

colN <- 2
cellColors <- scPalette(nrow(IdentLabel))
names(cellColors) <- IdentLabel$cellClusters

f1 <- ggplot(plotData, aes(x = umap_1, y = umap_2, colour = cellClusters)) +
  geom_point_rast(size = 0.15, stroke = 0, shape = 16, raster.dpi = 600) +
  geom_point(data = IdentLabel, color = cellColors, fill = "white", size = 7, alpha = 0.6, shape = 21) +
  geom_text(data = IdentLabel, aes(label = Label), colour = "black", size = 2.5) +
  scale_colour_manual(labels = IdentLabel$cellClusters, values = cellColors) +
  guides(colour = guide_legend(
    ncol = colN, override.aes = list(size = 3),
    label.theme = element_text(size = 6),
    keyheight = 0.5
  )) +
  labs(x = "umap_1", y = "umap_2", title = "") +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.box.background = element_blank(),
    legend.position = "right",
    legend.spacing.x = unit(0, "cm")
  )
ggsave(sprintf("%s/Figure1A.pdf", figdir), f1, width = 5, height = 4)


# Violin plot for fibroblast signatures
load("./data/sigVar.RData")
sigScore <- readRDS("./data/all_sigscore.rds")
sigScore_long <- as.data.frame(sigScore) %>%
  rownames_to_column(var = "sig") %>%
  pivot_longer(!sig, names_to = "cellName", values_to = "score")
sigScore_meta <- left_join(sigScore_long, metaData, by = "cellName")

sigScore_meta_sub <- sigScore_meta %>% filter(sig %in% sigVar)
sigScore_meta_sub$celltype <- factor(sigScore_meta_sub$celltype, levels = cluster.order)

cluster.color <- scPalette(20)
names(cluster.color) <- sort(unique(metaData$cluster))
g <- ggplot(sigScore_meta %>% filter(sig %in% sigVar), aes(x = cluster, y = score, fill = cluster, color = cluster)) +
  geom_violin(
    scale = "width",
    adjust = 1,
    trim = TRUE
  ) +
  facet_wrap(~sig, nrow = 3, scales = "free_y", strip.position = "top", dir = "h") +
  theme_classic() +
  ylab("Gene score") +
  xlab("") +
  scale_color_manual(values = cluster.color) +
  scale_fill_manual(values = cluster.color) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "none", strip.text = element_text(size = 10)
  )
ggsave(sprintf("%s/Figure1D.pdf", figdir), g, width = 7, height = 7)


# Heatmap for top marker genes
load("./data/fib_DEgenes.RData")
load("./data/avgExp_DEgenes.RData")
load("./data/clusterOrder.RData")

expr <- avgExp[topGenes, cluster.order]
mat2 <- t(scale(t(expr)))
mat2 <- MinMax(mat2, max = 2.5, min = (-1) * 2.5)
colnames(mat2) <- colnames(expr)
plotData <- mat2

col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))

cols <- col.heat(100)
breaks <- seq(-2.5, 2.5, length = 100)
col1 <- colorRamp2(breaks, cols)

ha <- HeatmapAnnotation(
  FCs = cluster.order,
  col = list(FCs = cluster.color)
)

markerGeneL <- rownames(plotData) %in% labelGenes
markerGenes <- rownames(plotData)[markerGeneL]

rowHa <- rowAnnotation(
  link = anno_mark(
    at = which(markerGeneL),
    labels = markerGenes, labels_gp = gpar(fontsize = 10), padding = 0.5
  ),
  width = unit(1, "cm") + max_text_width(markerGenes, gp = gpar(fontsize = 12))
)

pdf(sprintf("%s/Figure1E.pdf", figdir), width = 8, height = 12)
Heatmap(plotData,
  name = "z-score",
  top_annotation = ha,
  right_annotation = rowHa,
  cluster_rows = F, cluster_columns = F,
  col = col1, border = F,
  show_row_names = FALSE, show_column_names = TRUE
)
dev.off()


# Dot plot of fibroblast abundance across tissues
metadata_summary <- metaData %>%
  group_by(tissue, cluster) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100)

metadata_summary$tissue <- factor(metadata_summary$tissue, levels = c("Synovium", "Colon", "Stomach", "Breast", "Kidney", "Lung", "Pancreas", "Ovary", "Liver", "Skin", "HNSC"))
metadata_summary$cluster <- factor(metadata_summary$cluster, levels = cluster.order)

mycols <- brewer.pal(9, "GnBu")[c(1, 3, 5, 7)]
names(mycols) <- c("<1%", "1%-5%", "5%-10%", ">10%")

p <- ggplot(data = metadata_summary, mapping = aes_string(x = "tissue", y = "cluster")) +
  geom_point_rast(mapping = aes_string(size = "clust_prop", color = "label")) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)
  ) +
  scale_color_manual(values = mycols)
ggsave(sprintf("%s/Figure1F_dot.pdf", figdir), p, width = 5, height = 8)

metadata_summary <- metadata_summary %>%
  mutate(label = ifelse(clust_prop > 10, ">10%", ifelse(clust_prop > 5, "5%-10%", ifelse(clust_prop > 1, "1%-5%", "<1%"))))

gg.ls <- list()
for (i in 1:length(cluster.order)) {
  count.sub <- metadata_summary %>% filter(cluster == rev(cluster.order)[i])
  count.sub <- count.sub %>%
    group_by(label) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(clust_total = sum(count)) %>%
    dplyr::mutate(clust_prop = count / clust_total * 100)

  mycols <- c("#D8D9DA", "#ECB884", "#E4E45F", "#4758A2")
  names(mycols) <- c("<1%", "1%-5%", "5%-10%", ">10%")
  mycols <- mycols[count.sub$label]
  gg.ls[[i]] <- ggplot(count.sub, aes(x = "", y = clust_prop, fill = label)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = mycols) +
    theme_void()
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 20, ncol = 1)
ggexport(multi.page, filename = sprintf("%s/Figure1F_pie.pdf", figdir), width = 1, height = 20)


# Phenotype enrichment of fibroblast clusters
metaData$pcat <- factor(metaData$pcat, levels = c("Healthy", "Non-inflamed", "Inflamed", "Normal", "Fetal", "Tumor", "Metastasis"))
res.chisq <- chisq.test(table(metaData[["cluster"]], metaData[["pcat"]]))
R.oe <- (res.chisq$observed) / (res.chisq$expected)
R.oe <- R.oe[cluster.order, ]

pdf(sprintf("%s/Figure1F.pdf", figdir), width = 5, height = 8)
pheatmap(R.oe,
  cluster_rows = FALSE,
  color = c("#FEE8C8", "#FDBB84", "#FC8D59", "#EF6548"),
  breaks = c(0, 1, 1.5, 3, max(R.oe)),
  cluster_cols = FALSE,
  angle_col = 45,
  fontsize = 18, border_color = "white",
  display_numbers = matrix(ifelse(R.oe > 3, "+++", ifelse(R.oe > 1.5, "++", ifelse(R.oe > 1, "+", "+/-"))), nrow(R.oe)),
  number_color = "black"
)
dev.off()


# Heatmap of signature expression
load("./data/fib_avgExp_sig_score.RData")
annotation_col <- fibSignature.df %>%
  dplyr::select(geneset, cluster_max, tissue, category) %>%
  arrange(category, cluster_max) %>%
  as.data.frame() %>%
  "rownames<-"(.[, "geneset"]) %>%
  mutate(Tissue = tissue, Category = category) %>%
  select(Tissue, Category)

TissueColor <- brewer.pal("Dark2", n = length(unique(annotation_col$Tissue)))
names(TissueColor) <- unique(annotation_col$Tissue)

CategoryColor <- tableau_color_pal(palette = "Tableau 10")(length(unique(annotation_col$Category)))
names(CategoryColor) <- unique(annotation_col$Category)

ann_colors <- list(
  Category = CategoryColor,
  Tissue = TissueColor
)

plot.data <- t(fib.avgExp)[, fibSignature.df$geneset]

expr <- as.matrix(plot.data[cluster.order, rownames(annotation_col)])
mat2 <- scale(expr)
thresh.col <- 2
mat2 <- MinMax(mat2, max = thresh.col, min = (-1) * thresh.col)
colnames(mat2) <- colnames(expr)

pdf(sprintf("%s/FigureS1G.pdf", figdir), width = 10, height = 8)
pheatmap(mat2,
  cluster_rows = FALSE,
  color = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)),
  breaks = seq(-2, 2, length = 101),
  annotation_col = annotation_col,
  show_colnames = TRUE, show_rownames = TRUE,
  annotation_colors = ann_colors,
  border_color = "white", border = "white",
  cluster_cols = FALSE, angle_col = 90, fontsize = 12
)
dev.off()

# ROGUE index
rogue_res <- readRDS("./data/rogue_index.rds")
myColor <- c(
  "#E41B1B", "#4376AC", "#48A75A", "#87638F", "#D87F32", "#737690", "#D690C6", "#B17A7D", "#847A74", "#4285BF",
  "#204B75", "#588257", "#B6DB7B", "#E3BC06", "#FA9B93", "#E9358B", "#A0094E", "#999999", "#6FCDDC", "#BD5E95"
)
names(myColor) <- paste0("c", stringr::str_pad(as.character(seq_len(20)), 2, side = "left", "0"))

plotData <- rogue_res %>%
  tidyr::gather(key = clusters, value = ROGUE) %>%
  filter(!is.na(ROGUE))

g <- ggplot(data = plotData, aes(clusters, ROGUE, color = clusters)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast(shape = 16, position = position_jitter(0.2)) +
  scale_color_manual(values = myColor) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")
  ) +
  labs(x = "", y = "ROGUE index") +
  ylim(0, 1)

ggsave(sprintf("%s/FigureS1H.pdf", figdir), g, width = 7, height = 5)


# Heatmap of correspondence between clusters
load("./data/metadata_cluster.RData")
cluster.order <- metaData_clusters %>%
  arrange(cluster) %>%
  pull(cluster) %>%
  unique()
metadata_summary <- metaData_clusters %>%
  dplyr::select(bbknn, scVI) %>%
  group_by(scVI, bbknn) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100)

metadata_summary$bbknn <- factor(metadata_summary$bbknn, levels = cluster.order)
metadata_summary$scVI <- factor(metadata_summary$scVI, levels = cluster.order)

p <- metadata_summary %>%
  ggplot(aes(x = bbknn, y = scVI, width = 1, height = 1)) +
  geom_tile(aes(fill = clust_prop)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_viridis(na.value = "#E1E1E1", option = "D", begin = 0, end = 1)

ggsave(sprintf("%s/FigureS1J.pdf", figdir), p, width = 5, height = 4)
