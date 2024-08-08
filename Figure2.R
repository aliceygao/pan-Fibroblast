
library(Seurat)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(readr)
library(grDevices)
library(ggraph)
library(igraph)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(rstatix)
library(ggpubr)

# ------------- Figure 2 --------------
# Boxplot of fibroblast cluster proportions

seu <- readRDS("./data/all_fib_cells.rds")
metaData <- seu@meta.data

plot.data <- metaData %>%
  group_by(sample, cluster) %>%
  summarize(count = n()) %>%
  mutate(clust_total = sum(count)) %>%
  mutate(clust_prop = count / clust_total * 100)
plot.data <- left_join(plot.data, metaData[, "sample", "tissue", "pcat"], by = "sample")

pcat.color <- c(
  "Healthy" = "#A6CEE3", "Non-inflamed" = "#E8B523", "Inflamed" = "#ac8412",
  "Normal" = "#EA3C8C", "Fetal" = "#b7d9b7", "Tumor" = "#56A255", "Metastasis" = "#E41A1C"
)

g <- ggboxplot(plot.data %>% filter(cluster == clusterVar, tissue == tissueVar),
  x = "pcat", y = "clust_prop",
  add = c("jitter"), alpha = 0.2, add.param = list(alpha = 0.5),
  color = "pcat", fill = "pcat",
  palette = pcat.color, outlier.shape = NA
) +
  ylab("Proportion (%)") + xlab("") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
ggsave(sprintf("%s/Figure2A_%s_%s.pdf", figdir, clusterVar, tissueVar), g, width = 1.5, height = 3)


# Barplot of fibroblast compositions across phenotypes
fibcluster.map <- readRDS("./data/fib_cluster_category.rds")
metaData$clusterCategory <- plyr::mapvalues(metaData$cluster, from = fibcluster.map$cluster, to = fibcluster.map$category)
clusterCategory.order <- rev(c("Fibroblast-like cells", "CTNNB1+ fibroblast", "Myofibroblasts", "STMN1+ fibroblast", "Progenitor-like fibroblasts", "MMP1+ fibroblast", "CD74+ fibroblast", "Inflammatory fibroblasts", "Tissue-resident fibroblasts"))
pcat.order <- c("Healthy", "Inflamed", "Tumor", "Metastasis")

color.category <- c("#4376AC", "#D690C6", "#BD5E95", "#6ECDDC", "#49A75A", "#A0094E", "#87638F", "#847A74", "#B17A7D")
names(color.category) <- clusterCategory.order

plot.data <- metaData %>%
  filter(pcat %in% pcat.order) %>%
  dplyr::select(pcat, clusterCategory) %>%
  dplyr::mutate(pcat = factor(pcat, levels = pcat.order), clusterCategory = factor(clusterCategory, levels = clusterCategory.order)) %>%
  group_by(pcat, clusterCategory) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100)

g <- ggplot(plot.data, aes(x = pcat, y = clust_prop, fill = clusterCategory)) +
  geom_bar(stat = "identity", color = "white", width = 1, size = 0.5) +
  scale_fill_manual(values = color.category) +
  facet_grid(. ~ pcat, scales = "free", space = "free") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank()
  ) +
  labs(
    x = "", y = "Proportion (%)"
  )
ggsave(sprintf("%s/Figure2B.pdf", figdir), g, width = 5, height = 5)


# Heatmap of fibroblast activation score
load("./data/avg_FAScore_pcat.RData")
mat <- t(scale(t(avgExp)))
thresh.col <- 2
mat2 <- MinMax(mat, max = thresh.col, min = (-1) * thresh.col)
colnames(mat2) <- colnames(avgExp)

pcat.order <- c("Healthy", "Normal", "Inflamed", "Tumor", "Metastasis")
plotData <- mat2[, pcat.order]

MarkGenes <- c("COL4A1", "COL1A1", "COL18A1", "ISG15", "HLA-A", "IFI6", "PDGFRB", "THBS2", "TIMP1", "MMP14", "SULF1", "TGFBI", "IGFBP4", "COL6A1")
markerGeneL <- rownames(plotData) %in% MarkGenes
markerGenes <- rownames(plotData)[markerGeneL]

rowHa <- rowAnnotation(
  link = anno_mark(
    at = which(markerGeneL),
    labels = markerGenes, labels_gp = gpar(fontsize = 8, fontface = "italic"), link_gp = gpar(lwd = 0.3), padding = 0.5
  ),
  width = unit(0.2, "cm") + max_text_width(markerGenes, gp = gpar(fontsize = 8))
)

ht_opt$TITLE_PADDING <- unit(c(2.5, 2.5), "points")
cols <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
breaks <- seq(-2, 2, length = 100)
col1 <- colorRamp2(breaks, cols)

pdf(sprintf("%s/Figure2C.pdf", figdir), width = 3, height = 4)
Heatmap(as.matrix(plotData),
  name = "z-score",
  right_annotation = rowHa,
  cluster_rows = T, cluster_columns = F,
  col = col1, border = F,
  show_row_names = FALSE, show_column_names = TRUE,
  use_raster = TRUE, raster_resize_mat = max
)
dev.off()


# Boxplot of fibroblast activation score
load("./data/all_FAScore.RData")
metadata_sub <- metaData %>% filter(tissue != "HNSC", pcat %in% c("Healthy", "Inflamed", "Normal", "Tumor", "Metastasis"))
metadata_sub$tissue_pcat <- paste0(metadata_sub$tissue, "_", metadata_sub$pcat)
metadata_sub$FAscore <- FAScore.df$score[match(metadata_sub$cellName, FAScore.df$cellName)]
data_median <- metadata_sub %>%
  group_by(tissue_pcat) %>%
  dplyr::summarize(mean = mean(FAscore))
data_median$tissue <- sapply(strsplit(as.character(data_median$tissue_pcat), split = "_"), "[", 1)
tissue_pcat.order <- data_median %>%
  arrange(tissue, mean) %>%
  pull(tissue_pcat)
metadata_sub$tissue_pcat <- factor(metadata_sub$tissue_pcat, levels = rev(tissue_pcat.order))

g <- ggplot(metadata_sub, aes(x = tissue_pcat, y = FAscore, fill = pcat)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, size = 0.2) +
  scale_fill_manual(values = pcat.color) +
  facet_grid(rows = vars(tissue), scales = "free", space = "free") +
  theme_bw() +
  xlab("") +
  ylab("Fibrobalst activation score") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()
ggsave(sprintf("%s/Figure2D.pdf", figdir), g, width = 4, height = 10)


# Volcano plot of DE genes
load("./data/Inflamed_Tumor_DEgenes.RData")
DEgenes$Topranked <- ifelse(DEgenes$gene %in% topGenes, "Top-ranked", "Others")
markerSig <- DEgenes %>%
  filter(abs(avg_log2FC) > 0.25, p_val_adj < 0.01) %>%
  pull(gene)
DEgenes$Significance <- ifelse(DEgenes$gene %in% markerSig, "TRUE", "FALSE")
DEgenes$Significance[DEgenes$avg_log2FC > 0.25 & DEgenes$Significance == "TRUE"] <- "UP"
DEgenes$Significance[DEgenes$avg_log2FC < -0.25 & DEgenes$Significance == "TRUE"] <- "DW"
DEgenes$pAdj.log <- -log10(DEgenes$p_val_adj)

g <- ggplot(DEgenes, aes(x = avg_log2FC, y = pAdj.log)) +
  geom_point_rast(aes(color = Significance)) +
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black")
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey") +
  geom_text_repel(
    data = subset(DEgenes, Topranked == "Top-ranked"),
    aes(label = gene),
    size = 4.6, fontface = "italic", color = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
  ) +
  ggtitle(titleName) +
  ylab("-log10(Padj)")
ggsave(sprintf("%s/Figure2E.pdf", figdir), g1, width = 6, height = 6)


# Barplot of enriched pathway
load("./data/GSEA_term_fibclusters.RData")

enrichgo.top <- gsea_all.df %>% filter(Term %in% topPathways)
term.order <- reshape2::dcast(enrichgo.top, Term ~ Cluster, value.var = "Logqvalue") %>%
  arrange(-c19_Inflamed, -c20_Inflamed, c19_Tumor, c20_Tumor) %>%
  pull(Term)

cluster.order <- c("c19_Inflamed", "c20_Inflamed", "c19_Tumor", "c20_Tumor")
colorSet <- c("#00a6ac", "#94d6da", "#fdb933", "#dec674")
names(colorSet) <- cluster.order
enrichgo.top$Term <- factor(enrichgo.top$Term, levels = rev(term.order))
enrichgo.top$Cluster <- factor(enrichgo.top$Cluster, levels = cluster.order)

g <- ggplot(enrichgo.top, aes(x = Term, y = Logqvalue, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.6) + 
  scale_fill_manual(values = colorSet) +
  facet_grid(~Cluster, scales = "free", space = "free_y") +
  theme_bw() +
  xlab("") +
  ylab("-log10(qvalue)") +
  theme(
    strip.text.x = element_text(colour = "black", angle = 360, size = 12, hjust = 0, face = "bold"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 12),
    axis.ticks = element_line(size = 1, colour = "black"),
    axis.line = element_line(size = 1, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )
g <- g + coord_flip()
ggsave(sprintf("%s/Figure2F.pdf", figdir), g, width = 12, height = 5)


load("./data/GSEA_term_InfvsT.RData")
GSEA_term_Sig$Logqvalue[which(GSEA_term_Sig$Cluster == "Inflamed")] <- -GSEA_term_Sig$Logqvalue[which(GSEA_term_Sig$Cluster == "Inflamed")]
GSEA_term_Sig$Logqvalue <- MinMax(GSEA_term_Sig$Logqvalue, min = -10, max = 10)

g <- ggbarplot(GSEA_term_Sig,
  x = "Description", y = "Logqvalue",
  fill = "Cluster",
  color = "white",
  palette = pcat.color,
  sort.val = "asc",
  sort.by.groups = FALSE,
  x.text.angle = 90,
  ylab = "-Log10(qvalue)",
  xlab = FALSE, width = 0.98,
  legend.title = "Tumor vs Inflamed"
) + coord_flip() + ylab("")
ggsave(sprintf("%s/Figure2G.pdf", figdir), g, width = 7, height = 4.1)


# Enrichment of fibroblast clusters
load("./fib_metaData.RData")
Roe.ls <- list()
tissues <- unique(metaData$tissue)
for (i in 1:length(tissues)) {
  metadata_sub <- metaData %>% filter(tissue == tissues[i])
  res.chisq <- chisq.test(table(metadata_sub[["cluster"]], metadata_sub[["pcat"]]))
  R.oe <- (res.chisq$observed) / (res.chisq$expected)
  R.oe <- R.oe[cluster.order, ]
  count.m <- table(metadata_sub$cluster, metadata_sub$pcat)
  cluster.pick <- rownames(count.m)[apply(count.m, 1, sum) > 20]

  Roe.df <- R.oe[cluster.pick, ] %>%
    as.data.frame.matrix() %>%
    rownames_to_column(var = "cluster")
  Roe_long <- reshape2::melt(Roe.df, id.vars = c("cluster"), measure.vars = colnames(Roe), variable.name = "pcat", value.name = "Roe")
  Roe_long$tissue <- tissues[i]
  Roe.ls[[tissues[i]]] <- Roe_long
}


plotData <- do.call(rbind, Roe.ls)
plotData$color <- NA
plotData$color[plotData$Roe > 1] <- "Enrichment"
plotData$color[plotData$Roe < 1] <- "Depletion"

pcat.order <- c("Fetal", "Healthy", "Normal", "Non-inflamed", "Inflamed", "Tumor", "Metastasis")
tissue.order <- c("Breast", "Colon", "Lung", "Pancreas", "Ovary", "Stomach", "Skin", "Kidney", "Liver", "HNSC", "Synovium")
plotData$pcat <- factor(plotData$pcat, levels = pcat.order)
plotData$tissue <- factor(plotData$tissue, levels = tissue.order)
color_Roe <- c(Enrichment = "#E41A1C", Depletion = "#377EB8")

g_panel <- ggplot(plotData, aes(x = tissue, y = pcat, size = Roe, color = color)) +
  geom_point(alpha = 0.9) +
  facet_wrap(~cluster, ncol = 5, scales = "free_y", strip.position = "top", dir = "h") +
  theme_classic() +
  scale_color_manual(values = color_Roe) +
  scale_size(range = c(0, 5)) +
  rotate_x_text(60) +
  rremove("xlab") +
  labs(y = "RO/E") +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line = element_line(size = 0.5), axis.ticks = element_line(size = 0.5),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2),
    strip.text.x = element_text(face = "bold", color = "white", size = 12)
  )

g <- ggplot_gtable(ggplot_build(g_panel))
stripr <- which(grepl("strip-", g$layout$name))
fills <- cluster.color[sort(unique(plotData$cluster))]
k <- 1
for (i in stripr) {
  j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- fills[k]
  g$grobs[[i]]$heights <- unit(0.25, "in")
  k <- k + 1
}
ggsave(plot = plot(g), filename = sprintf("%s/FigureS2A.pdf", figdir), width = 15, height = 6.5)


# Heatmap of regulon activity score
load("./data/fib_regulons.RData")
auc_mtx <- readRDS("./data/fib_auc_mtx.rds")
rss.m <- read_csv(sprintf("%s/cluster_RSS.csv", figdir)) %>%
  column_to_rownames(var = "X1") %>%
  as.matrix() %>%
  t()

auc_mtx$cluster <- metaData$cluster[match(auc_mtx$cellName, metaData$cellName)]
auc_mtx_avg.ls <- lapply(
  X = unique(x = auc_mtx$cluster),
  FUN = function(cluster) {
    data.use <- auc_mtx[auc_mtx$cluster == cluster, 2:(ncol(x = auc_mtx) - 1), drop = FALSE]
    avg.exp <- apply(
      X = data.use,
      MARGIN = 2,
      FUN = mean
    )
    return(avg.exp)
  }
)
auc_mtx_avg.m <- do.call(rbind, auc_mtx_avg.ls)
rownames(auc_mtx_avg.m) <- unique(x = auc_mtx$cluster)

top10regulons.ls <- list()
for (i in 1:ncol(rss.m)) {
  regulons <- names(sort(rss.m[, i], decreasing = TRUE))[1:10]
  top10regulons.ls[[i]] <- regulons
}
names(top10regulons.ls) <- colnames(rss.m)
topGenes <- unique(unlist(top10regulons.ls))
expr <- t(auc_mtx_avg.m)[topGenes, cluster.order]
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

pdf(sprintf("%s/FigureS2B.pdf", figdir), width = 7, height = 8)
Heatmap(plotData,
  name = "z-score",
  top_annotation = ha,
  right_annotation = rowHa,
  cluster_rows = T, cluster_columns = T,
  col = col1, border = F,
  show_row_names = FALSE, show_column_names = TRUE
)
dev.off()


# TF regulon target gene network
tf.df <- readRDS("./data/regulon_target.rds")
gene.DE.tb <- readRDS("./data/gene_DE_list.rds")

network.plot.tb <- as.data.table(tf.df)[source %in% TF.list, ]
network.plot.tb[, TF := gsub("\\(\\+\\)", "", source)]
TF.rename <- gsub("\\(\\+\\)", "", TF.list)
g.info.tb <- gene.DE.tb[cluster == clusterVar & geneID %in% unique(c(TF.rename, network.plot.tb$target)), c("geneID", "comb.ES")]
node.tb <- network.plot.tb
node.tb <- left_join(node.tb, g.info.tb, by = c("target" = "geneID"))
node.tb$type <- "target gene"
node.tb[target %in% TF.rename, type := "TF"]
node.tb[, name := target]
node.tb.top <- node.tb %>%
  filter(!is.na(comb.ES), comb.ES > 0) %>%
  top_n(40, wt = comb.ES) %>%
  arrange(-comb.ES)
node.tb.top[, source := TF]

edge.list <- node.tb.top %>%
  filter(type != "TF") %>%
  select(source, target)
node.list <- node.tb.top %>% select(name, type, comb.ES)
ES.breaks <- seq(0, 0.6, by = 0.1)

if (!TF.rename %in% node.list$name) {
  node.list <- rbind(node.list, data.table(name = TF.rename, type = "TF", comb.ES = max(ES.breaks)))
} else {
  node.list[type == "TF", comb.ES := max(ES.breaks)]
}
node.list[, ES := comb.ES]
node.list[comb.ES > ES.breaks[length(ES.breaks)], ES := ES.breaks[length(ES.breaks)]]
node.list[comb.ES < ES.breaks[1], ES := ES.breaks[1]]
node.list$nodeSize <- 8 + (node.list$ES - min(node.list$ES)) / (max(node.list$ES) - min(node.list$ES)) * 10
flareGraph <- graph_from_data_frame(edge.list, vertices = node.list)
palette <- brewer.pal(9, "PuRd")

g <- ggraph(flareGraph, "circlepack", circular = TRUE) +
  geom_edge_link(aes(alpha = ..index..)) +
  scale_edge_alpha("Direction", guide = "none") +
  geom_node_point(aes(color = ES, size = nodeSize)) +
  geom_node_text(aes(label = name), color = "black", size = 3.2, nudge_x = 0.3, nudge_y = 0) +
  scale_color_gradient(low = palette[1], high = palette[6], breaks = ES.breaks) +
  ggforce::theme_no_axes() + ggforce::theme_no_axes()
ggsave(sprintf("%s/FigureS2D.pdf", figdir), g, width = 7, height = 6)


# Dot plot of fibroblast proportions
metaData <- readRDS("./data/fib_metaData.RData")
load("./data/tissueDataset_color.RData")

Dataset.pcat <- table(metaData$Dataset, metaData$pcat)
Dataset.pick <- rownames(Dataset.pcat)[apply(Dataset.pcat[, c("Normal", "Tumor")], 1, function(x) sum(x > 20)) == 2] 
metadata <- metaData %>%
  filter(Dataset %in% Dataset.pick, pcat %in% c("Normal", "Tumor")) %>%
  as.data.frame()
metadata$Dataset_pcat <- paste0(metadata$Dataset, "_", metadata$pcat)

plot.data <- as.data.frame(metadata) %>%
  dplyr::select(cluster, Dataset_pcat) %>%
  group_by(Dataset_pcat, cluster) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::mutate(clust_total = sum(count)) %>%
  dplyr::mutate(clust_prop = count / clust_total * 100)

plot.data$Dataset <- sapply(strsplit(plot.data$Dataset_pcat, split = "_"), "[", 1)
plot.data$pcat <- sapply(strsplit(plot.data$Dataset_pcat, split = "_"), "[", 2)
plot.data$tissue <- metadata$tissue[match(plot.data$Dataset, metadata$Dataset)]

cluster.order <- c("c04", "c11", "c16", "c19", "c20")

tissueDataset.color <- tissueDataset_color.df$color
tissueDataset.shape <- tissueDataset_color.df$shape
names(tissueDataset.color) <- tissueDataset_color.df$tissue_Dataset
names(tissueDataset.shape) <- tissueDataset_color.df$tissue_Dataset

gg.ls <- list()
for (i in 1:length(cluster.order)) {
  plot.data.sub <- plot.data %>% filter(cluster == cluster.order[i])
  plot.data.sub$tissue <- metadata$tissue[match(plot.data.sub$Dataset, metadata$Dataset)]
  plot.data.sub$tissue_Dataset <- paste0(plot.data.sub$tissue, "_", plot.data.sub$Dataset)

  gg.ls[[i]] <- ggplot(plot.data.sub, aes(x = pcat, y = clust_prop, group = Dataset)) +
    geom_line(size = 0.1) +
    geom_point(size = 3, aes(color = tissue_Dataset, shape = tissue_Dataset)) +
    scale_color_manual(values = tissueDataset.color) +
    scale_shape_manual(values = tissueDataset.shape) +
    scale_x_discrete("") +
    theme_classic() +
    ggtitle(cluster.order[i]) +
    theme(
      legend.position = "top",
      axis.line.y = element_line(size = .5),
      axis.line.x = element_line(size = .5),
      plot.title = element_text(hjust = 0.5)
    )
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = 3, common.legend = TRUE, legend = "right")
ggexport(multi.page, filename = sprintf("%s/FigureS2E.pdf", figdir), width = 7, height = 6, useDingbats = F)


# Violin plot for pathway signature expression
load("./data/sigVar_metabolic.RData")
load("./data/myoF_pcat_colors.RData")
sigScore <- readRDS("./data/all_Metabolic_sigscore.rds")
sigScore_long <- as.data.frame(sigScore) %>%
  rownames_to_column(var = "sig") %>%
  pivot_longer(!sig, names_to = "cellName", values_to = "score")
sigScore_meta <- left_join(sigScore_long, metaData, by = "cellName")
sigScore_meta_sub <- sigScore_meta %>% filter(sig %in% sigVar, cluster %in% cluster.order, pcat %in% c("Inflamed", "Tumor"))
sigScore_meta_sub$cluster_pcat <- paste0(sigScore_meta_sub$cluster, "_", ssScore_meta_sub$pcat)
sigScore_meta_sub$cluster_pcat <- factor(sigScore_meta_sub$cluster_pcat, levels = paste0(rep(cluster.order, each = 2), "_", rep(c("Inflamed", "Tumor"), length(cluster.order))))

g <- ggplot(sigScore_meta_sub, aes(x = cluster_pcat, y = score, fill = cluster_pcat, color = cluster_pcat)) +
  geom_violin(
    scale = "width",
    adjust = 1,
    trim = TRUE,
  ) +
  stat_summary(
    fun = "mean",
    geom = "crossbar",
    width = 0.5, linewidth = 0.15,
    colour = "black"
  ) +
  facet_wrap(~sig, nrow = 4, scales = "free_y", strip.position = "top", dir = "h") +
  theme_classic() +
  ylab("Gene score") +
  xlab("") +
  scale_color_manual(values = cluster_pcat.color) +
  scale_fill_manual(values = cluster_pcat.color) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "none", strip.text = element_text(size = 10)
  )
ggsave(sprintf("%s/FigureS2F.pdf", figdir), g, width = 8, height = 10)
