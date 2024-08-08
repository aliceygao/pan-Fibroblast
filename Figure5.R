library(Seurat)
library(plyr)
library(tibble)
library(data.table)
library(ggplot2)
library(ggrastr)
library(ggpubr)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(grDevices)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(compositions)

# ------------- Figure 5 --------------
# Heatmap of fibroblast cluster correlations
load("./data/cellFreq_all.RData")
load("./data/cellModule.RData")

cellfreq.m <- reshape2::dcast(cellFreq.df, celltype ~ sample, value.var = "freq") %>%
    column_to_rownames(var = "celltype") %>%
    as.matrix()

corTable <- data.frame(cor(t(cellfreq.m)[pcatSample.ls[["Tumor"]], ]), check.names = F)
colorSets <- tableau_color_pal("Tableau 20")(20)
TMEC.color <- rep(colorSets[1:length(table(cellModule.df$Module))], as.numeric(table(cellModule.df$Module)))
names(TMEC.color) <- cellModule.df$Celltype
labelEach <- rownames(corTable)
labelColor <- TMEC.color[rownames(corTable)]

col.heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))
cols <- col.heat(100)
breaks <- c(seq(min(corTable), 0, length = 50), seq(0.01, max(corTable), length = 50))
col1 <- colorRamp2(breaks, cols)
if (is.null(names(labelColor)) || names(labelColor) != labelEach) {
    names(labelColor) <- labelEach
}
ha <- HeatmapAnnotation(
    Cluster = labelEach,
    col = list(Cluster = labelColor),
    show_legend = F
)
ht_opt$TITLE_PADDING <- unit(c(2.5, 2.5), "points")
pdf(sprintf("%s/Figure5A.pdf", figdir), width = 8, height = 8)
ht <- Heatmap(as.matrix(corTable),
    name = "Cor",
    top_annotation = ha,
    cluster_rows = T, cluster_columns = T,
    col = col1, border = F,
    row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
    show_row_dend = T, show_column_dend = T # ,
)
ht <- draw(ht)
dev.off()


# Dot plot of marker gene
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
load("./data/CM_marker_genes.RData")
load("./data/avgExp_CM.RData")
scale.func <- switch(
    EXPR = "radius",
    "size" = scale_size,
    "radius" = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
)

plot <- ggplot(data = data.plot, mapping = aes_string(x = "gene", y = "CM")) +
    geom_point_rast(mapping = aes_string(size = "pct.exp", fill = "avg.exp.scaled"), stroke = 0.5, shape = 21, color = "black") +
    scale_x_discrete(position = "top") +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale.func(range = c(0, 6)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = "Percent Expressed")) +
    theme_bw() +
    theme(
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, face = "italic"),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid = element_line(colour = "grey", linetype = "dashed"),
        panel.grid.major = element_line(colour = "grey", linetype = "dashed", size = 0.2)
    )
plot <- plot + scale_fill_gradientn(colours = SpatialColors)
ggsave(sprintf("%s/Figure5C.pdf", figdir, height = 2.5, width = 8))


# Density plot of KL divergence

load("./data/ST/KL_divergence.RData")

kl_pair.ls <- list()
for (j in 1:nrow(cellSinglePair)) {
    kl_pair.ls[[j]] <- klpair_all.df %>% filter(cellType1 == cellSinglePair$cellType1[j], cellType2 == cellSinglePair$cellType2[j])
}
plotData <- do.call(rbind, kl_pair.ls)
plotData$colors <- cellSinglePair$color

xlim1 <- ifelse(min(plotData$x) < min(klboot), min(plotData$x) - 0.2, min(klboot) - 0.2)
xlim2 <- ifelse(max(plotData$x) > max(klboot), max(plotData$x) + 0.2, max(klboot) + 0.2)
ylim1 <- density(klboot)
ylim1 <- ylim1$y[which.max(ylim1$y)] + 0.1

g <- ggplot(plotData) +
    geom_density(aes(x = permutT), alpha = 0.6, color = "lightgray", fill = "lightgray", data = data.frame(permutT = klboot)) +
    geom_point_rast(aes(x = x, y = y), size = 0.1, colour = plotData$colors) +
    geom_text_repel(aes(x = x, y = y, label = label),
        color = plotData$colors, segment.color = plotData$colors,
        size = 3,
        force_pull = 0,
        nudge_y = 0.3,
        direction = "x",
        angle = 90,
        hjust = 0,
        segment.size = 0.5,
        max.iter = 1e4, max.time = 1
    ) +
    xlim(xlim1, xlim2) +
    ylim(0, ylim1) +
    xlab("KL-divergence") +
    theme_classic() +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10), legend.position = "None"
    )
ggsave(sprintf("%s/Figure5D_KL_divergence.pdf", firdir), g, width = 4, height = 3, dpi = 300)


# Scatter plot of frequency
load("./data/colorSet.RData")
Cell_freq.tb <- readRDS("./data/fib_freq.rds")
data_wide <- reshape2::dcast(Cell_freq.tb %>% filter(pcat == "Tumor"), cluster ~ sample, value.var = "freq")
data_wide.m <- data_wide %>%
    column_to_rownames(var = "cluster") %>%
    as.matrix()
progenitor.v <- apply(data_wide.m[c("c05", "c03"), ], 2, sum)
data_wide_v2.m <- rbind(data_wide.m, progF = progenitor.v)
data <- as.data.frame(t(data_wide_v2.m))
data$cancerType <- Cell_freq.tb$tissue[match(rownames(data), Cell_freq.tb$sample)]
g <- ggscatter(data, x = "progF", y = "c04", size = 1.5, shape = 21, color = "cancerType") +
    coord_trans(clip = "off") +
    geom_smooth(method = "lm") +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson") +
    scale_color_manual(values = cancerType.color)
ggsave(sprintf("%s/FigureS5A_pct.pdf", figdir), g, width = 3, height = 2.6)

baseILR <- ilrBase(
    x = data_wide_v2.m,
    method = "basic"
)
cell_ilr <- as.matrix(ilr(data_wide_v2.m, baseILR))
colnames(cell_ilr) <- paste0("ILR_", 1:ncol(cell_ilr))

data <- as.data.frame(t(cell_ilr))
g <- ggscatter(data, x = "progF", y = "c04", size = 2, shape = 21, color = "black") +
    coord_trans(clip = "off") +
    geom_smooth(method = "lm") +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "pearson")
ggsave(sprintf("%s/FigureS5A_ILR.pdf", figdir), g, width = 3, height = 2.6)


# Boxplot of signature score
load("./data/TCGA_signature_score.RData")
gg.ls <- list()
for (i in 1:length(signature.order)) {
    gg.ls[[i]] <- ggboxplot(metaData, x = "stromalType", y = signature.order[i], outlier.shape = NA, color = "stromalType") +
        geom_point_rast(aes(color = stromalType), alpha = 0.3, size = 0.2, shape = 16, raster.dpi = 300, position = position_jitter(width = 0.2, height = 0)) +
        theme_bw() +
        stat_compare_means(paired = F, size = 2.8, method = "wilcox.test") +
        scale_colour_manual(values = scPalette(length(unique(metaData$stromalType)))) +
        guides(color = guide_legend(override.aes = list(size = 0.5)))
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 1, ncol = 4, common.legend = TRUE)
ggexport(multi.page, filename = sprintf("%s/FigureS5B.pdf", figdir), width = 8, height = 3)



# Cell module abundance
load("./data/moduleAbundance.RData")
sampleModule.df <- data.frame(
    sample = colnames(moduleAbundance.m),
    group = rownames(moduleAbundance.m)[apply(moduleAbundance.m, 2, which.max)],
    pct = apply(moduleAbundance.m, 2, max)
)
sampleModule.df <- left_join(sampleModule.df, metaData %>% select(sample, patient, tissue, Dataset, pcat) %>% distinct(), by = "sample")
module.color <- tableau_color_pal("Tableau 20")(20)[1:nrow(moduleAbundance.m)]
names(module.color) <- rownames(moduleAbundance.m)
sample.order <- sampleModule.df %>%
    arrange(group, -pct) %>%
    pull(sample)
ann_colors <- list(
    module = module.color
)
annCol <- data.frame(module = sampleModule.df$group[match(sample.order, sampleModule.df$sample)])
rownames(annCol) <- sample.order
plotData <- MinMax(moduleAbundance.m[, sample.order], max = 1, min = 0)
cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))(100)
breaks <- seq(0, 1, length = 101)

pdf(sprintf("%s/FigureS5C.pdf", figdir), height = 3, width = 10)
pheatmap(plotData,
    cluster_rows = FALSE, color = cols, breaks = breaks,
    annotation_col = annCol, annotation_colors = ann_colors,
    border_color = "white", border = "white", show_colnames = FALSE,
    cluster_cols = FALSE, angle_col = "45", fontsize = 10
)
dev.off()


# Correlation between fibroblast and immune cells
load("./data/celltype_list.RData")
cor.df <- read_csv("./data/cell_freq_correlation.csv")
cor.df$cor <- as.numeric(as.matrix(cor.df)[, "c04"])
cor_sub.df <- cor.df %>%
    filter(celltype %in% celltype.ls[["immunecells"]]) %>%
    arrange(-cor) %>%
    filter(cor >= 0)
celltype.category.color <- scPalette(nrow(cor_sub.df))
names(celltype.category.color) <- cor_sub.df$celltype
g <- ggdotchart(cor_sub.df,
    x = "celltype", y = "cor",
    color = "celltype", palette = celltype.category.color, sorting = "descending",
    add = "segments", add.params = list(color = "lightgray", size = 0.5), dot.size = 6, title = "",
    label = round(cor_sub.df$cor, 2), font.label = list(color = "white", size = 7, vjust = 0.5), ggtheme = theme_pubr()
) + xlab("") + ylab("") +
    geom_hline(yintercept = 0.15, linetype = 2, color = "lightgray") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
ggsave(sprintf("%s/FigureS5D.pdf", figdir), g, width = 4, height = 5)

metaData <- readRDS("./data/TCGA_sig_score.rds")
data.ls <- split(metaData[, c(celltype.ls[["immunecells"]], celltype.ls[["fibroblast"]])], metaData$cancerType)
M.ls <- lapply(data.ls, function(x) cor(x, method = "pearson"))
M.df <- data.frame(cor = unlist(lapply(M.ls, function(x) x[, "c04"])), celltype = unlist(lapply(M.ls, rownames)), cancerType = rep(names(M.ls), unlist(lapply(M.ls, nrow))))
celltype_order <- M.df %>%
    group_by(celltype) %>%
    summarize(mean = mean(cor)) %>%
    arrange(-mean) %>%
    pull(celltype)
M.df$celltype <- factor(M.df$celltype, levels = celltype_order)
y.mean <- M.df %>%
    group_by(celltype) %>%
    summarize(mean = mean(cor)) %>%
    pull(mean) %>%
    mean()
cancerType.color <- scPalette(length(unique(M.df$cancerType)))
names(cancerType.color) <- unique(M.df$cancerType)

g <- ggplot(cor.df %>% filter(celltype %in% celltype.ls[["immunecells"]]), aes(x = celltype, y = cor)) +
    geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.1, size = 0.2) +
    geom_jitter(aes(color = cancerType), alpha = 1, size = 1, shape = 20, position = position_jitter(w = 0.1, h = 0.1)) +
    theme_classic() +
    ylab("Correlation with c04") +
    geom_hline(yintercept = y.mean, linetype = "dashed", color = "blue") +
    scale_color_manual(values = cancerType.color) +
    guides(colour = guide_legend(
        ncol = 2, override.aes = list(size = 3),
        label.theme = element_text(size = 6),
        keyheight = 0.5
    )) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
    )
ggsave(sprintf("%s/FigureS5E.pdf", figdir), g, width = 8, height = 4.1)

cellPair.df <- data.frame(cluster1 = c("c04", "c19", "c19", "c05"), cluster2 = c("Macro_SPP1", "CD4.c08.Treg.FOXP3", "cDC3_LAMP3", "CD8.c06.Temra.CX3CR1"), stringsAsFactors = FALSE)
for (i in 1:nrow(cellPair.df)) {
    cluster1 <- cellPair.df$cluster1[i]
    cluster2 <- cellPair.df$cluster2[i]
    gg.ls[[i]] <- ggscatter(metaData, x = cluster1, y = cluster2, size = 0.1, color = "cancerType") +
        labs(x = paste0(cluster1, " Score"), y = paste0(cluster2, " Score"), title = "") +
        geom_smooth(method = "lm", color = "black") +
        guides(
            color = guide_legend(order = 1, override.aes = list(size = 5), ncol = 4),
            legend.direction = "horizontal"
        ) +
        coord_trans(clip = "off") +
        theme(legend.position = "bottom", legend.title = element_blank()) +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), ) +
        stat_regline_equation(aes(label = paste(..eq.label..)), vjust = 2.5) +
        scale_color_manual(values = cancerType.color)
}

multi.page <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = 2, common.legend = TRUE)
ggexport(multi.page, filename = sprintf("%s/FigureS5F.pdf", figdir), width = 10, height = 8)

results_list <- readRDS("./data/scaden_result.rds")
for (j in 1:nrow(cellPair.df)) {
    ct1 <- cellPair.df$cluster1[j]
    ct2 <- cellPair.df$cluster2[j]
    cm1 <- list()
    for (i in seq_along(results_list)) {
        if (ct1 %in% colnames(results_list[[i]]) & ct2 %in% colnames(results_list[[i]])) {
            print(unique(results_list[[i]]$ct_abbr))
            cm1[[i]] <- results_list[[i]] %>% select(as.name(ct1), as.name(ct2), ct_abbr)
        }
    }
    cm <- do.call(rbind, cm1)
    cm$ct_abbr <- factor(cm$ct_abbr)
    cm.ls <- split(cm, cm$ct_abbr)
    cor.v <- unlist(lapply(cm.ls, function(x) cor(x[, 1], x[, 2])))
    p.v <- unlist(lapply(cm.ls, function(x) cor.test(x[, 1], x[, 2])$p.value))
    ct_pick <- names(cm.ls)[p.v < 0.05 & cor.v > 0.15]
    cm_sub <- cm %>% filter(ct_abbr %in% ct_pick)
    sp <- ggscatter(cm_sub, x = ct1, y = ct2, color = "ct_abbr", shape = 21, add = "reg.line", conf.int = F) + scale_color_manual(values = cancerType.color)
    sp <- facet(sp, facet.by = "ct_abbr", scales = "free", nrow = 1)
    g <- sp + stat_cor(aes(color = "black"), method = "pearson")
    ggsave(sprintf("%s/FigureS5G.pdf", figdir), g, width = 3 * length(unique(cm_sub$ct_abbr)), height = 3.5)
}
