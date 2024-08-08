library(data.table)
library(dplyr)
library(gridExtra)
library(grid)
library(ggrepel)
library(readr)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(grDevices)
library(circlize)

# ------------- Figure 6 --------------
# Chord diagram of ligand-receptor interactions
small_legend <- function(fontsize = 5, keysize = .1, marginsize = c(-.1, 0, 0, 0), ...) {
    small_legend_theme <- theme(
        legend.title = element_text(size = fontsize),
        legend.text = element_text(size = fontsize),
        legend.key.size = unit(keysize, "lines"),
        legend.margin = margin(marginsize[1], marginsize[2], marginsize[3], marginsize[4], unit = "cm"),
        ...
    )
    return(small_legend_theme)
}

load("./data/ligandreceptor.RData")
colorF <- c("c05" = "#FDB462", "c04" = "#87638F", "c19" = "#84cbda")
colorTM <- c("#FB8072", "#8DD3C7", "#BEBADA", "#80B1D3", "#80B1D3", "#ed995a", "#FB8072", "#f3c89d", "#8cc598", "#3d8ea9")
names(colorTM) <- c("CD8.c06.Temra.CX3CR1", "Macro_LYVE1", "Macro_FN1", "Mono_CD14", "Macro_SPP1", "CD4.c08.Treg.FOXP3", "Macro_INHBA", "Macro_NLRP3", "pDC_LILRA4", "cDC3_LAMP3")

plotData <- ligandReceptor.df %>%
    filter(clusterA %in% fibClusterVar, clusterB %in% TMClusterVar) %>%
    mutate(source_ligand = paste0(ligand_complex, "_", source), target_receptor = paste0(receptor_complex, "_", target)) %>%
    dplyr::select(source_ligand, target_receptor, ligand = "ligand_complex", receptor = "receptor_complex", source, target) %>%
    dplyr::arrange(target, ligand, receptor)
color_use <- c(colorF, colorTM)[unique(plotData$target)]

sep <- ">@<"
vertices <- data.frame(name = c("root", unique(plotData$source), unique(plotData$target)), type = NA, celltype = NA, label = NA)
for (i in c(unique(plotData$source))) {
    vertices <- rbind(vertices, data.frame(name = paste0(i, sep, "ligand"), type = NA, celltype = i, label = NA))
}
for (i in c(unique(plotData$target))) {
    vertices <- rbind(vertices, data.frame(name = paste0(i, sep, "receptor"), type = NA, celltype = i, label = NA))
}

vertices2 <- data.frame()
for (i in seq_len(nrow(plotData))) {
    vertices2 <- rbind(vertices2, data.frame(
        name = paste0(plotData$source[i], sep, plotData$ligand[i]),
        type = "ligand", celltype = plotData$source[i], label = plotData$ligand[i]
    ))
}
for (i in seq_len(nrow(plotData))) {
    vertices2 <- rbind(vertices2, data.frame(
        name = paste0(plotData$target[i], sep, plotData$receptor[i]),
        type = "receptor", celltype = plotData$target[i], label = plotData$receptor[i]
    ))
}
vertices2 <- vertices2 %>% dplyr::arrange(celltype)
vertices <- rbind(vertices, vertices2)
vertices <- vertices[!duplicated(vertices, fromLast = T), ]

edges <- data.frame(
    from = c(rep("root", length(c(unique(plotData$source), unique(plotData$target)))), unique(plotData$source), unique(plotData$target)),
    to = c(
        unique(plotData$source), unique(plotData$target), paste0(unique(plotData$source), sep, "ligand"),
        paste0(unique(plotData$target), sep, "receptor")
    ), name = NA, edgeColor = NA
)

edges2 <- data.frame()
for (i in seq_len(nrow(plotData))) {
    edges2 <- rbind(edges2, data.frame(
        from = paste0(plotData$source[i], sep, "ligand"),
        to = paste0(plotData$source[i], sep, plotData$ligand[i]),
        name = paste0(plotData$source[i], sep, plotData$ligand[i], sep, plotData$target[i], sep, plotData$receptor[i]),
        edgeColor = plotData$target[i]
    ))
}

for (i in seq_len(nrow(plotData))) {
    edges2 <- rbind(edges2, data.frame(
        from = paste0(plotData$target[i], sep, "receptor"),
        to = paste0(plotData$target[i], sep, plotData$receptor[i]),
        name = paste0(plotData$source[i], sep, plotData$ligand[i], sep, plotData$target[i], sep, plotData$receptor[i]),
        edgeColor = plotData$target[i]
    ))
}
edges2 <- edges2 %>% dplyr::arrange(edgeColor)
edges <- rbind(edges, edges2)
edges <- edges[!duplicated(edges, fromLast = T), ]

from <- match(paste0(plotData$source, sep, plotData$ligand), vertices$name)
to <- match(paste0(plotData$target, sep, plotData$receptor), vertices$name)

vertices$id <- NA
myleaves <- which(is.na(match(vertices$name, edges$from)))
nleaves <- length(myleaves)
vertices$id[myleaves] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$hjust <- ifelse(vertices$angle < -90, 1, 0)
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)

mygraph <- igraph::graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
pl <- ggraph(mygraph, layout = "dendrogram", circular = TRUE) +
    geom_conn_bundle(data = get_con(from = from, to = to), tension = 0.5) +
    scale_edge_color_manual(values = color_use) +
    geom_node_point(aes(filter = leaf, x = x * 1.01, y = y * 1.01, size = 5, color = celltype, shape = type)) +
    theme_void() + coord_fixed() + scale_shape_manual(values = c(ligand = 19, receptor = 1)) +
    geom_node_text(aes(x = x * 1.05, y = y * 1.05, filter = leaf, label = label, angle = angle, hjust = hjust, color = celltype), size = 3, alpha = 1) +
    scale_colour_manual(values = color_use) +
    small_legend(keysize = 0.5) +
    theme(
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
ggsave(sprintf("%s/Figure6B.pdf", figdir), pl, width = 6, height = 6)


# Boxplot of ligand-receptor pair co-expression
load("./data/ST_ligand_receptor.RData")
ligandReceptor.df <- ligandReceptor.df %>%
    filter(ligand %in% rownames(ST.exp), receptor %in% rownames(ST.exp)) %>%
    distinct(ligand, receptor)
ligandReceptor.df$lrpair <- paste0(ligandReceptor.df$ligand, "_", ligandReceptor.df$receptor)
ligandReceptor.ls <- df2list(ligandReceptor.df[, c("ligand", "receptor")])

LigandR_mean.ls <- lapply(ligandReceptor.ls, function(ct) {
    apply(ST.exp[ct, ], 2, mean)
})

LigandR_mean.m <- do.call(rbind, LigandR_mean.ls)
rownames(LigandR_mean.m) <- unlist(lapply(ligandReceptor.ls, function(x) paste0(x[1], "_", x[2]))) # step 2
avgExp <- as.data.frame(t(LigandR_mean.m)) %>% rownames_to_column(var = "cellName") # step 3
data.plot <- left_join(avgExp, as.data.frame(spot_abundance_category[, fibCluster]) %>% rownames_to_column(var = "cellName"), by = "cellName")
colnames(data.plot)[ncol(data.plot)] <- fibCluster
data_long <- reshape2::melt(data.plot, id.vars = c("cellName", fibCluster), measure.vars = ligandReceptor.df$lrpair, variable.name = "LigandReceptor", value.name = "LRmean")
data_long[, fibCluster] <- factor(data_long[, fibCluster])

g <- ggplot(data_long, aes_string(x = "LigandReceptor", y = "LRmean", fill = fibCluster)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, lwd = 0.3, position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = ST.color) +
    theme_bw() +
    xlab("") +
    ylab("Normalized ligand-receptor\naverage co-expression") +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust = 1, size = 8), panel.grid = element_blank(), axis.title = element_text(size = 8), axis.text.y = element_text(size = 8))
ggsave(sprintf("%s/Figure6E.pdf", figdir), width = 8, height = 3)


# Dot plot of signature score
load("./data/bulk_exp_stim.RData")
Exp_bygCluster <- lapply(Signature.ls, function(genes) colMeans(normData[genes, df$SeqNo]))
Exp.m <- do.call(rbind, Exp_bygCluster)
CTColor <- c("NSCLC" = "#CE6094", "CRC" = "#6480BA")
sigmarker <- rownames(Exp.m)
data.df <- sampleMeta %>%
    filter(Stim %in% c("Ctrl", "SPP1")) %>%
    arrange(Donor, Stim)
data.df$Stim <- factor(data.df$Stim, levels = c("Ctrl", "SPP1"))

gg.ls <- list()
for (i in 1:length(sigmarker)) {
    data.df$Exp <- Exp.m[sigmarker[i], data.df$SeqNo]
    stat.test <- data.df %>%
        t_test(Exp ~ Stim, paired = TRUE, alternative = "less") %>%
        add_significance()
    stat.test <- stat.test %>% add_xy_position(x = "Stim")
    gg.ls[[i]] <- ggplot(data.df, aes(x = Stim, y = Exp)) +
        geom_boxplot() +
        ggtitle(sigmarker[i]) +
        geom_line(aes(group = Donor, color = CT), position = position_dodge(0)) +
        geom_point(aes(fill = CT, group = Donor), size = 2, shape = 21, position = position_dodge(0)) +
        scale_fill_manual(values = CTColor) +
        scale_colour_manual(values = CTColor) +
        ylab(sigmarker[i]) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        stat_pvalue_manual(stat.test, label = "{p}{p.signif}") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
}
multiPage <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")
ggexport(multiPage, filename = sprintf("%s/Figure6F.pdf", figdir), width = 8, height = 6)


# Correlation of ST frequency
load("./data/ST/spotFreq.RData")
cancerType.color <- scPalette(length(unique(pct.df$cancerType)))
names(cancerType.color) <- unique(pct.df$cancerType)

gg.ls <- list()
for (i in 1:nrow(cellpair.df)) {
    gg.ls[[i]] <- ggscatter(pct.df, x = cellpair.df$clusterA[i], y = cellpair.df$clusterB[i], size = 2, shape = 21, color = "cancerType") +
        coord_trans(clip = "off") +
        theme(legend.position = "bottom", legend.title = element_blank()) +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), method = "spearman") +
        scale_color_manual(values = cancerType.color)
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = length(gg.ls) / 2, common.legend = TRUE, legend = "right")
ggexport(multi.page, filename = sprintf("%s/FigureS6A.pdf", figdir), width = nrow(cellpair.df), height = 4)


# Cell proportion in different groups
TwoGroup.pct <- function(data.df, thresh = 15) {
    distance <- data.df$`Distance (μm)`
    return(c(sum(distance <= thresh) / length(distance), sum(distance > thresh) / length(distance)))
}
data.df <- readRDS("./data/mIHC_nearest_neighbor.csv")
data.ls <- split(data.df, f = data.df$`Analysis Region`)
distPct <- lapply(data.ls, TwoGroup.pct, thresh = 15)
distPct.df <- do.call(rbind, distPct)
colnames(distPct.df) <- c("Close to LRRC15-fib cells", "Not close to LRRC15-fib cells")
distPct.df <- as.data.frame(distPct.df) %>% rownames_to_column(var = "region")
plotData <- reshape2::melt(distPct.df, id.vars = c("region"), measure.vars = setdiff(colnames(distPct.df), "region"), variable.name = "group", value.name = "percentage")
plotData$group <- factor(plotData$group, levels = rev(c("Close to LRRC15-fib cells", "Not close to LRRC15-fib cells")))

g <- ggbarplot(plotData, x = "group", y = "percentage", add = c("mean_se", "jitter"), alpha = 0.2, color = "group", fill = "group", palette = "jco", position = position_dodge(0.8)) +
    ylab("% of SPP1-Mac total") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black")) +
    stat_compare_means(comparisons = list(c("Close to LRRC15-fib cells", "Not close to LRRC15-fib cells")), method = "wilcox.test", hide.ns = TRUE)
ggsave(sprintf("%s/FigureS6C.pdf", figdir), g, width = 3.5, height = 5)

regions <- list.files("./data/mIHC_regions/")
data.ls <- list()
pct.ls <- list()
for (i in 1:length(regions)) {
    LRRC15neg <- read_csv(sprintf("%s/LRRC15-FAP+ to SPP1+CD68 distance.csv", regions[i])) %>% pull(`Distance (μm)`)
    LRRC15pos <- read_csv(sprintf("%s/LRRC15+FAP+ to SPP1+CD68 distance.csv", regions[i])) %>% pull(`Distance (μm)`)
    pct.ls[[regions[i]]] <- data.frame(LRRC15neg = sum(LRRC15neg < 15) / length(LRRC15neg), LRRC15pos = sum(LRRC15pos < 15) / length(LRRC15pos), region = regions[i])
}

group.color <- c("#0073C2FF", "#EFC000FF")
names(group.color) <- c("FAP+LRRC15- Fib", "FAP+LRRC15+ Fib")

data.df <- do.call(rbind, pct.ls)
plotData <- reshape2::melt(data.df, id.vars = c("region"), measure.vars = setdiff(colnames(data.df), "region"), variable.name = "group", value.name = "percentage")
plotData$group <- factor(plotData$group, levels = c("FAP+LRRC15- Fib", "FAP+LRRC15+ Fib"))
plotData$percentage <- plotData$percentage * 100

g <- ggplot(plotData, aes(x = group, y = percentage, group = region)) +
    geom_line(size = 0.1, color = "grey") +
    geom_point(aes(color = group), size = 1, shape = 21) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
    ) +
    ylab("% of Fib close to SPP1-Mac") +
    xlab("") +
    scale_color_manual(values = group.color) +
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test", hide.ns = FALSE, paired = TRUE, label.y = 30)
ggsave(sprintf("%s/FigureS6D.pdf", figdir), g, width = 4, height = 3.5)

data.df <- readRDS("./data/ST_fib_group_pct.rds")
plotData <- reshape2::melt(data.df, id.vars = c("slide", "cancerType"), measure.vars = c("LRRC15neg", "LRRC15pos"), variable.name = "group", value.name = "pct")
plotData$group <- plyr::mapvalues(plotData$group, from = c("LRRC15neg", "LRRC15pos"), to = c("FAP+LRRC15- Fib", "FAP+LRRC15+ Fib"))
plotData$group <- factor(plotData$group, levels = c("FAP+LRRC15- Fib", "FAP+LRRC15+ Fib"))
group.color <- c("#0073C2FF", "#EFC000FF")
names(group.color) <- c("FAP+LRRC15- Fib", "FAP+LRRC15+ Fib")

g <- ggplot(plotData, aes(x = group, y = pct, color = group)) +
    geom_quasirandom(size = 0.5, shape = 21) +
    stat_summary(
        fun = "mean",
        geom = "crossbar",
        width = 0.5, linewidth = 0.15,
        colour = "black"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 8)
    ) +
    ylab("% of Fib close to SPP1-Mac") +
    xlab("") +
    facet_grid(. ~ cancerType, scales = "free", space = "free") +
    scale_color_manual(values = group.color) +
    stat_compare_means(aes(group = group), label = "p.format", method = "wilcox.test", hide.ns = FALSE, paired = TRUE, label.y = 0.9, size = 3)
ggsave(sprintf("%s/FigureS6E.pdf", figdir), g, width = 5, height = 3.1)


# Dot plot of ligand-receptor interactions
load("./data/CPDB_tumor_cellpairs.RData")
ctpair.df <- read_csv("./data/celltype_pair.csv") %>% arrange(ct1, ct2)
ctpair.df$ctpair <- paste0(ctpair.df$ct1, "-", ctpair.df$ct2)
plotData <- allcpdbF.df %>% filter(clusterA %in% ctpair.df$ct1, clusterB %in% ctpair.df$ct2, lrpair_order %in% genepairs_pick)
plotData$ctpair <- paste0(plotData$clusterA, "-", plotData$clusterB)
plotData <- plotData %>% filter(ctpair %in% ctpair.df$ctpair)
plotData$lrpair_order <- factor(plotData$lrpair_order, levels = rev(genepairs_pick))
plotData$ctpair <- factor(plotData$ctpair, levels = unique(ctpair.df$ctpair))
plotData$lr_means[plotData$lr_means == 0] <- 0.001
plotData$mean <- log2(plotData$lr_means)
plotData$mean <- MinMax(plotData$mean, max = 2, min = -4)
my_palette <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n = 20)
g <- ggplot(plotData, aes(x = lrpair_order, y = ctpair)) +
    geom_point(aes(size = -log10(cellphone_pvals + 0.0001), color = mean)) +
    scale_color_gradientn("lr means", colors = my_palette) +
    theme_bw() +
    theme(
        strip.text = element_text(size = 15, face = "bold"),
        panel.background = element_rect(colour = "black", fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "lightgray", size = 0.1, linetype = 2),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
    ) +
    rotate_x_text()

ggsave(sprintf("%s/FigureS6G", figdir), g, width = 14, height = 8)
