library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggrastr)
library(dendextend)
library(sscClust)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(grDevices)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(survival)
library(eks)
library(sscVis)


# ------------- Figure 4 --------------
run.cutree <- function(dat, method.hclust = "ward.D2", method.distance = "spearman", k = NULL, ...) {
    ret <- list()
    branch <- FALSE
    obj.hclust <- NULL
    if (method.distance == "spearman" || method.distance == "pearson") {
        tryCatch(
            {
                obj.distM <- as.dist(1 - cor.BLAS((dat), method = method.distance, nthreads = 1))
                obj.hclust <- stats::hclust(obj.distM, method = method.hclust)
            },
            error = function(e) {
                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n")
            }
        )
    } else if (method.distance == "cosine") {
        tryCatch(
            {
                sim <- dat / sqrt(rowSums(dat * dat))
                sim <- sim %*% t(sim)
                obj.distM <- as.dist(1 - sim)
                obj.hclust <- stats::hclust(obj.distM, method = method.hclust)
            },
            error = function(e) {
                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n")
            }
        )
    }
    if (is.null(obj.hclust)) {
        obj.distM <- stats::dist(dat)
        obj.hclust <- stats::hclust(obj.distM, method.hclust)
    }
    obj.dend <- as.dendrogram(obj.hclust)
    cluster <- cutree(obj.hclust, k = k, ...)
    colSet.cls <- auto.colSet(length(unique(cluster)), "Paired")
    branch <- dendextend::color_branches(obj.dend, clusters = cluster[order.dendrogram(obj.dend)], col = colSet.cls)
    branch <- dendextend::set(branch, "branches_lwd", 1.5)
    ret[["hclust"]] <- obj.hclust
    ret[["dist"]] <- obj.distM
    ret[["cluster"]] <- cluster
    ret[["branch"]] <- branch
    return(ret)
}

c68 <- c(
    "#ed3d30", "#55d151", "#0e09b7", "#f69a95", "#aa329a", "#1d96a9", "#EDE816", "#48dc88", "#b15b0a", "#ac66d2", "#dc1c6f",
    "#CB00FC", "#af3292", "#5da3ac", "#c6ab48", "#0c9986", "#2f52d1", "#aa152e", "#729e0c", "#be5793", "#3cc0c0", "#9ebe71", "#9a221c",
    "#FD9A97", "#0c9c54", "#C0AFE0", "#D283FA", "#A30097", "#FDC44D", "#737D58", "#9e2d71", "#5A2AB4", "#0D6094", "#9F4700", "#16A795",
    "#BF658B", "#59d185", "#FED8A3", "#7E9800", "#F82263", "#A35FFD", "#2285FE", "#F65FFE", "#75ADCD", "#959099", "#FB5A8A", "#BD856E",
    "#92628B", "#C5D0B7", "#C4B43B", "#E89B62", "#0DA8E7", "#822E35", "#5DE9EB", "#FEC4EC", "#8F00B3", "#1CB232", "#F668DE", "#9878CE",
    "#A61C79", "#FC7216", "#76AE86", "#F87AA3", "#E3DDEF", "#CEE769", "#E7ABFD", "#53475A", "#0D56B2"
)

# Volcano plot for marker genes
load("./data/c05_c03_DEgenes.RData")
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
ggsave(sprintf("%s/Figure4A.pdf", figdir), g1, width = 6, height = 6)


# Barplot of enriched pathways
geneset.df <- readRDS("./data/c05vsc03_enriched_Term.rds")
g <- ggplot(geneset.df, aes(x = label, y = logFDR, fill = cluster)) +
    geom_bar(stat = "identity", width = 0.6) + # geom_boxplot(alpha=.9) +
    scale_fill_manual(values = c("#49A75A", "#D77F32")) +
    facet_grid(rows = vars(cluster), scales = "free", space = "free") +
    theme_bw() +
    xlab("") +
    ylab("-log10(FDR)") +
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
ggsave(sprintf("%s/Figure4B.pdf", figdir), g, width = 7, height = 3)


# Boxplot of fibroblast proportions
library(data.table)
load("./data/fib_freq_sample.RData")
clusterVar <- c("c05", "c03", "c04")
pcatVar <- c("Normal", "Tumor")
tissueVar <- c("Colon", "Stomach", "Ovary", "Lung", "Kidney", "Pancreas")
min.NTotal <- 30

fib_freq_sub.tb <- fib_freq.tb[cluster %in% clusterVar & tissue %in% tissueVar & pcat %in% pcatVar, ]
fib_freq_sub.tb <- fib_freq_sub.tb[NTotal >= min.NTotal, ]
fib_freq_sub.tb[, pcat := factor(pcat, levels = pcatVar)]

gg.ls <- list()
for (i in 1:length(tissueVar)) {
    fib_freq_tissue <- fib_freq_sub.tb[tissue == tissueVar[i], ]
    fib_freq_tissue$cluster <- factor(fib_freq_tissue$cluster, levels = clusterVar)
    gg.ls[[i]] <- ggboxplot(fib_freq_tissue,
        x = "cluster",
        y = "freq",
        add = "jitter",
        outlier.shape = NA,
        color = "pcat",
        size = 0.5,
        add.params = list(size = 1.5, alpha = 0.8)
    ) +
        scale_color_manual(values = pcat.color) +
        rotate_x_text(45) +
        labs(x = "", y = "Percentage (%)", title = tissueVar[i]) +
        theme(legend.position = "None", plot.title = element_text(hjust = 0.5))
    gg.ls[[i]] <- g
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = length(tissueVar) / 2, common.legend = FALSE)
ggexport(multi.page, filename = sprintf("%s/Figure4C.pdf", figdir), width = length(tissueVar), height = 5)

cellpct.df <- readRDS("./data/staining_cell_pct.rds")
g <- ggplot(cellpct.df, aes(x = pcat, y = clust_prop, color = pcat)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = pcat), alpha = 1, size = 1, position = position_jitter(w = 0.1, h = 0.1)) +
    facet_grid(cluster ~ cancerType, scales = "free_y", space = "fixed") +
    theme_classic() +
    ylab("Cell Proportion (%)") +
    xlab("") +
    scale_color_manual(values = pcat.color) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
    ) +
    stat_compare_means(comparisons = list(c("Normal", "Nadj"), c("Nadj", "Tumor")), label = "p.format", method = "t.test", hide.ns = TRUE, size = 2.5)
ggsave(sprintf("%s/Figure4E.pdf", figdir), g, width = 5, height = 4)


# Heatmap of fibroblast frequency
load("./data/fib_freq_tumor.RData")
data.T <- metaData %>%
    filter(pcat == "Tumor") %>%
    select(sample, cluster) %>%
    group_by(sample, cluster) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(clust_total = sum(count)) %>%
    dplyr::mutate(clust_prop = count / clust_total)
plot.data <- reshape2::dcast(data.T, cluster ~ sample, value.var = "clust_prop") %>%
    column_to_rownames(var = "cluster") %>%
    as.matrix()

cancerType.color <- c("#725e82", "#8c4356", "#db5a6b", "#21a675", "#d9b611", "#e29c45", "#4c8dae", "#003371", "#a1afc9", "#88ada6", "#c93756")
names(cancerType.color) <- c("Ovary", "Kidney", "Breast", "Synovium", "Colon", "Stomach", "Lung", "Liver", "Skin", "HNSC", "Pancreas")

mapping.tb <- metaData %>%
    dplyr::filter(pcat == "Tumor") %>%
    select(sample, tissue) %>%
    distinct()
mapping.vec <- structure(mapping.tb$tissue, names = mapping.tb$sample)

ha.col <- HeatmapAnnotation(
    cancerType = mapping.vec[colnames(plot.data)],
    FibClusters = anno_barplot(t(plot.data),
        gp = gpar(fill = cluster.color, col = NA),
        bar_width = 1, height = unit(2.5, "cm")
    ),
    col = cancerType.color,
    show_legend = T,
    simple_anno_size = unit(1, "cm")
)

obj.hclust.col <- run.cutree(t(plot.data), k = 8, method.distance = "")
obj.hclust.row <- run.cutree(plot.data, k = 2)
obj.hclust.col$branch <- dendextend::set(obj.hclust.col$branch, "branches_lwd", 2)
obj.hclust.row$branch <- dendextend::set(obj.hclust.row$branch, "branches_lwd", 2.5)
plot.data.bin <- floor(floor(plot.data * 100) / 5)
plot.data.bin[plot.data.bin > 9] <- 9
plot.data.bin[1:3, 1:4]

ha.row <- rowAnnotation(overall = anno_boxplot(plot.data,
    outline = FALSE, width = unit(4, "cm"),
    gp = gpar(fill = cluster.color[rownames(plot.data)])
))

sscVis::plotMatrix.simple(
    plot.data.bin,
    col.ht = rev(structure(viridis::viridis(10), names = 0:9)),
    par.legend = list(ncol = 1, labels = rev(sprintf("%s%%~%s%%", 5 * (0:9), c(5 * (1:9), 100)))),
    show.dendrogram = T,
    out.prefix = sprintf("%s/Figure4F", fig.dir),
    show.number = F,
    clust.row = obj.hclust.row$branch,
    clust.column = obj.hclust.col$branch,
    exp.name = expression(Freq),
    top_annotation = ha.col,
    column.split = 8,
    par.heatmap = list(
        row_names_gp = gpar(fontsize = 12),
        right_annotation = ha.row,
        column_gap = unit(0.8, "mm"), border = FALSE,
        row_dend_width = unit(1.5, "cm"),
        column_dend_height = unit(1.5, "cm")
    ),
    pdf.width = 16, pdf.height = 10
)

plot.data.long <- reshape2::melt(as.data.frame(plot.data), id.vars = c("cluster"), measure.vars = colnames(plot.data), variable.name = "sample", value.name = "pct")
plot.data.long$cluster <- factor(plot.data.long$cluster, levels = cluster.order)

g <- ggplot(plot.data.long, aes(cluster, pct, color = cluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = cluster), alpha = 0.3, size = 1) +
    scale_color_manual(values = cluster.color) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    labs(fill = NULL) +
    ylab("") +
    ylim(c(0, 0.8)) +
    theme_light() +
    theme(
        strip.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.2, "cm"),
        legend.background = element_rect(color = "black", size = 0.4),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 7, color = "black", face = "bold", angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_blank()
    )
ggsave(sprintf("%s/Figure4F_boxplot.pdf", figdir), g, width = 6, height = 2)


# Survival forest plot
do.cox.group <- function(df, group.order) {
    cox <- c()
    df$stromalType <- factor(df$stromalType, levels = group.order)
    if (length(unique(df$gender)) > 1) {
        cox <- coxph(Surv(OS.time, OS) ~ stromalType + gender + age + stage, data = df)
    } else {
        cox <- coxph(Surv(OS.time, OS) ~ stromalType + age + stage, data = df)
    }
    cox.summary <- summary(cox)
    nvar <- length(unique(df$stromalType)) - 1
    HR <- round(cox.summary$conf.int[1:nvar, 1], 2)
    HR.lower <- round(cox.summary$conf.int[1:nvar, 3], 2)
    HR.upper <- round(cox.summary$conf.int[1:nvar, 4], 2)
    HR.range <- sprintf("(%.1f-%.1f)", HR.lower, HR.upper)
    coef <- cox.summary$coefficients[1:nvar, 1]
    coef.se <- cox.summary$coefficients[1:nvar, 3]
    Pval <- round(cox.summary$coefficients[1:nvar, 5], 4)
    stromalType <- gsub("stromalType", "", rownames(cox.summary$conf.int)[1:nvar])
    return(data.frame(stromalType = stromalType, HR = HR, HR.range = HR.range, coef = coef, coef.se = coef.se, Pval = Pval))
}

survivalInfo <- readRDS("./data/TCGA_fib_group.rds")
plot.data <- survivalInfo[, c("gender", "age", "stage", "OS.time", "OS", "cancerType", "stromalType")]
HR.df <- ddply(plot.data, c("cancerType"), do.cox.group, group.order = c("myoFloprogFhi", "myoFhiprogFlo"))
HR.df <- HR.df %>%
    group_by(stromalType) %>%
    mutate(adj.Pval = p.adjust(Pval, method = "BH")) %>%
    arrange(HR) %>%
    filter(stromalType == "myoFloprogFhi")

meta <- meta::metagen(
    TE = HR.df$coef, seTE = HR.df$coef.se, studlab = HR.df$cancerType,
    comb.fixed = F, comb.random = T, prediction = F, sm = "HR"
)
pdf(sprintf("%s/Figure4H.pdf", figdir), width = 14, height = 10)
meta::forest(meta,
    test.overall.random = T, digits.pval = 4,
    colgap.forest.left = "5cm", zero.pval = T
)
dev.off()

seu <- readRDS("./data/all_fib_cells.rds")
genes <- c("PI16", "MFAP5")
geneData <- t(as.matrix(seu@assays$RNA@data[genes, ]))
geneData <- as.data.frame(geneData) %>% rownames_to_column(var = "cellName")
geneData$tissue_pcat <- seu$tissue_pcat[match(geneData$cellName, seu$cellName)]

tissue_pcat.order <- seu@meta.data %>%
    filter(tissue %in% c("Stomach", "Colon", "Lung", "Kidney"), pcat %in% c("Healthy", "Normal", "Tumor")) %>%
    arrange(tissue, pcat) %>%
    pull(tissue_pcat) %>%
    unique()

gg.ls <- list()
for (i in 1:length(tissue_pcat.order)) {
    gg.ls[[i]] <- DoISFACSPlot(
        plot.data = geneData %>% dplyr::filter(tissue_pcat %in% tissue_pcat.order[i]),
        genes = genes, xlims = c(-0.5, 6), ylims = c(-0.5, 5), color.by = "None", color.panel = "YlOrRd",
        dot.size = 0.1, x.cutoff = 1, y.cutoff = 1, bins = 30, bandwidth = 2, nColors = 100, title = tissue_pcat.order[i]
    )
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 4, ncol = 3, common.legend = FALSE)
ggexport(multi.page, filename = sprintf("%s/FigureS4A.pdf", figdir), width = 15, height = 20)


# Scatter plot of DE genes
load("./data/HealthyvsNormal_avgExp_genes.RData")
markers_DEs.ls <- readRDS("./data/HealthyvsNormal_DE.rds")
tissues.order <- c("Colon", "Kidney", "Lung", "Stomach")
category.colors <- c("#E41A1C", "#E8E8E8", "#377EB8")
names(category.colors) <- c("Normal+", "NoDiff", "Healthy+")

gg.ls <- list()
for (i in tissues.order) {
    tissue <- tissues.order[i]
    avgExp_genes <- avgExp[, c(paste0(tissue, "_Healthy"), paste0(tissue, "_Normal"))] %>% rownames_to_column(var = "gene")
    colnames(avgExp_genes) <- c("gene", "Healthy", "Normal")
    avgExp_genes$Diff <- avgExp_genes$Normal - avgExp_genes$Healthy
    avgExp_genes$Group <- "NoDiff"
    avgExp_genes$Group[which(avgExp_genes$gene %in% c(markers_DEs.ls[[tissue]] %>% filter(avg_log2FC > 0.25, p_val_adj < 0.01) %>% pull(gene)))] <- "Normal+"
    avgExp_genes$Group[which(avgExp_genes$gene %in% c(markers_DEs.ls[[tissue]] %>% filter(avg_log2FC < -0.25, p_val_adj < 0.01) %>% pull(gene)))] <- "Healthy+"
    avgExp_genes$Group <- factor(avgExp_genes$Group, levels = c("Normal+", "Healthy+", "NoDiff"))
    labLabel <- avgExp_genes %>% filter(gene %in% upgenes, Group == "Normal+")
    gg.ls[[i]] <- ggplot(avgExp_genes, aes(x = Healthy, y = Normal, colour = Group)) +
        geom_point_rast(size = 0.5) +
        xlim(c(0, 6)) +
        ylim(c(0, 6)) +
        geom_text_repel(
            data = labLabel, aes(x = Healthy, y = Normal, label = gene), colour = "black", box.padding = 0.5,
            segment.size = 0.1, size = 4.3, fontface = "italic"
        ) +
        scale_colour_manual(values = category.colors) +
        labs(title = "", x = "Healthy", y = "Adj-Normal") +
        ggtitle(tissue) +
        theme_light() +
        theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 15),
            plot.title = element_text(hjust = 0.5, size = 20),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            legend.title = element_blank(), legend.position = "None"
        )
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 2, ncol = 2)
ggexport(multi.page, filename = sprintf("%s/FigureS4B.pdf", figdir), width = 10, height = 10)

# Heatmap of ligand-induced gene expression
load("./data/Stim_induced_DE_genes.RData")
Donor.v <- metaData %>%
    filter(SeqNo %in% colnames(normData)) %>%
    pull(Donor) %>%
    unique()
data.ls <- list()
for (i in 1:length(Donor.v)) {
    samplesTest <- metaData %>%
        filter(Donor == Donor.v[i], Stim != "Ctrl") %>%
        pull(SeqNo)
    samplesCtrl <- metaData %>%
        filter(Donor == Donor.v[i], Stim == "Ctrl") %>%
        pull(SeqNo)
    data.ls[[Donor.v[i]]] <- normData[, samplesTest] - normData[, samplesCtrl]
}
dataStim.m <- do.call(cbind, data.ls)

dataStim.df <- as.data.frame(dataStim.m) %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(!gene, names_to = "SeqNo", values_to = "zscore")
dataStim.df <- left_join(dataStim.df, metaData, by = "SeqNo")

dataStim_avgExp.df <- dataStim.df %>%
    filter(gene %in% DEgenes) %>%
    group_by(gene, Stim) %>%
    summarize(mean = mean(zscore))

dataStim_avgExp.m <- reshape2::dcast(dataStim_avgExp.df, gene ~ Stim, value.var = "mean") %>%
    column_to_rownames(var = "gene") %>%
    as.matrix()

expr <- as.matrix(dataStim_avgExp.m)
mat2 <- t(scale(t(expr)))
plotData <- MinMax(mat2, max = 2, min = (-1) * 2)

MarkGenes <- c(
    "ADH1B", "MFAP5", "ACKR3", "OGN", "CLU", "CXCL12", "MFAP4", "COL15A1", "NFIX", "CFD", "C7", "APOD", "PI16", "DPT",
    "CD55", "NFKBIA", "SOD2", "CXCL5", "IL1B", "CCL2", "CXCL8", "IL6", "CXCL1", "CXCL6"
)
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
cols <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "RdBu")))(100)
breaks <- seq(-2, 2, length = 100)
col1 <- colorRamp2(breaks, cols)

pdf(sprintf("%s/FigureS4F", figdir), width = 3, height = 5)
Heatmap(as.matrix(plotData),
    name = "z-score",
    right_annotation = rowHa,
    cluster_rows = T, cluster_columns = T,
    col = col1, border = F,
    show_row_names = FALSE, show_column_names = TRUE,
    use_raster = TRUE, raster_resize_mat = max
)
dev.off()


# Violin plot of minimum distances
load("./data/ST/sampleFiles.RData")
for (dataFile in dataList) {
    c04CancerDist <- fread(paste0("./data/ST/", dataFile, "_c04_cancer_dist.csv"))
    c05CancerDist <- fread(paste0("./data/ST/", dataFile, "_c05_cancer_dist.csv"))
    data <- data.frame(cat = "c04", values = c04CancerDist$min_val)
    data <- rbind(data, data.frame(cat = "c05", values = c05CancerDist$min_val))
    print(paste0(dataFile, ": Mean value of c04, ", mean(c04CancerDist$min_val), " and c05,", mean(c05CancerDist$min_val)))
    g <- ggplot(data, aes(x = cat, y = values, fill = cat, color = cat)) +
        geom_violin(scale = "width") +
        geom_boxplot(outlier.shape = NA, width = 0.05, color = "grey") +
        scale_fill_manual(values = c("#87638F", "#D77F32")) +
        scale_color_manual(values = c("#87638F", "#D77F32")) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(x = "", y = "Minimum distance between cancer cells") +
        stat_compare_means()
    ggsave(sprintf("%s/FigureS4G_%s.pdf", figdir, dataFile), g, width = 2, height = 3)
}


# Cluster proportion in different ST spot groups
load("./data/ST/ST_tumor_neighbor.RData")
score_df$group <- factor(score_df$group, levels = c("tumor", "juxtatumor", "distant"))
CT.color <- scPalette(length(unique(score_df$cancerType)))
names(CT.color) <- unique(score_df$cancerType)

gg.ls <- list()
for (i in 1:length(cluster.order)) {
    gg.ls[[cluster.order[i]]] <- ggplot(score_df %>% filter(cluster == cluster.order[i]), aes(x = group, y = score, color = cancerType)) +
        geom_quasirandom(size = 0.5, shape = 21) +
        stat_summary(
            fun = "mean",
            geom = "crossbar",
            width = 0.5, linewidth = 0.15,
            colour = "black"
        ) +
        theme_classic() +
        ylab(cluster.order[i]) +
        scale_color_manual(values = CT.color) +
        theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
        ) +
        stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test", ref.group = "tumor", hide.ns = FALSE, size = 2.5)
}
multi.page <- ggarrange(plotlist = gg.ls, nrow = 1, ncol = length(cluster.order), common.legend = TRUE, legend = "right")
ggexport(multi.page, filename = sprintf("%s/FigureS4H.pdf", figdir), width = 3 * length(cluster.order), height = 5)
