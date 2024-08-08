library(dplyr)
library(Seurat)
library(grDevices)
library(ggplot2)
library(ggrastr)
library(ggthemes)
library(cowplot)
library(ggpubr)
library(progeny)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(survival)
library(maxstat)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(data.table)
library(RColorBrewer)
library(nichenetr)
library(tidyverse)
library(magrittr)
library(sscVis)
library(R.utils)
library(plyr)
library(grid)

RhpcBLASctl::omp_set_num_threads(1)
doParallel::registerDoParallel(cores = 12)

# ------------- Figure 3 --------------
# Dotplot for marker genes
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
`%||%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
        return(lhs)
    } else {
        return(rhs)
    }
}
scale.func <- switch(
    EXPR = scale.by,
    "size" = scale_size,
    "radius" = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
)
myDotPlot <- function(object,
                      assay = NULL,
                      features,
                      cols = c("lightgrey", "red"),
                      col.min = -2.5,
                      col.max = 2.5,
                      dot.min = 0,
                      dot.scale = 6,
                      group.by = NULL,
                      split.by = NULL,
                      scale.by = "radius",
                      scale.min = NA,
                      scale.max = NA, scale.pct = NULL,
                      keepFeaturesOrder = F, stroke = 0.5, color_pal = SpatialColors(20), returnDF = FALSE) {
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    data.features <- FetchData(object = object, vars = features)
    data.features <- data.features[colnames(object), ]
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    } else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }

    data.plot <- lapply(
        X = unique(x = data.features$id),
        FUN = function(ident) {
            data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
            avg.exp <- apply(
                X = data.use,
                MARGIN = 2,
                FUN = function(x) {
                    return(mean(x = expm1(x = x)))
                }
            )
            pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, threshold = 0)
            return(list(avg.exp = avg.exp, pct.exp = pct.exp))
        }
    )
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(
        X = names(x = data.plot),
        FUN = function(x) {
            data.use <- as.data.frame(x = data.plot[[x]])
            data.use$features.plot <- rownames(x = data.use)
            data.use$id <- x
            return(data.use)
        }
    )
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(
        X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
            data.use <- scale(x = data.use)
            data.use <- MinMax(data = data.use, min = col.min, max = col.max)
            return(data.use)
        }
    )
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(
            X = strsplit(x = as.character(x = data.plot$id), split = "_"),
            FUN = "[[",
            FUN.VALUE = character(length = 1L),
            2
        )
        data.plot$colors <- mapply(
            FUN = function(color, value) {
                return(colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(20)[value])
            },
            color = cols[splits.use],
            value = avg.exp.scaled
        )
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    geneGroup <- data.plot %>%
        group_by(features.plot) %>%
        filter(avg.exp == max(avg.exp)) %>%
        dplyr::select(features.plot, id) %>%
        as.data.frame() %>%
        dplyr::distinct(features.plot, .keep_all = TRUE) %>%
        dplyr::select(features.plot, group = id)

    data.plotT <- data.plot %>%
        dplyr::left_join(geneGroup, by = c("features.plot" = "features.plot")) %>%
        dplyr::group_by(group) %>%
        dplyr::arrange(desc(pct.exp), .by_group = TRUE)

    if (keepFeaturesOrder) {
        data.plotT$group <- levels(Idents(object = object))[1]
        data.plotT$features.plot <- factor(data.plotT$features.plot, levels = unique(features))
    } else {
        data.plotT$features.plot <- factor(data.plotT$features.plot, levels = unique(data.plotT$features.plot))
    }
    if (!is.null(scale.pct)) {
        data.plotT$pct.exp[data.plotT$pct.exp > scale.pct] <- scale.pct
    }

    plot <- ggplot(data = data.plotT, mapping = aes_string(x = "features.plot", y = "id")) +
        geom_point_rast(mapping = aes_string(size = "pct.exp", fill = color.by), stroke = stroke, shape = 21, color = "black") +
        scale_x_discrete(position = "top") +
        facet_grid(. ~ group, scales = "free", space = "free") +
        scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        guides(size = guide_legend(title = "Percent Expressed")) +
        labs(
            x = "Features",
            y = ifelse(test = is.null(x = split.by), yes = "Identity", no = "Split Identity")
        ) +
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
    plot <- plot + scale_fill_gradientn(colours = color_pal)
    if (returnDF) {
        return(data.plotT)
    } else {
        return(plot)
    }
}

load("./data/myoF_signature_genes.RData")
seu <- readRDS("./data/all_fib_cells.rds")
Idents(seu) <- "cluster"
seu.myoF <- subset(seu, idents = myoF.clusters)
Idents(seu.myoF) <- "cluster"
g <- myDotPlot(seu.myoF, features = myoF_markers, keepFeaturesOrder = FALSE, stroke = 0.5, color_pal = brewer.pal(8, "OrRd"))
ggsave(sprintf("%s/Figure3A.pdf", figdir), g, height = 2, width = 10, limitsize = FALSE)


# PROGENy pathway activity
seu <- readRDS("./data/all_fib_cells.rds")
metaData <- seu@meta.data
clusterColor <- readRDS("./data/cluster_colors.rds")

CellsClusters <- data.frame(
    Cell = colnames(seu),
    CellType = as.character(seu@meta.data$cluster),
    stringsAsFactors = FALSE
)

seu <- progeny(seu, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
seu <- Seurat::ScaleData(seu, assay = "progeny")
progeny_scores_df <- as.data.frame(t(GetAssayData(seu, slot = "scale.data", assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell)

progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
sub_progeny_scores <- summarized_progeny_scores %>% filter(Pathway %in% c("TGFb", "WNT", "PI3K", "Hypoxia", "TNFa", "NFkB"))

cluster.lineages <- c("c04", "c11", "c16", "c19", "c09", "c17", "c06", "c08", "c05", "c03")
sub_progeny_scores <- sub_progeny_scores %>% filter(CellType %in% cluster.lineages)

sub_progeny_scores$CellType <- factor(sub_progeny_scores$CellType, levels = rev(cluster.lineages))
sub_progeny_scores$Pathway <- factor(sub_progeny_scores$Pathway, levels = c("TGFb", "WNT", "PI3K", "Hypoxia", "TNFa", "NFkB"))
sub_progeny_scores <- sub_progeny_scores %>% arrange(Pathway)
sub_progeny_scores$group <- "Terminal"
sub_progeny_scores$group[sub_progeny_scores$CellType %in% c("c11", "c16", "c19")] <- "Transit1"
sub_progeny_scores$group[sub_progeny_scores$CellType %in% c("c09", "c17", "c06")] <- "Transit2"
sub_progeny_scores$group[sub_progeny_scores$CellType %in% c("c08", "c05", "c03")] <- "Progenitor"
sub_progeny_scores$group <- factor(sub_progeny_scores$group, levels = c("Terminal", "Transit1", "Transit2", "Progenitor"))

groupMean <- sub_progeny_scores %>%
    dplyr::group_by(Pathway) %>%
    dplyr::mutate(MeanV = mean(avg)) %>%
    dplyr::distinct(Pathway, .keep_all = TRUE) %>%
    dplyr::select(Pathway, MeanV)

g1 <- ggplot(sub_progeny_scores, aes(x = CellType, y = avg)) +
    geom_point(size = 2) +
    geom_hline(data = groupMean, aes(yintercept = MeanV), linetype = "dashed") +
    coord_flip() +
    theme_gray(base_size = 12) +
    theme(axis.text.x = element_text(size = 7)) +
    facet_grid(group ~ Pathway, scales = "free")
ggsave(sprintf("%s/Figure3D.pdf", figdir), g1, width = 6, height = 4, useDingbats = F)


# Survival forest plot

survivalInfo <- readRDS("./data/TCGA_cluster_score.rds")
do.cox <- function(df, type) {
    cox <- c()
    if (type == "maxstat") {
        aa_try <- try(aa <- maxstat.test(Surv(OS.time, OS) ~ Score, data = df, smethod = "LogRank"), silent = TRUE)
        if (!is(aa_try, "try-error") & length(table(df$Score >= aa$estimate)) == 2) {
            df$stromalType <- ifelse(df$Score >= aa_try$estimate, "High", "Low")
            print("Use maxstat.test to get cutting point")
        } else {
            df$stromalType <- ifelse(df$Score >= median(df$Score, na.rm = TRUE), "High", "Low")
            print("Use median as cutting point")
        }
    } else (type == "median") {
        df$stromalType <- ifelse(df$Score >= median(df$Score, na.rm = TRUE), "High", "Low")
    }

    df$stromalType <- factor(df$stromalType, levels = c("Low", "High"))
    if (length(unique(df$gender)) > 1) {
        cox <- coxph(Surv(OS.time, OS) ~ stromalType + gender + age + stage, data = df)
    } else {
        cox <- coxph(Surv(OS.time, OS) ~ stromalType + age + stage, data = df)
    }
    cox.summary <- summary(cox)
    nvar <- length(unique(df$stromalType)) - 1
    #
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
plot.data <- survivalInfo[, c("gender", "age", "stage", "OS.time", "OS", "cancerType", myoF.clusters)]
HR.ls <- list()
for (i in 1:length(myoF.clusters)) {
    plot.data$Score <- plot.data[, myoF.clusters[i]]
    HR.df <- ddply(plot.data, c("cancerType"), do.cox, type = "maxstat")
    HR.df <- HR.df %>%
        group_by(stromalType) %>%
        mutate(adj.Pval = p.adjust(Pval, method = "BH")) %>%
        arrange(HR)
    HR.df$cluster <- myoF.clusters[i]
    HR.ls[[i]] <- HR.df
}
HR_all <- do.call(rbind, HR.ls)

plyr::d_ply(HR_all, .(cluster), function(.df) {
    meta <- meta::metagen(
        TE = .df$coef, seTE = .df$coef.se, studlab = .df$cancerType,
        comb.fixed = F, comb.random = T, prediction = F, sm = "HR"
    )
    pdf(sprintf("%s/Figure3E_%s.pdf", figdir, unique(.df$cluster)), width = 14, height = 10)
    meta::forest(meta,
        test.overall.random = T, digits.pval = 4,
        colgap.forest.left = "5cm", zero.pval = T
    )
    dev.off()
}, .parallel = F)


# Boxplot of proportions
metadata_summary <- metaData %>%
    group_by(sample, cluster) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(clust_total = sum(count)) %>%
    dplyr::mutate(clust_prop = count / clust_total * 100)

metadata_summary <- left_join(metadata_summary, metaData[, c("sample", "tissue", "pcat")], by = "sample")
metadata_summary <- metadata_summary %>% filter(clust_total > 20)
metadata_summary$tissue_pcat <- paste0(metadata_summary$tissue, "_", metadata_summary$pcat)

tissue.order <- metadata_summary %>%
    filter(pcat == "Tumor", cluster == "c04") %>%
    group_by(tissue) %>%
    summarize(mean = median(clust_prop)) %>%
    arrange(-mean) %>%
    pull(tissue) %>%
    rev()
tissue.order <- c(tissue.order, "Synovium")
pcat.order <- c("Inflamed", "Tumor")
metadata_sub <- metadata_summary %>% filter(pcat %in% pcat.order, cluster %in% myoF.clusters)

g <- ggplot(metadata_sub, aes(tissue_pcat, clust_prop, fill = cluster)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = colSet$cluster.name.full) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    facet_grid(cluster ~ pcat, space = "free_x", scales = "free") +
    labs(fill = NULL) +
    ylab("Proportions %") +
    theme_light() +
    theme(
        strip.background = element_rect(colour = "black", fill = "white", size = 1.5, linetype = "solid"),
        strip.text.x = element_text(size = 10, color = "black", face = "bold"),
        strip.text.y = element_text(size = 10, color = "black", face = "bold"),
        legend.position = "top",
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.3, "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.title.y = element_text(face = "bold")
    ) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 0.5)))
ggsave(sprintf("%s/Figure3F.pdf", figdir), g, width = 4.1, height = 7)


# Heatmap of regulon activity score
load("./data/fib_regulons.RData")

auc_mtx <- readRDS("./data/fib_auc_mtx.rds")
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
topGenes <- unique(unlist(topregulons.ls[myoF.clusters]))
expr <- t(auc_mtx_avg.m)[topGenes, myoF.clusters]
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

pdf(sprintf("%s/Figure3G.pdf", figdir), width = 5, height = 5)
Heatmap(plotData,
    name = "z-score",
    top_annotation = ha,
    right_annotation = rowHa,
    cluster_rows = F, cluster_columns = F,
    col = col1, border = F,
    show_row_names = FALSE, show_column_names = TRUE
)
dev.off()

cluster.order <- sort(rownames(auc_mtx_avg.m))
avg.exp <- t(auc_mtx_avg.m)

ann_colors <- list(
    cluster = cluster.color
)
annCol <- data.frame(cluster = colnames(avg.exp))
rownames(annCol) <- colnames(avg.exp)

expr <- as.matrix(avg.exp[intersect(markerGenes, rownames(avg.exp)), ])
mat2 <- t(scale(t(expr)))
mat2 <- MinMax(mat2, max = 2, min = (-1) * 2)
colnames(mat2) <- colnames(expr)
plotData <- mat2[markerGenes, cluster.order]

pdf(sprintf("%s/Figure3H.pdf", figdir), width = 6, height = 3)
pheatmap(plotData,
    cluster_rows = TRUE, cluster_cols = TRUE,
    color = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(100)),
    breaks = seq(-2.5, 2.5, length = 101),
    annotation_col = annCol,
    annotation_colors = annColors,
    show_colnames = TRUE, show_rownames = TRUE,
    border_color = "white",
    border = "white",
    cluster_cols = TRUE,
    angle_col = 90, fontsize = 5
)
dev.off()



# Heatmap of regulatory potential for top-ranked ligands
weighted_networks <- readRDS("./data/Nichnet/weighted_networks.rds")
ligand_target_matrix <- readRDS("./data/Nichnet/ligand_target_matrix.rds")
lr_network <- readRDS("./data/Nichnet/lr_network.rds")
load("./data/fib_markers.RData")

getNichenet <- function(gene.oi, gene.bg, out.prefix) {
    str(ligand_target_matrix)
    lr_network_expressed <- lr_network %>% filter(to %in% gene.bg)
    potential_ligands <- lr_network_expressed %>%
        pull(from) %>%
        unique()
    ligand_activities <- predict_ligand_activities(
        geneset = gene.oi,
        background_expressed_genes = gene.bg,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
    )
    write_csv(ligand_activities, sprintf("%s.ligand_activities.csv", out.prefix))
    return(ligand_activities)
}

plotNichenet <- function(gene.oi, gene.bg, out.prefix, ligand_activities,
                         width.target = 7, width.pearson = 1.3, width.comb = 8, height = 7,
                         rel_widths = c(0.19, 0.81), comb.rel.height = c(10, 2),
                         n.top = 20, pearson.max = 0.2) {
    # ligand target
    best_upstream_ligands <- ligand_activities %>%
        top_n(n.top, pearson) %>%
        arrange(-pearson) %>%
        pull(test_ligand)
    active_ligand_target_links_df <- best_upstream_ligands %>%
        lapply(get_weighted_ligand_target_links,
            geneset = gene.oi,
            ligand_target_matrix = ligand_target_matrix, n = 250
        ) %>%
        bind_rows()
    active_ligand_target_links_df <- active_ligand_target_links_df %>% filter(!is.na(weight))
    active_ligand_target_links <- prepare_ligand_target_visualization(
        ligand_target_df = active_ligand_target_links_df,
        ligand_target_matrix = ligand_target_matrix,
        cutoff = 0.25
    )
    order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
    order_targets <- active_ligand_target_links_df$target %>% unique()
    order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
    active_ligand_target_links.debug <<- active_ligand_target_links
    order_targets.debug <<- order_targets
    order_ligands.debug <<- order_ligands
    vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
    colnames(vis_ligand_target) <- make.names(colnames(vis_ligand_target))
    gene.unexplained <- setdiff(gene.oi, colnames(active_ligand_target_links))
    vis_ligand_target.debug <<- vis_ligand_target
    vis_ligand_target[vis_ligand_target > 0.01] <- 0.01

    # ligand pearson
    ligand_pearson_matrix <- ligand_activities %>%
        select(pearson) %>%
        as.matrix() %>%
        magrittr::set_rownames(ligand_activities$test_ligand)

    vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, ] %>%
        as.matrix(ncol = 1) %>%
        magrittr::set_colnames("Pearson")
    vis_ligand_pearson.debug <<- vis_ligand_pearson
    vis_ligand_pearson[vis_ligand_pearson > pearson.max] <- pearson.max

    # ligand target
    lr_network_top <- lr_network %>%
        filter(from %in% best_upstream_ligands & to %in% gene.bg) %>%
        distinct(from, to)
    best_upstream_receptors <- lr_network_top %>%
        pull(to) %>%
        unique()
    lr_network_top_df <- weighted_networks$lr_sig %>%
        filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    lr_network_top_df <- lr_network_top_df %>% spread("from", "weight", fill = 0)
    lr_network_top_matrix <- lr_network_top_df %>%
        select(-to) %>%
        as.matrix() %>%
        magrittr::set_rownames(lr_network_top_df$to)

    dist_receptors <- dist(lr_network_top_matrix, method = "binary")
    hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
    order_receptors <- hclust_receptors$labels[hclust_receptors$order]

    dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

    vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]

    p_ligand_pearson <- vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands", "Ligand activity",
        color = brewer.pal(9, "Oranges")[6], legend_position = "top",
        x_axis_position = "top",
        legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)"
    )
    p_ligand_pearson <- p_ligand_pearson + theme(legend.key.width = unit(2, "cm"), axis.text.y = element_text(face = "italic", size = 12, color = "black"))
    ggsave(file = sprintf("%s.ligand_pearson.pdf", out.prefix), width = width.pearson, height = height)

    # ligand target
    p_ligand_target_network <- vis_ligand_target %>%
        make_heatmap_ggplot("Prioritized ligands", "genes in receiver cells",
            color = "#602ca3", legend_position = "top",
            x_axis_position = "top",
            legend_title = "Regulatory potential"
        ) +
        theme(axis.text.x = element_text(face = "italic", size = 12), legend.key.width = unit(2, "cm"))
    ggsave(file = sprintf("%s.ligand_target.pdf", out.prefix), width = width.target, height = height)

    # ligand combine
    figures_without_legend <- plot_grid(
        p_ligand_pearson +
            theme(legend.position = "none", axis.ticks = element_blank()) +
            theme(axis.title.x = element_text()),
        p_ligand_target_network + theme(
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank()
        ),
        align = "h", nrow = 1,
        rel_widths = rel_widths,
        rel_heights = c(
            nrow(vis_ligand_pearson),
            nrow(vis_ligand_target) + 3
        )
    )
    legends <- plot_grid(as_ggplot(get_legend(p_ligand_pearson)),
        as_ggplot(get_legend(p_ligand_target_network)),
        nrow = 2, align = "h"
    )
    pp <- plot_grid(figures_without_legend, legends, rel_heights = comb.rel.height, nrow = 2, align = "hv")
    ggsave(file = sprintf("%s.comb.pdf", out.prefix), pp, width = width.comb, height = height)
}

for (i in 1:length(myoF.clusters)) {
    gene.oi <- DElist.ls[[myoF.clusters[i]]]
    gene.bg <- DElist.ls[["background"]]

    ligand_activities <- getNichenet(
        gene.oi = gene.oi,
        gene.bg = gene.bg,
        out.prefix = sprintf("%s/%s.nichenet", figdir, myoF.clusters[i])
    )
    plotNichenet(gene.obj, gene.bg,
        out.prefix = sprintf("%s/%s.nichenet", figdir, myoF.clusters[i]),
        ligand_activities = ligand_activities, n.top = 20,
        width.target = 13.1,
        width.pearson = 1.3,
        width.comb = 19, height = 6
    )
}



# Dot plot for enriched pathway
load("./data/Enriched_pathways.RData")
dot.scale.min <- 2
dot.scale.max <- 6
scale.by <- "radius"
scale.func <- switch(
    EXPR = scale.by,
    "size" = scale_size,
    "radius" = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
)

enrichgo.top <- gsea_all %>%
    mutate(pathway = Term) %>%
    filter(pathway %in% pathway.pick)
enrichgo.top$Logqvalue <- MinMax(-enrichgo.top$Logqvalue, min = 2, max = 15)
enrichgo.top$geneCount <- MinMax(enrichgo.top$geneCount, min = 0, max = 30)

colP <- c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8")
enrichgo.top$pathway <- factor(x = enrichgo.top$pathway, levels = rev(pathway.order))
enrichgo.top$Cluster <- factor(enrichgo.top$Cluster, levels = cluster.order)
p <- ggplot(data = enrichgo.top) +
    geom_point(aes(x = Cluster, y = pathway, fill = Logqvalue, size = geneCount), pch = 21) +
    scale_fill_gradientn(colours = colP) +
    scale_size_continuous("geneCount") +
    xlab("Clusters") +
    ylab("") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = size, color = "black"),
        axis.text.y = element_text(size = size, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "left"
    )

ggsave(sprintf("%s/FigureS3A.pdf", figdir), width = 6.5, height = 7)

# PAGA connectivities
paga_prob <- read_csv("./data/paga_connectivity_all.csv")
paga_prob.m <- paga_prob %>%
    column_to_rownames(var = "X1") %>%
    as.matrix()

COL2 <- function(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200) {
    diverging <- match.arg(diverging)
    colors <- switch(diverging,
        RdBu = c(
            "#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF",
            "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"
        ),
        BrBG = c(
            "#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#FFFFFF",
            "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30"
        ),
        PiYG = c(
            "#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", "#FFFFFF",
            "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419"
        ),
        PRGn = c(
            "#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#FFFFFF",
            "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B"
        ),
        PuOr = c(
            "#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#FFFFFF",
            "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B"
        ),
        RdYlBu = c(
            "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFFF",
            "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"
        )
    )
    return(rev(colorRampPalette(colors)(n)))
}


pdf(sprintf("%s/FigureS3C.pdf", figdir), width = 8, height = 8)
g <- corrplot::corrplot(paga_prob.m[cluster.order, cluster.order],
    type = "lower", diag = FALSE, col = COL2("RdBu", 200),
    method = "pie", cl.lim = c(0, 1)
)
dev.off()

paga_prob <- read_csv("./data/paga_connectivity_progenitor_paths.csv")
paga_prob_cluster <- paga_prob %>%
    filter(X1 != "c04") %>%
    select("X1", "c04")
g <- ggbarplot(paga_prob_cluster, "X1", "c04",
    fill = "X1", palette = cluster.color, color = "white", sort.val = c("asc"), sort.by.groups = FALSE,
    x.text.angle = 90, legend = "right", width = 0.9,
    rotate = FALSE, xlab = "", ylab = "PAGA connectivities", ggtheme = theme_classic(), title = ""
) +
    ylim(0, max(paga_c04$c04)) +
    scale_y_continuous(breaks = seq(0, max(paga_c04$c04), 0.2)) +
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black"))
ggsave(sprintf("%s/FigureS3F.pdf", figdir), g, width = 3, height = 3)


# Pseudotime
dpt_path <- read_csv("./data/pseudotime_progenitor_paths.csv")
g <- ggboxplot(dpt_path,
    x = "cluster", y = "dpt_pseudotime", color = "cluster", palette = cluster.color,
    add = "jitter",
    order = rev(c("c05", "c03", "c16", "c06", "c19", "c04")),
    add.params = list(size = 0.03, alpha = 0.1), outlier.shape = NA
) + coord_flip() +
    theme(legend.position = "None") +
    labs(x = "", y = "Pseudotime")
ggsave(sprintf("%s/FigureS3E.pdf", figdir), g, width = 5, height = 4)


# Barplot of TF effect size
load("./data/TF_effectsize.RData")

dat.fig.list <- llply(gene.to.plot, function(gene) {
    vline.x <- dat.plot[geneID == gene & cluster == "c04", min(as.integer(tissue))]
    p <- ggplot(dat.plot, aes_string("tissue", "EffectSize")) +
        geom_col(aes(fill = adj.P.Val), color = NA, position = "dodge2", width = 0.8) +
        scale_fill_distiller(
            palette = "Purples", breaks = c(0, 0.002, 0.004, 0.006, 0.008, 0.01),
            limits = c(0, 0.01), na.value = "lightgray",
            guide = guide_colorbar(
                label.theme = element_text(angle = 45, hjust = 1),
                direction = "horizontal"
            )
        ) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        geom_hline(yintercept = 0.15, linetype = "dashed", alpha = 0.8, color = "lightgray") +
        theme_classic() +
        labs(x = "", y = "Effect Size", title = a.gene) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(
            axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 10)
        )
    if (!is.null(vline.x) && is.finite(vline.x)) {
        p <- p + geom_vline(xintercept = vline.x - 0.5, linetype = "dashed", color = "red", alpha = 0.8)
    }
    return(p)
}, .parallel = T)
names(dat.fig.list) <- gene.to.plot

p <- plot_grid(
    plot_grid(
        plotlist = llply(
            dat.fig.list,
            function(x) {
                x + theme(legend.position = "none")
            }
        ),
        align = "hv", ncol = 4
    ),
    get_legend(dat.fig.list[[1]]),
    ncol = 1,
    rel_heights = c(3, 1)
)
ggsave(sprintf("%s/FigureS3I.pdf", figdir), width = 8, height = 5)
