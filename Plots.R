library(ggplot2)
library(Seurat)
library(gridExtra)
library(cowplot)

## -------- QC plots ------- ##

vln_metrics_pre_sample <- VlnPlot(seurat,
             layer = "counts",
             features=c("nUMI",
                        "nGene",
                        "percent.MT",
                        "percent.RB",
                        "percent.Largest.Gene"),
             log = T,
             alpha = 0.006,
             fill.by = "condition",
             pt.size = 0.01, ncol = 2)

vln_metrics_pre_trait <- VlnPlot(seurat,
             layer = "counts",
             features=c("nUMI",
                        "nGene",
                        "percent.MT",
                        "percent.RB",
                        "percent.Largest.Gene"),
             log = T,
             alpha = 0.005,
             fill.by = "condition",
             pt.size = 0.01,
             group.by = "condition")

histo_numcells_pre <- seurat@meta.data %>%
    ggplot(aes(x = seq_folder, fill = condition)) +
    geom_bar() +
    facet_wrap(~ batch) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold")) +
    ggtitle("NCells per batch by condition")+
    scale_y_continuous(breaks = seq(0,
                                    30000,
                                    by = 5000))

percent_LargGene_scatter <- FeatureScatter(seurat,
                                           feature1 = "percent.Largest.Gene",
                                           feature2 = "nUMI")+
    geom_vline(xintercept = 3,
               color = "red",
               linetype = "dashed",
               size = 0.5) 

# complexity -- use seurat initial matrices (before pre-processing)
qc.metrics %>%
    mutate(complexity=log10(nGene) / log10(nUMI))  -> qc.metrics

lm(log10(qc.metrics$nGene)~log10(qc.metrics$nUMI)) -> complexity.lm

qc.metrics %>%
    mutate(
        complexity_diff = log10(nGene) - ((log10(qc.metrics$nUMI)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
    ) -> qc.metrics

complexity_diff_plot <- qc.metrics %>%
    ggplot(aes(x=complexity_diff)) +
    geom_density(fill="yellow")

min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) -> complexity_scale

complexity_diff_plot2 <- qc.metrics %>%
    mutate(complexity_diff=replace(complexity_diff,complexity_diff< -0.1,-0.1)) %>%
    ggplot(aes(x=log10(nUMI), y=log10(nGene), colour=complexity_diff)) +
    geom_point(size=0.5) +
    geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
    scale_colour_gradient2(low="blue2",mid="grey",high="red2")

Complex_diff_largGene <- qc.metrics %>%
    ggplot(aes(x=complexity_diff,
               y=percent.Largest.Gene,
               col=seq_folder,
               shape=batch)) +
    geom_point(alpha=0.6)

nUMI_sample_distrib <- metadata %>%
    ggplot(aes(color=condition, x=nUMI, fill=seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell Density") + ggtitle("Number of UMIs per cell") +
    geom_vline(xintercept = 1000)

nGene_sample_distrib <- metadata %>%
    ggplot(aes(color=condition, x=nGene, fill=seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() + ggtitle("Number of genes per cell") +
    geom_vline(xintercept = 500)

nGenes_distrib_cond <- metadata %>% 
    ggplot(aes(x=seq_folder, y=log10(nGene), fill=condition)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("nGenes per sample") 

nGenes_distrib_batch <- metadata %>% 
    ggplot(aes(x=seq_folder, y=log10(nGene), fill=batch)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("nGenes per sample") 

percentLargGene_distrib_aggregate <- metadata %>%
    ggplot(aes(percent.Largest.Gene)) + 
    geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
    ggtitle("Distribution of Percentage largest gene") +
    geom_vline(xintercept = 5)

percentLargGene_distrib_sample <- metadata %>%
    ggplot(aes(x = percent.Largest.Gene)) + 
    geom_histogram(binwidth = 0.5, fill = "yellow", colour = "black") +
    ggtitle("Distribution of Percentage Largest Gene by Sample") +
    geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
    facet_wrap(~ seq_folder, scales = "free_y") + scale_y_log10() + # Crea un pannello per ogni campione
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

mt_per_cell_ditrib <- metadata %>%
    ggplot(aes(color=batch, x=percent.MT, fill=seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    geom_vline(xintercept = 10) + ggtitle("%MT per cell")

rb_per_cell <- metadata %>%
    ggplot(aes(color=batch, x=percent.RB, fill=seq_folder)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ggtitle("%RB per cell")

novelty_score_sample <- metadata %>%
    ggplot(aes(x=log10GenesPerUMI,color=condition, fill=seq_folder)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) + ggtitle("Novelty score")

nUMI_per_cell <- metadata %>%
    ggplot(aes(x=nUMI, color=batch, fill=seq_folder)) +
    geom_density(alpha = 0.1) +
    theme_classic() +
    geom_vline(xintercept = 5000) +
    geom_vline(xintercept = 20000)


PCA_plot <- DimPlot(seurat,
                    reduction = "pca",
                    group.by = c("condition","batch")) 

UMAP_plot <- DimPlot(seurat,
                     reduction = "umap",
                     group.by = c("condition", "seq_folder", "batch"))

UMAP_splitted <- DimPlot(seurat,
        reduction = "umap",
        group.by = c("condition", "batch"), split.by = "condition")
## --------------------------------------------------------------------##

## ------------------------ Cell type annotation --------------------- ##

p3 <- DimPlot(seurat,
              reduction = "ref.umap",
              group.by = c("predicted.celltype.l2"),
              label = T,
              label.size = 5,
              repel = T) + NoLegend()
p4 <- FeaturePlot(seurat,
                  features = "predicted.celltype.l2.score",
                  reduction = "ref.umap")
p3+p4

umap_mono <- DimPlot(cluster_mono,
                     reduction = "umap",
                     label = TRUE, pt.size = 0.5) +
    DimPlot(cluster_mono,
            reduction = "pca",
            label = TRUE, pt.size = 0.5)

Idents(seurat) <- "ann"
umap_aggr_filtered <- DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.5, repel = T) 

seurat2 <- subset(seurat, idents=c("23", "24", "21", "22"))
p1 <- DotPlot(seurat2, features = rev(cluster_21), cols = c("blue", "red"), dot.scale = 8, 
              split.by = "condition", assay = "RNA") + RotatedAxis()
p2 <- DotPlot(seurat2, features = rev(cluster_22), cols = c("blue", "red"), dot.scale = 8, 
              split.by = "condition", assay = "RNA") + RotatedAxis()
p3 <- DotPlot(seurat2, features = rev(cluster_23), cols = c("blue", "red"), dot.scale = 8, 
              split.by = "condition", assay = "RNA") + RotatedAxis()
p4 <- DotPlot(seurat2, features = rev(unique(cluster_24)), cols = c("blue", "red"), dot.scale = 8, 
              split.by = "condition", assay = "RNA") + RotatedAxis()
p1|p2|p3|p4

# markers for every cluster

for (i in clusters) {
    cluster_name <- paste0("cluster_", i)
    assign(cluster_name, top10$gene[top10$cluster_id == i])
}

p <- list()

# Genera un DotPlot per ogni cluster
for (i in clusters) {
    cluster_name <- paste0("cluster_", i)
    genes <- get(cluster_name)  
    
    p[[length(p) + 1]] <- DotPlot(seurat, 
                                  features = rev(genes), 
                                  cols = c("blue", "red"), 
                                  dot.scale = 8, 
                                  split.by = "condition", 
                                  assay = "RNA") + 
        RotatedAxis() + 
        ggtitle(paste("DotPlot for Cluster", i))  
}

## --------------------------- Compositional analysis --------------- ##
permutation_plot(prop_test)

## -------------------------------------------------------------------##

## -------------------------- DE ------------------------------------- ##
a <- ggplot(pca_df, aes(x = PC_1, y = PC_2, color = cell.type, shape=condition)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PCA of expr. Pseudobulk Data",
         x = "1st PC",
         y = "2nd PC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

b <- ggplot(pca_df, aes(x = PC_1, y = PC_2, color = batch, shape = condition)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PCA of Pseudobulk expr. Data",
         x = "1st PC",
         y = "2nd PC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

a + b

types <- unique(seurat$ann)
plots_a <- list()

for (i in types) {
    temp_df <- pca_df[pca_df$cell.type == i, ]
    
    a <- ggplot(temp_df, aes(x = PC_1, y = PC_2, color = batch, shape = condition)) +
        geom_point(size = 3, alpha = 0.7) +
        labs(title = paste("Cell Type:", i),
             x = "PC1",
             y = "PC2") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right")  
    
    plots_a[[i]] <- a
}

for (i in seq_along(plots_a)){
    filename <- paste0("/home/francesco/Desktop/plotsdrive/pca_pseudobulk_", i, ".png")
    
    ggsave(filename = filename,
           plot = plots_a[[i]],
           width = 8,   # Larghezza del file in pollici
           height = 6,  # Altezza del file in pollici
           dpi = 300)   
}

# DEGs

sets <- lapply(names(degs), function(cell_type) {
    degs[[cell_type]]
})
names(sets) <- names(degs) 

filtered_degs <- lapply(sets, function(x) {
    filtered <- x[x$p_val_adj <= 0.1, ]
    na.omit(filtered)
})

filtered_degs <- filtered_degs[sapply(filtered_degs, nrow) > 0]

filtered_degs <- lapply(names(filtered_degs), function(cell_type) {
    data <- filtered_degs[[cell_type]]
    data$type <- cell_type  
    return(data)
})

data <- bind_rows(filtered_degs)

data$avg_log2FC <- as.numeric(data$avg_log2FC)
data$p_val_adj <- as.numeric(data$p_val_adj)
data$absFC <- abs(data$avg_log2FC)

data <- data %>%
    mutate(logFC_sign = ifelse(avg_log2FC > 0, "positive", "negative"))

shapes <- c("positive" = 24, "negative" = 25) 

degs_plot <- ggplot(data, aes(x = type, y = rownames(data), size = absFC, color = -log10(p_val_adj), shape = logFC_sign)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(0.5, 5)) +  
    scale_color_gradient(low = "blue", high = "red", limits = c(0, max(-log10(data$p_val_adj)))) +
    scale_shape_manual(values = shapes) +  
    theme_bw() +
    labs(title = "PIRA vs no PIRA",
         x = "Cell Type",
         y = "DEG",
         size = "Average Log2 Fold Change",
         color = "adjusted Log10 Pval",
         shape = "Direction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
    guides(shape = guide_legend(override.aes = list(size = 3)),
           size = guide_legend(override.aes = list(shape = 24)))
## ------------------------------------------------------------------ ##

## -------------------------- Enrichment ----------------------------- ##
path_regroup <- function(results_list) {
    
    pathway_grouped <- list()
    
    pathway_types <- c("custom", "metabolic", "biocarta", "hallmarks", 
                       "kegg", "reactome", "wiki", "BP", "MF", "CC", 
                       "immune", "pos", "pid", "mir", "tft")
    
    for (name in names(results_list)) {
        pathway_type <- NULL
        for (type in pathway_types) {
            if (grepl(type, name)) {
                pathway_type <- type
                break
            }
        }
        
        if (!is.null(pathway_type)) {
            if (pathway_type %in% names(pathway_grouped)) {
                pathway_grouped[[pathway_type]] <- c(pathway_grouped[[pathway_type]], list(results_list[[name]]))
            } else {
                pathway_grouped[[pathway_type]] <- list(results_list[[name]])
            }
        }
    }
    
    return(pathway_grouped)
}

results_grouped_by_type <- path_regroup(results_list)

############################ HEATMAP #############################
plots <- list()

# library(text)
# library(torch)

for (i in names(results_grouped_by_type)) {
    
    data <- results_grouped_by_type[[i]]
    
    combined <- bind_rows(data)
    
    filtered_data <- combined %>%
        group_by(Pathway) %>%
        mutate(all_qval_below_2 = all(qval <= 2)) %>%  
        filter(!all_qval_below_2) %>%  
        ungroup()  
    
    wide_df <- filtered_data %>%
        dplyr::select(Pathway, qval, cell_type) %>%  
        pivot_wider(names_from = cell_type, values_from = qval)
    
    wide_df2 <- wide_df %>%
        column_to_rownames("Pathway")
    
    heatmap_data <- wide_df2 %>% 
        as.matrix()
    
    qval_colors <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 3, 6))
    
    immunity_keywords <- c("immune",
                           "immunity",
                           "lymphocyte",
                           "antigen",
                           "adaptive",
                           "cytokine",
                           "inflammation",
                           "inflammatory",
                           "interferon",
                           "b_cell",
                           "leukocyte",
                           "mhc",
                           "adhesion",
                           "migration",
                           "chemotaxis",
                           "apoptotic",
                           "interleukin",
                           "microglial",
                           "apoptosis",
                           "signaling",
                           "macrophage")
    
    if (i %in% c("custom", "BP", "immune")) {
        
        filtered_rows <- grep(paste(immunity_keywords, collapse = "|"), rownames(heatmap_data), ignore.case = TRUE)
        
        filtered_heatmap_data <- heatmap_data[filtered_rows, , drop = FALSE]
    } else {
        filtered_heatmap_data <- heatmap_data
    }
    
    plot <- ComplexHeatmap::Heatmap(
        filtered_heatmap_data, 
        name = "Qval",  
        col = qval_colors, 
        row_title = "Gene Set",  
        show_row_names = TRUE,  
        row_names_gp = gpar(fontsize = 6), 
        column_names_gp = gpar(fontsize = 10), 
        border = TRUE,  
        cluster_rows = TRUE,    
        cluster_columns = FALSE,  
        row_names_side = "right",  
        column_labels = cell_types,
        column_names_rot = 45    
    )
    
    plots[[i]] <- plot
}

plots

for (i in seq_along(plots)) {
    file_name <- paste0("/home/francesco/Desktop/enrich/plot_", i, ".png")
    png(filename = file_name, width = 3000, height = 3000, res = 330)
    print(plots[[i]])
    dev.off()
}

################################## DOTPLOT ###############################

create_dot_plot <- function(data, title) {
    
    ggplot(data, aes(x = reorder(Pathway, -FC), y = FC)) +
        geom_point(aes(color = qval), alpha = 0.7) +  
        scale_color_gradient(low = "blue", high = "red") + 
        labs(title = title, x = "Gene Set", y = "Fold Change (FC)",
             color = "Q-value") +
        theme_minimal() +
        coord_flip() +
        theme(axis.text.y = element_text(size = 7))  
    
}

if (!exists("dot_plots")) {
    dot_plots <- list()
} else {
    dot_plots <- list()  
}

immunity_keywords <- c("immune",
                       "immunity",
                       "lymphocyte",
                       "antigen",
                       "adaptive",
                       "cytokine",
                       "inflammation",
                       "inflammatory",
                       "interferon",
                       "b_cell",
                       "leukocyte",
                       "mhc",
                       "adhesion",
                       "migration",
                       "chemotaxis",
                       "apoptotic",
                       "interleukin",
                       "microglial",
                       "apoptosis",
                       "signaling",
                       "macrophage")

for (group in names(results_grouped_by_type)) {
    data <- results_grouped_by_type[[group]]
    
    combined <- bind_rows(data)
    
    for (cell_type in unique(combined$cell_type)) {
        cell_data <- combined %>%
            filter(cell_type == !!cell_type)
        
        filtered_cell_data <- cell_data %>%
            group_by(Pathway) %>%
            mutate(all_qval_below_2 = all(qval <= 2)) %>% 
            filter(!all_qval_below_2) %>%  
            ungroup()
        
        if (length(filtered_cell_data) == 0 || nrow(filtered_cell_data) == 0) {
            next  
        } else {
            
            if (group %in% c("custom", "BP", "immune")) {
                filtered_pathways <- filtered_cell_data %>%
                    filter(sapply(Pathway, function(p) any(sapply(immunity_keywords, function(k) grepl(k, p, ignore.case = TRUE)))))
            } else {
                filtered_pathways <- filtered_cell_data
            }
            df_positive <- filtered_pathways %>% filter(FC >= 0.5)
            df_negative <- filtered_pathways %>% filter(FC <= -0.5)
            
            if (nrow(df_positive) > 0) {
                dot_plot_positive <- create_dot_plot(df_positive, paste0(group, "_", cell_type, " Enrichment Plot (p-value < 0.05, FC > 0.5)"))
                print(dot_plot_positive)
                
                dot_plots[[paste0(group, "_", cell_type, "_positive")]] <- dot_plot_positive
            } else {
                print(paste0("Nessun pathway con p-value aggiustato < 0.05 e FC positivo trovato per ", group, "_", cell_type))
            }
            
            if (nrow(df_negative) > 0) {
                dot_plot_negative <- create_dot_plot(df_negative, paste0(group, "_", cell_type, " Enrichment Plot (p-value < 0.05, FC < -0.5)"))
                print(dot_plot_negative)
                dot_plots[[paste0(group, "_", cell_type, "_negative")]] <- dot_plot_negative
            } else {
                print(paste0("Nessun pathway con p-value aggiustato < 0.05 e FC negativo trovato per ", group, "_", cell_type))
            }
        }
    }
}

dot_plots

for (i in seq_along(dot_plots)) {
    file_name <- paste0("/home/francesco/Desktop/enrich/dot/plot_", i, ".png")
    png(filename = file_name, width = 3500, height = 2000, res = 330)
    print(dot_plots[[i]])
    dev.off()
}
## ------------------------------------------------------------------- ##

## ------------------------------ Network ----------------------------- ##
## ------------------------------------------------------------------ ##

## -------------------------------- cell-cell comm. -----------------------##

## ----------------------------------------------------------------- ##
