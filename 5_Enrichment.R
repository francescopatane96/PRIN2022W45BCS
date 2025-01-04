library(Seurat)
library(ggplot2)
library(msigdbr)
library(SCPA)
library(dplyr)
library(tibble)
library(circlize)
library(grid)
library(tidyr)

Idents(seurat) <- "ann"

pathway_biocarta <- msigdbr("Homo sapiens", "C2", subcategory = "CP:BIOCARTA") %>% 
    format_pathways()

pathway_hallmarks <- msigdbr("Homo sapiens", "H") %>% 
    format_pathways()

# C2: curated gene sets
# CP: Canonical pathways
pathway_kegg <- msigdbr("Homo sapiens", "C2", subcategory = "CP:KEGG") %>% 
    format_pathways()
pathway_reactome <- msigdbr("Homo sapiens", "C2", subcategory = "CP:REACTOME") %>% 
    format_pathways()
pathway_wiki <- msigdbr("Homo sapiens", "C2", subcategory = "CP:WIKIPATHWAYS") %>% 
    format_pathways()
pathway_BP <- msigdbr("Homo sapiens", "C5", subcategory = "GO:BP") %>% 
    format_pathways()
pathway_MF <- msigdbr("Homo sapiens", "C5", subcategory = "GO:MF") %>% 
    format_pathways()
pathway_CC <- msigdbr("Homo sapiens", "C5", subcategory = "GO:CC") %>% 
    format_pathways()
pathway_immune <- msigdbr("Homo sapiens", "C7", subcategory = "IMMUNESIGDB") %>% 
    format_pathways()
pos_gene_sets <- msigdbr("Homo sapiens", "C1") %>% 
    format_pathways()

# Canonical Pathways gene sets derived from the PID pathway database.
# CP: Canonical pathways
pathway_pid <- msigdbr("Homo sapiens", "C2", subcategory = "CP:PID") %>% 
    format_pathways()
# C3: regulatory target gene sets
pathway_mir <- msigdbr("Homo sapiens", "C3", subcategory = "MIR:MIRDB") %>% 
    format_pathways()
pathway_tft <- msigdbr("Homo sapiens", "C3", subcategory = "TFT:GTRD") %>% 
    format_pathways()

to_do <- list(pathway_biocarta,
              pathway_hallmarks,
              pathway_kegg,
              pathway_reactome,
              pathway_wiki,
              pathway_BP,
              pathway_MF,
              pathway_CC,
              pathway_immune,
              pos_gene_sets,
              pathway_pid,
              pathway_mir,
              pathway_tft)

names(to_do) <- c("pathway_biocarta", 
                  "pathway_hallmarks", "pathway_kegg", "pathway_reactome", 
                  "pathway_wiki", "pathway_BP", "pathway_MF", "pathway_CC", 
                  "pathway_immune", "pos_gene_sets", "pathway_pid", 
                  "pathway_mir", "pathway_tft")

cell_types <- unique(seurat$ann)

seurat_condition <- SplitObject(seurat, split.by = "condition")

results_list <- list()

for (i in cell_types) {
    for (enr_name in names(to_do)) {  
        enr <- to_do[[enr_name]]  
        RRMS <- seurat_extract(seurat_condition$RR, 
                               meta1 = "ann", value_meta1 = i)
        
        PIRA <- seurat_extract(seurat_condition$SP, 
                               meta1 = "ann", value_meta1 = i)
        print(paste("Processing pathway:", enr_name, "for cell type:", i))
        
        result <- compare_pathways(
            list(RRMS, PIRA),
            enr, 
            downsample = 1000,
            min_genes = 15,
            max_genes = 500,
            parallel = TRUE,
            cores = 24
        )
        
        if (!is.null(result)) {
            df <- result %>%
                dplyr::select(Pathway, qval, FC, adjPval)
            # %>%
            #     rename(
            #         
            #         FC = FC,
            #         pval_adj = adjPval
            #     )
            
            df$cell_type <- i
            df$pathway <- enr_name  
            
            df_name <- paste(i, enr_name, sep = "_")  
            results_list[[df_name]] <- df  
        } else {
            print(paste("No result for", i, "and", enr_name))
        }
    }
}

path_regroup <- function(results_list) {
    
    pathway_grouped <- list()
    
    pathway_types <- c("custom", "metabolic", "biocarta", "hallmarks", 
                       "kegg", "reactome", "wiki", "BP", "MF", "CC", 
                       "immune", "pos", "pid", "mir", "tft")
    
    for (name in names(results_list)) {
        # Trova la tipologia del pathway (una delle parole chiave) nel nome dell'oggetto
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

plots <- list()

for (i in names(results_grouped_by_type)) {
    
    data <- results_grouped_by_type[[i]]
    
    combined <- bind_rows(data)
    
    # #creazione embeddings per poi clusterizzare tramite similaritá semantica
    # textrpp_initialize(save_profile = TRUE)
    # embeddings <- textEmbed(
    #     texts = combined$Pathway,
    #     model = "bert-base-uncased"
    # )
    
    # filtro i dati per migliore visualizzazione
    # tengo solo gene sets con pval<0.05
    # poi visualizzo i 50 con qval piú alto
    filtered_data <- combined %>%
        group_by(Pathway) %>%
        mutate(all_qval_below_2 = all(qval <= 2)) %>%  # Verifica se tutti i qval per un pathway sono sotto 2
        dplyr::filter(!all_qval_below_2) %>%  # Rimuovi i pathway dove tutti i qval sono < 2
        ungroup()  # Rimuovi il raggruppamento
    
    wide_df <- filtered_data %>%
        dplyr::select(Pathway, qval, cell_type) %>%  
        pivot_wider(names_from = cell_type, values_from = qval)
    
    wide_df2 <- wide_df %>%
        column_to_rownames("Pathway")
    
    heatmap_data <- wide_df2 %>% 
        as.matrix()
    
    # heatmap_data[is.na(heatmap_data)] <- 0
    # heatmap_data <- heatmap_data[rowSums(heatmap_data != 0) > 0, ]
    
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
        
        # Seleziona solo i pathway relativi all'immunità
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
        geom_point(aes(color = qval), alpha = 0.7) +  # Solo il colore, senza la dimensione
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
    
    # Combina i dati in un unico dataframe
    combined <- bind_rows(data)
    
    # Scorri ogni tipo di cella unico presente nel gruppo
    for (cell_type in unique(combined$cell_type)) {
        # Filtra i dati per il tipo di cella corrente
        cell_data <- combined %>%
            dplyr::filter(cell_type == !!cell_type)
        
        # Filtro per rimuovere i pathway con qval < 2 in tutte le righe
        filtered_cell_data <- cell_data %>%
            group_by(Pathway) %>%
            mutate(all_qval_below_2 = all(qval <= 2)) %>%  # Verifica se tutti i qval per un pathway sono sotto 3
            dplyr::filter(!all_qval_below_2) %>%  # Rimuovi i pathway dove tutti i qval sono < 3
            ungroup()
        
        if (length(filtered_cell_data) == 0 || nrow(filtered_cell_data) == 0) {
            # Se è vuota, salta tutto il resto del codice
            next  # Salta l'iterazione corrente del ciclo (se presente)
        } else {
            
            if (group %in% c("custom", "BP", "immune")) {
                # Filtra i pathway in base ai termini chiave di immunità
                filtered_pathways <- filtered_cell_data %>%
                    dplyr::filter(sapply(Pathway, function(p) any(sapply(immunity_keywords, function(k) grepl(k, p, ignore.case = TRUE)))))
            } else {
                filtered_pathways <- filtered_cell_data
            }
            # Filtra i dati positivi e negativi in base a FC
            df_positive <- filtered_pathways %>% dplyr::filter(FC >= 0.5)
            df_negative <- filtered_pathways %>% dplyr::filter(FC <= -0.5)
            
            # Crea e salva il dot plot per i pathway positivi
            if (nrow(df_positive) > 0) {
                dot_plot_positive <- create_dot_plot(df_positive, paste0(group, "_", cell_type, " Enrichment Plot (p-value < 0.05, FC > 0.5)"))
                print(dot_plot_positive)
                
                # Salva il dot plot positivo nella lista
                dot_plots[[paste0(group, "_", cell_type, "_positive")]] <- dot_plot_positive
            } else {
                print(paste0("Nessun pathway con p-value aggiustato < 0.05 e FC positivo trovato per ", group, "_", cell_type))
            }
            
            # Crea e salva il dot plot per i pathway negativi
            if (nrow(df_negative) > 0) {
                dot_plot_negative <- create_dot_plot(df_negative, paste0(group, "_", cell_type, " Enrichment Plot (p-value < 0.05, FC < -0.5)"))
                print(dot_plot_negative)
                
                # Salva il dot plot negativo nella lista
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
