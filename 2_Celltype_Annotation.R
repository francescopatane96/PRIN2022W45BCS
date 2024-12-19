
library(Seurat)
library(Azimuth)
library(AnnotationHub)
library(purrr)
library(dplyr)
library(tibble)

meta <- seurat@meta.data
meta <- meta[meta$condition != "CTR",]

celltokeep <- rownames(meta)
seurat <- subset(seurat, cells=celltokeep)

seurat <- RunAzimuth(seurat, reference = "pbmcref")

Idents(seurat) <- "predicted.celltype.l2"

cells_to_remove <- WhichCells(seurat, ident = c("Platelet", "ILC", "MAIT", "HSPC"))
seurat <- subset(seurat, cells = setdiff(Cells(seurat), cells_to_remove))

seurat$predicted.celltype.l2 <- ifelse(
    seurat$predicted.celltype.l2 %in% c("CD16 Mono", "CD14 Mono"),
    "Mono",
    seurat$predicted.celltype.l2
)

Idents(seurat) <- "predicted.celltype.l2"

cluster_mono <- subset(seurat, idents = "Mono")

cluster_mono <- FindNeighbors(cluster_mono, dims = 1:15)
cluster_mono <- FindClusters(cluster_mono, resolution = 0.2)

Idents(cluster_mono) <- "RNA_snn_res.0.2"

new_names <- c("Mono1", "Mono2", "Mono3", "Mono4")
names(new_names) <- levels(cluster_mono)

cluster_mono <- RenameIdents(object = cluster_mono, new_names)

Idents(seurat) <- "predicted.celltype.l2"
cells_to_remove <- WhichCells(seurat, idents = "Mono")

# Crea un nuovo oggetto Seurat escludendo le cellule identificate
seurat <- subset(seurat, cells = setdiff(Cells(seurat), cells_to_remove))

Idents(seurat) <- "predicted.celltype.l2"

cells_to_remove <- WhichCells(seurat, ident = c("CD8 Proliferating", "ASDC", "CD4 Proliferating"))
seurat <- subset(seurat, cells = setdiff(Cells(seurat), cells_to_remove))

seurat <- merge(seurat, y = cluster_mono)
print(seurat)
table(Idents(seurat))
table(Idents(seurat), seurat$condition)
seurat$ann <- Idents(seurat)

Idents(seurat) <- "ann"

seurat <- JoinLayers(seurat)
seurat <- RunPCA(seurat, npcs = 15, verbose = FALSE)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:15)



# Crea una lista unica dei tipi cellulari presenti in seurat$ann
unique_celltypes <- unique(seurat$ann)

# Assegna un numero a ogni tipo cellulare
celltype_numbers <- set_names(seq_along(unique_celltypes), unique_celltypes)

# Usa map per aggiungere i numeri ai metadati del Seurat object
seurat$ann_number <- map_chr(seurat$ann, ~ celltype_numbers[.x])

# Controlla i risultati
head(seurat@meta.data)

# Imposta la nuova colonna come identità per l'oggetto Seurat
Idents(seurat) <- "ann_number"

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
    dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)


Idents(seurat) <- "ann_number"                                 
get_conserved <- function(cluster){
    FindConservedMarkers(seurat,
                         ident.1 = cluster,
                         grouping.var = "condition",
                         only.pos = TRUE,
                         logfc.threshold = 0.25, assay = "RNA") %>%
        rownames_to_column(var = "gene") %>%
        left_join(y = unique(annotations[, c("gene_name", "description")]),
                  by = c("gene" = "gene_name")) %>%
        cbind(cluster_id = cluster, .)
}

conserved_markers <- map_dfr(c(1:24), get_conserved)

top10 <- conserved_markers %>% 
    mutate(avg_fc = (RR_avg_log2FC + SP_avg_log2FC) /2) %>% 
    group_by(cluster_id) %>% 
    top_n(n = 10, 
          wt = avg_fc)

clusters <- unique(top10$cluster_id)
for (i in clusters) {
    cluster_name <- paste0("cluster_", i)
    assign(cluster_name, top10$gene[top10$cluster_id==i])
}


