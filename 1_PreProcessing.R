library(Seurat)
library(ggplot2)
library(tidyverse)
library(parallel)
library(BiocSingular)
library(scDblFinder)
library(RCurl)
library(AnnotationHub)
library(celldex)
library(dplyr)
set.seed(100)
# -------------------

seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA)/log10(seurat$nCount_RNA)
seurat[["percent.RB"]] <- PercentageFeatureSet(seurat,pattern="^RP[LS]") 
seurat[["percent.MT"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[rownames(seurat) != "MALAT1",] -> seurat.nomalat
apply(
    LayerData(seurat.nomalat, layer="counts"),
    2,
    max
) -> seurat.nomalat$largest_count
apply(
    LayerData(seurat.nomalat, layer="counts"),
    2,
    which.max
) -> seurat.nomalat$largest_index
rownames(seurat.nomalat)[seurat.nomalat$largest_index] -> seurat.nomalat$largest_gene
100 * seurat.nomalat$largest_count / seurat.nomalat$nCount_RNA -> seurat.nomalat$percent.Largest.Gene
seurat.nomalat$largest_gene -> seurat$largest_gene
seurat.nomalat$percent.Largest.Gene -> seurat$percent.Largest.Gene

rm(seurat.nomalat)

metadata <- seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>% 
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
seurat@meta.data <- metadata

seurat <- subset(x=seurat,
                 subset = (nUMI >= 1000) &
                     (nGene >= 1000) &
                     (log10GenesPerUMI >= 0.8) &
                     (percent.MT <= 10) &
                     (percent.Largest.Gene <= 3))

counts <- GetAssayData(object = seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes,]

seurat <- CreateSeuratObject(filtered_counts, meta.data = seurat@meta.data)
filtered_singleCell <- as.SingleCellExperiment(seurat)

dbl <- scDblFinder(filtered_singleCell, BPPARAM=BiocParallel::SerialParam(), samples = "sample")
seurat[["dblFinder"]] <- dbl$scDblFinder.class
seurat[["dblScore"]] <- dbl$scDblFinder.score
rm(filtered_singleCell, dbl)
gc()

seurat <- subset(seurat, subset=dblFinder=="singlet")
seurat <- NormalizeData(seurat,
                        normalization.method = "LogNormalize")
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$sample)

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
seurat <- ScaleData(seurat, vars.to.regress = "percent.MT")
seurat <- RunPCA(seurat)

seurat <- IntegrateLayers(
  object = seurat,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

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

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
    dplyr::filter(phase == "S") %>%
    pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
    dplyr::filter(phase == "G2/M") %>%
    pull("gene_name")

seurat <- CellCycleScoring(seurat,
                           g2m.features = g2m_genes,
                           s.features = s_genes)



