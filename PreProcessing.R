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


# -------------------
sc_1743 <- Read10X(data.dir = "/mnt/calabrese/francesco/scRNAseq/analysis/1743_1_RF_RR/filtered_feature_bc_matrix")
sc_1743 <- CreateSeuratObject(counts = sc_1743, project = '1743', min.cells = 3, min.features = 200)

sc_1744 <- Read10X(data.dir = "/mnt/calabrese/francesco//scRNAseq/analysis/1744_2_GD_RR/filtered_feature_bc_matrix")
sc_1744 <- CreateSeuratObject(counts = sc_1744, project = '1744', min.cells = 3, min.features = 200)

sc_2344 <- Read10X(data.dir = "/mnt/calabrese/francesco/scRNAseq/analysis/2344_1_BM_SP/filtered_feature_bc_matrix")
sc_2344 <- CreateSeuratObject(counts = sc_2344, project = '2344', min.cells = 3, min.features = 200)

sc_2345 <- Read10X(data.dir = "/mnt/calabrese/francesco/scRNAseq/analysis/2345_2_PMG_SP/filtered_feature_bc_matrix")
sc_2345 <- CreateSeuratObject(counts = sc_2345, project = '2345', min.cells = 3, min.features = 200)

sc_2355 <- Read10X(data.dir = "/mnt/calabrese/francesco/scRNAseq/analysis/2355_3_PF_RR/filtered_feature_bc_matrix")
sc_2355 <- CreateSeuratObject(counts = sc_2355, project = '2355', min.cells = 3, min.features = 200)

sc_2804 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2804/outs/filtered_feature_bc_matrix/")
sc_2804 <- CreateSeuratObject(counts = sc_2804, project = "2804", min.cells = 3, min.features = 200)

sc_2803 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2803/outs/filtered_feature_bc_matrix/")
sc_2803 <- CreateSeuratObject(counts = sc_2803, project = "2803", min.cells = 3, min.features = 200)

sc_2777 <-  Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2777/outs/filtered_feature_bc_matrix/")
sc_2777 <- CreateSeuratObject(counts = sc_2777, project = "2777", min.cells = 3, min.features = 200)

sc_2775 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2775/outs/filtered_feature_bc_matrix/")
sc_2775 <- CreateSeuratObject(counts = sc_2775, project = "2775", min.cells = 3, min.features = 200)

sc_2774 <-  Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2774/outs/filtered_feature_bc_matrix/")
sc_2774 <- CreateSeuratObject(counts = sc_2774, project = "2774", min.cells = 3, min.features = 200)

sc_2772 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2772/outs/filtered_feature_bc_matrix/")
sc_2772 <- CreateSeuratObject(counts = sc_2772, project = "2772", min.cells = 3, min.features = 200)

sc_2652 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2652/outs/filtered_feature_bc_matrix/")
sc_2652 <- CreateSeuratObject(counts = sc_2652, project = "2652", min.cells = 3, min.features = 200)

sc_2651 <- Read10X(data.dir = "/mnt/calabrese/cellranger_counts/240508_A01083_0325_AH555NDRX5/2651/outs/filtered_feature_bc_matrix/")
sc_2651 <- CreateSeuratObject(counts = sc_2651, project = "2651", min.cells = 3, min.features = 200)

sc_2653 <- Read10X(data.dir = "~/Desktop/mnt2/bcl_fastq/RUN374NS/7_EHZ/outs/filtered_feature_bc_matrix/")
sc_2653 <- CreateSeuratObject(counts = sc_2653, project = "2653", min.cells = 3, min.features = 200)

sc_2654 <- Read10X(data.dir = "~/Desktop/mnt2/bcl_fastq/RUN374NS/8_SC/outs/filtered_feature_bc_matrix/")
sc_2654 <- CreateSeuratObject(counts = sc_2654, project = "2654", min.cells = 3, min.features = 200)

sc_2773 <- Read10X(data.dir = "~/Desktop/mnt2/bcl_fastq/RUN374NS/10_VG/outs/filtered_feature_bc_matrix/")
sc_2773 <- CreateSeuratObject(counts = sc_2773, project = "2773", min.cells = 3, min.features = 200)

sc_2776 <- Read10X(data.dir = "~/Desktop/mnt2/bcl_fastq/RUN374NS/13_FS_SCS/outs/filtered_feature_bc_matrix/")
sc_2776 <- CreateSeuratObject(counts = sc_2776, project = "2776", min.cells = 3, min.features = 200)

sc_2805 <- Read10X(data.dir = "~/Desktop/mnt2/bcl_fastq/RUN374NS/17_EI/outs/filtered_feature_bc_matrix/")
sc_2805 <- CreateSeuratObject(counts = sc_2805, project = "2805", min.cells = 3, min.features = 200)


sc_1743@meta.data$sample <- "RR_1"
sc_1743@meta.data$batch <- "1"
sc_1743@meta.data$condition <- "RR"
sc_1743@meta.data$treatment <- "mavenclad"
sc_1743@meta.data$age <- "29"
sc_1743@meta.data$sex <- "M"
sc_1743@meta.data$edss <- 1.5

sc_1744@meta.data$sample <- "RR_2"
sc_1744@meta.data$batch <- "1"
sc_1744@meta.data$condition <- "RR"
sc_1744@meta.data$treatment <- "ofatutumab"
sc_1744@meta.data$age <- "21"
sc_1744@meta.data$sex <- "F"
sc_1744@meta.data$edss <- 3.5

sc_2355@meta.data$sample <- "RR_3"
sc_2355@meta.data$batch <- "2"
sc_2355@meta.data$condition <- "RR"
sc_2355@meta.data$treatment <- "tec"
sc_2355@meta.data$age <- "50"
sc_2355@meta.data$sex <- "F"
sc_2355@meta.data$edss <- 3

sc_2344@meta.data$sample <- "SP_1"
sc_2344@meta.data$batch <- "2"
sc_2344@meta.data$condition <- "SP"
sc_2344@meta.data$treatment <- "sipo"
sc_2344@meta.data$age <- "56"
sc_2344@meta.data$sex <- "F"
sc_2344@meta.data$edss <- 4.5

sc_2345@meta.data$sample <- "SP_2"
sc_2345@meta.data$batch <- "2"
sc_2345@meta.data$condition <- "SP"
sc_2345@meta.data$treatment <- "sipo"
sc_2345@meta.data$age <- "56"
sc_2345@meta.data$sex <- "F"
sc_2345@meta.data$edss <- 3.5

sc_2652@meta.data$sample <- "RR_4"
sc_2652@meta.data$batch <- "3"
sc_2652@meta.data$condition <- "RR"
sc_2652@meta.data$treatment <- "tec"
sc_2652@meta.data$age <- "35"
sc_2652@meta.data$sex <- "F"
sc_2652@meta.data$edss <- 0

sc_2651@meta.data$sample <- "SP_3"
sc_2651@meta.data$batch <- "3"
sc_2651@meta.data$condition <- "SP"
sc_2651@meta.data$treatment <- "sipo"
sc_2651@meta.data$age <- "39"
sc_2651@meta.data$sex <- "F"
sc_2651@meta.data$edss <- 4.5

sc_2772@meta.data$sample <- "SP_4"
sc_2772@meta.data$batch <- "3"
sc_2772@meta.data$condition <- "SP"
sc_2772@meta.data$treatment <- "sipo"
sc_2772@meta.data$age <- "54"
sc_2772@meta.data$sex <- "F"
sc_2772@meta.data$edss <- 4.5

sc_2774@meta.data$sample <- "RR_5"
sc_2774@meta.data$batch <- "3"
sc_2774@meta.data$condition <- "RR"
sc_2774@meta.data$treatment <- "sipo"
sc_2774@meta.data$age <- "52"
sc_2774@meta.data$sex <- "F"
sc_2774@meta.data$edss <- 1

sc_2775@meta.data$sample <- "CTR_1"
sc_2775@meta.data$batch <- "3"
sc_2775@meta.data$condition <- "CTR"
sc_2775@meta.data$treatment <- "eculizumab"
sc_2775@meta.data$age <- "43"
sc_2775@meta.data$sex <- "F"
sc_2775@meta.data$edss <- 3

sc_2777@meta.data$sample <- "SP_5"
sc_2777@meta.data$batch <- "3"
sc_2777@meta.data$condition <- "SP"
sc_2777@meta.data$treatment <- "sipo"
sc_2777@meta.data$age <- "50"
sc_2777@meta.data$sex <- "F"
sc_2777@meta.data$edss <- 6

sc_2803@meta.data$sample <- "SP_6"
sc_2803@meta.data$batch <- "3"
sc_2803@meta.data$condition <- "SP"
sc_2803@meta.data$treatment <- "sipo"
sc_2803@meta.data$age <- "55"
sc_2803@meta.data$sex <- "F"
sc_2803@meta.data$edss <- 2

sc_2804@meta.data$sample <- "SP_7"
sc_2804@meta.data$batch <- "3"
sc_2804@meta.data$condition <- "SP"
sc_2804@meta.data$treatment <- "sipo"
sc_2804@meta.data$age <- "54"
sc_2804@meta.data$sex <- "F"
sc_2804@meta.data$edss <- 4.5

sc_2653@meta.data$sample <- "RR_6"
sc_2653@meta.data$batch <- "4"
sc_2653@meta.data$condition <- "RR"
sc_2653@meta.data$treatment <- "oza"
sc_2653@meta.data$age <- "30"
sc_2653@meta.data$sex <- "F"
sc_2653@meta.data$edss <- 3.5 

sc_2654@meta.data$sample <- "RR_7"
sc_2654@meta.data$batch <- "4"
sc_2654@meta.data$condition <- "RR"
sc_2654@meta.data$treatment <- "oza"
sc_2654@meta.data$age <- "46"
sc_2654@meta.data$sex <- "F"
sc_2654@meta.data$edss <- 2.5 

sc_2773@meta.data$sample <- "SP_8"
sc_2773@meta.data$batch <- "4"
sc_2773@meta.data$condition <- "SP"
sc_2773@meta.data$treatment <- "oza"
sc_2773@meta.data$age <- "48"
sc_2773@meta.data$sex <- "F"
sc_2773@meta.data$edss <- 3 

sc_2776@meta.data$sample <- "RR_8"
sc_2776@meta.data$batch <- "4"
sc_2776@meta.data$condition <- "RR"
sc_2776@meta.data$treatment <- "oza"
sc_2776@meta.data$age <- "37"
sc_2776@meta.data$sex <- "F"
sc_2776@meta.data$edss <- 1 

sc_2805@meta.data$sample <- "RR_9"
sc_2805@meta.data$batch <- "4"
sc_2805@meta.data$condition <- "RR"
sc_2805@meta.data$treatment <- "oza"
sc_2805@meta.data$age <- "29"
sc_2805@meta.data$sex <- "F"
sc_2805@meta.data$edss <- 1


y <- c(sc_1744,
       sc_2355,
       sc_2344,
       sc_2345,
       sc_2652,
       sc_2651,
       sc_2772,
       sc_2774,
       sc_2775,
       sc_2777,
       sc_2803,
       sc_2804,
       sc_2805,
       sc_2776,
       sc_2654,
       sc_2653,
       sc_2773)

seurat <- merge(x=sc_1743,
                y=y)

rm(sc_1743, sc_1744, sc_2355, sc_2344, sc_2345, sc_2652, sc_2651, sc_2772, sc_2774, sc_2775, sc_2777, sc_2803, sc_2804)
gc()

seurat <- JoinLayers(seurat)

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

set.seed(100)

dbl <- scDblFinder(filtered_singleCell, BPPARAM=BiocParallel::SerialParam(), samples = "sample")

seurat[["dblFinder"]] <- dbl$scDblFinder.class
seurat[["dblScore"]] <- dbl$scDblFinder.score

rm(filtered_singleCell, dbl)
gc()

seurat <- subset(seurat, subset=dblFinder=="singlet")

seurat <- NormalizeData(seurat,
                        normalization.method = "LogNormalize")

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
seurat <- ScaleData(seurat, vars.to.regress = "percent.MT")
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims=1:15)

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



