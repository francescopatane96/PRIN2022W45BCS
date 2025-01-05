library(hdWGCNA)
library(cowplot)
library(Seurat)
library(patchwork)
library(UCell)
enableWGCNAThreads(nThreads = 20)

theme_set(theme_cowplot())

## -------- 1. Pre-processing and Network creation ---------- ##

seurat_obj <- SetupForWGCNA(
  seurat,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "MS" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("ann", "condition"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'ann' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Mono1"),
  group.by='ann'
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Mono1' # name of the topological overlap matrix written to disk
)

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="condition"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'ann', group_name = 'Mono'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Mono1-M"
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

## 1. Differential module eigengene (DME) analysis

group1 <- seurat_obj@meta.data %>% subset(ann == 'Mono1' & condition == "SP") %>% rownames
#ref
group2 <- seurat_obj@meta.data %>% subset(ann == 'Mono1' & condition == "RR") %>% rownames

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  wgcna_name = 'MS', 
  test.use='wilcox'
)

#one-versus-all DME analysis

group.by = 'ann'

DMEs_all <- FindAllDMEs(
  seurat_obj,
  group.by = 'ann',
  wgcna_name = 'MS'
)

seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

## 2.  Module Trait Correlation
seurat_obj$condition <- as.factor(seurat_obj$condition)
seurat_obj$sex <- as.factor(seurat_obj$sex)
seurat_obj$batch <- as.factor(seurat_obj$batch)

seurat_obj$condition <- as.factor(seurat_obj$condition)
seurat_obj$batch <- as.factor(seurat_obj$batch)

# list of traits to correlate
cor_traits <- c("batch", "condition")

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cor_traits,
  group.by='ann'
)
## 3. Enrichment

