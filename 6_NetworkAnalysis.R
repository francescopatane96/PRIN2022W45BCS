library(hdWGCNA)
library(cowplot)
library(Seurat)
library(patchwork)
enableWGCNAThreads(nThreads = 20)

theme_set(theme_cowplot())

seurat_obj <- SetupForWGCNA(
    seurat,
    gene_select = "fraction", 
    fraction = 0.05, 
    wgcna_name = "MS-PRIN" 
)

#construct metacells

seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("ann", "condition"), 
    reduction = 'umap', 
    k = 25, 
    max_shared = 10, 
    ident.group = 'ann' 
)

seurat_obj <- NormalizeMetacells(seurat_obj)

# co-expression analysis

types <- unique(seurat_obj$ann)

# Rimuovi i tipi specifici di cellule
types <- setdiff(types, c("cDC1", "B naive", "CD4 Naive", "NK Proliferating", "Mono4"))

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = types,
  group.by='ann'
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# construct the network
seurat_obj <- ConstructNetwork(
    seurat_obj,
    tom_name = 'Mono1' # name of the topoligical overlap matrix written to disk
)

# PlotDendrogram(seurat_obj, main='Mono1 hdWGCNA Dendrogram')

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
    group.by = 'ann', group_name = 'Mono1'
)

# rename the modules
seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = "Mono1-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='UCell'
)

# Differential module eigengene (DME) analysis

group1 <- seurat_obj@meta.data %>% subset(ann == 'Mono1' & condition == "SP") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(ann == 'Mono1' & condition == "RR") %>% rownames

DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group1,
    barcodes2 = group2,
    wgcna_name = 'tutorial', 
    test.use='wilcox'
)

PlotDMEsLollipop(
    seurat_obj, 
    DMEs,
    wgcna_name='tutorial', 
    pvalue = "p_val_adj"
)


# set up an empty dataframe for the DMEs
DMEs <- data.frame()

seurat2 <- seurat$condition
seurat2 <- gsub("RR", 0, seurat2)
seurat2 <- gsub("SP", 1, seurat2)
seurat_obj@meta.data$condition <- as.numeric(seurat2)
seurat_obj@meta.data$condition <- as.list(seurat_obj@meta.data$condition)

# loop through the clusters
for(cur_cluster in unique(seurat$ann)){
    
    # identify barcodes for group1 and group2 in eadh cluster
    group1 <- seurat_obj@meta.data %>% subset(annotation == cur_cluster & condition == 0) %>% rownames
    group2 <- seurat_obj@meta.data %>% subset(annotation == cur_cluster & condition == 1) %>% rownames
    
    # run the DME test
    cur_DMEs <- FindDMEs(
        seurat_obj,
        barcodes1 = group1,
        barcodes2 = group2,
        test.use='wilcox',
        pseudocount.use=0.01, # we can also change the pseudocount with this param
        wgcna_name = 'Mono1'
    )
    
    # add the cluster info to the table
    cur_DMEs$cluster <- cur_cluster
    
    # append the table
    DMEs <- rbind(DMEs, cur_DMEs)
}

# Differential regulon analysis
