
library(Seurat)

pseudoBulk <- function(clusters, directory = "./", seurat = seurat) {
    
    pseudo_bulk <- AggregateExpression(seurat,
                                       assays = "RNA",
                                       return.seurat = TRUE,
                                       group.by = c("condition",
                                                    "sample",
                                                    clusters),
                                       normalization.method = "LogNormalize")
    pseudo_bulk$celltype.condition <- paste(pseudo_bulk$ann,
                                            pseudo_bulk$condition,
                                            sep = "_")
    Idents(pseudo_bulk) <- pseudo_bulk$celltype.condition
    types <- unique(Idents(pseudo_bulk))
    types <- gsub("_SP", "", types)
    types <- gsub("_RR", "", types)
    types <- unique(types)
    
    results <- list()
    
    for (type in types) {
        tryCatch({
            result <- FindMarkers(object = pseudo_bulk,
                                  ident.1 = paste(type, "_SP", sep = ""),
                                  ident.2 = paste(type, "_RR", sep = ""),
                                  test.use = "DESeq2", # wilcox, wilcox_limma, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2
                                  min.cells.group = 5)
            
            name <- paste(type, sep="_")
            results[[name]] <- result
            
            file_name <- paste(name, ".csv", sep="_")
            pathdest <- file.path(directory, file_name)
            result <- cbind(gene = rownames(result), result)
            write.csv(result, file=pathdest, row.names = F)
        }, error = function(e) {
            cat("error")
        })
    }
    
    return(results)
}

degs <- pseudoBulk(clusters = "ann",
                   directory = "/.",
                   seurat = seurat)