library(scProportionTest)
library(Seurat)

Idents(seurat) <- seurat$condition
prop_test <- sc_utils(seurat)
cluster <- "ann"

prop_test <- permutation_test(prop_test,
                              cluster_identity = cluster,
                              sample_1 = "RR",
                              sample_2 = "SP",
                              sample_identity = "condition", n_permutations = 1000)

