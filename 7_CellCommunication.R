library(CellChat)
library(Seurat)
options(future.globals.maxSize = 10 * 1024^3)  # 2 GiB


Idents(seurat) <- "condition"
seurat2 <- SplitObject(seurat)
RR <- seurat2$RR
SP <- seurat2$SP

data.input.RR <- RR[["RNA"]]$data # normalized data matrix
data.input.SP <- SP[["RNA"]]$data

# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
Idents(RR) <- "ann"
Idents(SP) <- "ann"

labels.RR <- Idents(RR)
labels.SP <- Idents(SP)

samples.RR <- RR$condition
samples.SP <- SP$condition

meta.RR <- data.frame(labels = labels.RR, row.names = names(labels.RR), samples = samples.RR) #samples=samples) # create a dataframe of the cell labels
meta.SP <- data.frame(labels = labels.SP, row.names = names(labels.SP), samples = samples.SP)

cellChat.RR <- createCellChat(object = RR, group.by = "labels", meta=meta.RR)
cellChat.SP <- createCellChat(object = SP, group.by = "labels", meta=meta.SP)

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB)
cellChat.RR@DB <- CellChatDB.use # aggiungere qui ogni oggetto 
cellChat.SP@DB <- CellChatDB.use

cellchat.RR <- subsetData(cellChat.RR)
cellchat.SP <- subsetData(cellChat.SP)

future::plan("multisession", workers = 5) # do parallel

gc()

cellchat.RR <- identifyOverExpressedGenes(cellchat.RR)
gc()
cellchat.SP <- identifyOverExpressedGenes(cellchat.SP)
gc()

cellchat.RR <- identifyOverExpressedInteractions(cellchat.RR)
gc()
cellchat.SP <- identifyOverExpressedInteractions(cellchat.SP)
gc()

cellchat.RR <- computeCommunProb(cellchat.RR, type = "triMean", population.size = T)
gc()
cellchat.SP <- computeCommunProb(cellchat.SP, type = "triMean", population.size = T) 
gc()

cellchat.RR <- filterCommunication(cellchat.RR, min.cells = 10)
cellchat.SP <- filterCommunication(cellchat.SP, min.cells = 10)

cellchat.RR <- computeCommunProbPathway(cellchat.RR)
cellchat.SP <- computeCommunProbPathway(cellchat.SP)

cellchat.RR <- aggregateNet(cellchat.RR)
cellchat.SP <- aggregateNet(cellchat.SP)

