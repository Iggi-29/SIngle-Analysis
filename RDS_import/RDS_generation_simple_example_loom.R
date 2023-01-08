###Library and WD
setwd("/home/vant/Escritorio/shinyng/data1/")
#Library
library(Seurat)
library(Matrix)
library(loomR)

##Creation of SeuratObject
mat <- Connect(filename = "./s_fca_biohub_head_10x.loom", mode = "r")

##SeuratObject
seuobj <- CreateSeuratObject(counts = mat)

##Save and read rds file
saveRDS(object = seuobj, file =  "/home/vant/Escritorio/shinyng/app/scRNA.rds")

##Read RDS data
seuobj <- readRDS(file = "/home/vant/Escritorio/shinyng/app/scRNA.rds")
Assays(seuobj)
dim(seuobj)
