###Library and WD
setwd("/home/vant/Escritorio/def_app/scripts data import//")
#Library
library(Seurat)
library(Matrix)
library(Matrix)
library(tictoc)

##Data import
#Path imputation
matrix <- ("./matrix.mtx")
features <- ("./genes.tsv")
barcodes <- ("./barcodes.tsv")

##Creation of SeuratObject
mat <- ReadMtx(mtx = matrix,
               features = features,
               cells = barcodes, strip.suffix = F)

##SeuratObject
seuobj <- CreateSeuratObject(counts = mat)

##Save and read rds file
saveRDS(object = seuobj, file =  "/home/vant/Escritorio/shinyng/app/scRNA2.rds")

##Read RDS data
seuobj <- readRDS(file = "/home/vant/Escritorio/shinyng/app/scRNA.rds")
Assays(seuobj)
DefaultAssay(seuobj)
