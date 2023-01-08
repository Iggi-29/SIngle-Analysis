###Library and WD
setwd("/home/vant/Escritorio/shinyng/data1/")
#Library
library(Seurat)
library(Matrix)

##Creation of SeuratObject
mat <- Read10X_h5(filename = "./filtered_feature_bc_matrix.h5",
               use.names = T,
               unique.features = T)

##SeuratObject
seuobj <- CreateSeuratObject(counts = mat)

##Save and read rds file
saveRDS(object = seuobj, file =  "/home/vant/Escritorio/shinyng/app/scRNA.rds")

##Read RDS data
seuobj <- readRDS(file = "/home/vant/Escritorio/shinyng/app/scRNA.rds")
Assays(seuobj)
dim(seuobj)
