##Library and WD
setwd("/home/vant/Escritorio/def_app/aa_def/data_raw/pbmc_human/")
#Library
library(Seurat)
library(Matrix)
library(Matrix)
library(RECODE)
library(tictoc)

##Data import
#Path imputation
matrix <- ("./matrix.mtx")
features <- ("./genes.tsv")
barcodes <- ("./barcodes.tsv")

#Creation of SeuratObject
mat <- ReadMtx(mtx = matrix,
               features = features,
               cells = barcodes, strip.suffix = T)

##SeuratObject
seuobj <- CreateSeuratObject(counts = mat)
data <- as.matrix(seuobj[["RNA"]]@counts)

##RECODE
tic()
data_RECODED <- RECODE(data)
toc()

#Enter denoised data into Seurat
seuobj[["RECODED"]] <- CreateAssayObject(Matrix(data_RECODED, sparse = T))
DefaultAssay(seuobj) <- "RECODED"

##Save and read rds file
saveRDS(object = seuobj, file =  "../../data/pbmc_recoded/pbmc_RECODED.rds")

##Read RDS data
seuobj <- readRDS(file = "../../data/pbmc_recoded/pbmc_RECODED.rds")
Assays(seuobj)
DefaultAssay(seuobj)
