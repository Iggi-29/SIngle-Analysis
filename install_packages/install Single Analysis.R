##Install all R packages required by Single Analysis
#remote packages
install.packages("devtools")
install.packages("remote")

#shiny
install.packages("shiny")
install.packages("shinythemes")
install.packages("shinydashboard")

install.packages("plyr")

#ggplot
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("RColorBrewer")

#data manipulation
install.packages("dplyr")
install.packages("stats")

#single cell data manipullation
install.packages("Seurat")
install.packages("Matrix")
remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github("dynverse/dyno")
devtools::install_github("sqjin/CellChat")

#bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("slingshot")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")