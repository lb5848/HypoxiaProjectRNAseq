rm(list = ls())

install_packages <- TRUE

if(install_packages){
  BiocManager::install("DOSE")
  BiocManager::install("clusterProfiler")
}

library(DOSE)
library(clusterProfiler)
library(tidyverse)
library(readr)
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(hypeR)
library(gprofiler2)
library(magrittr)
library(data.table)
library(biomaRt)
library(readxl)

library(org.Hs.eg.db)

