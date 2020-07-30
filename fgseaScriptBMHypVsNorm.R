rm(list = ls())

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

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# load pathway
pathway.HALLMARK <- gmtPathways(file.choose())
head(pathway.HALLMARK)

# generate custom pathways from 
# 1- Immunity 2020, Beltra et al
# 2- PNAS, IL21 paper

# convert mouseGenes to human SYMBOL
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = TRUE)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
xls.file <- "geneSignatures.xlsx"
xls.file <- paste(PrimaryDirectory, "Beltra et al Immunity 2020 gene list", xls.file, sep = "/")
Texprog1 <- read_excel(xls.file, sheet = 1)
Texprog2 <- read_excel(xls.file, sheet = 2)
Tint <- read_excel(xls.file, sheet = 3)
Tterm <- read_excel(xls.file, sheet = 4)

Texprog1VSTexprog2_UP <- convertMouseGeneList(Texprog1$...1[Texprog1$`Texprog1 vs Texprog2` == 1])
Texprog1VSTint_UP <- convertMouseGeneList(Texprog1$...1[Texprog1$`Texprog1 vs Texint` == 1])
Texprog1VSTterm_UP <- convertMouseGeneList(Texprog1$...1[Texprog1$`Texprog1 vs Texterm` == 1])

Texprog2VSTexprog1_UP <- convertMouseGeneList(Texprog2$...1[Texprog2$`Texprog2 vs Texprog1` == 1])
Texprog2VSTint_UP <- convertMouseGeneList(Texprog2$...1[Texprog2$`Texprog2 vs Texint` == 1])
Texprog2VSTterm_UP <- convertMouseGeneList(Texprog2$...1[Texprog2$`Texprog2 vs Texterm` == 1])

TintVSTexprog1_UP <- convertMouseGeneList(Tint$...1[Tint$`Texint vs Texprog1` == 1])
TintVSTexprog2_UP <- convertMouseGeneList(Tint$...1[Tint$`Texint vs Texprog2` == 1])
TintVSTterm_UP <- convertMouseGeneList(Tint$...1[Tint$`Texint vs Texterm` == 1])

TtermVSTexprog1_UP <- convertMouseGeneList(Tterm$...1[Tterm$`Texterm vs Texprog1` == 1])
TtermVSTexprog2_UP <- convertMouseGeneList(Tterm$...1[Tterm$`Texterm vs Texprog2` == 1])
TtermVSTint_UP <- convertMouseGeneList(Tterm$...1[Tterm$`Texterm vs Texint` == 1])

xls.pnas <- "pnas.1920413117.sd01.xls"
xls.pnas <- paste(PrimaryDirectory, "IL21 PNAS paper", xls.pnas, sep = "/")
pnas.ALL <- read_excel(xls.pnas, sheet = 2)
