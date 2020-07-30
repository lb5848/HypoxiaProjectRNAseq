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

workingDirectory <- paste(PrimaryDirectory, "workingDirectory", sep = "/")
dir.create(workingDirectory)

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

IL2vsNC_UP <- read_excel(xls.pnas, sheet = 2, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2vsNC_UP <- IL2vsNC_UP$Symbol
IL2vsNC_UP <- convertMouseGeneList(IL2vsNC_UP)

IL21vsNC_UP <- read_excel(xls.pnas, sheet = 3, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL21vsNC_UP <- IL21vsNC_UP$Symbol
IL21vsNC_UP <- convertMouseGeneList(IL21vsNC_UP)

IL2LDHivsIL2_UP <- read_excel(xls.pnas, sheet = 4, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2LDHivsIL2_UP <- IL2LDHivsIL2_UP$Symbol
IL2LDHivsIL2_UP <- convertMouseGeneList(IL2LDHivsIL2_UP)
  
IL21LDHivsIL21_UP <- read_excel(xls.pnas, sheet = 5, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL21LDHivsIL21_UP <- IL21LDHivsIL21_UP$Symbol
IL21LDHivsIL21_UP <- convertMouseGeneList(IL21LDHivsIL21_UP)

IL2vsIL21_UP <- read_excel(xls.pnas, sheet = 6, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2vsIL21_UP <- IL2vsIL21_UP$Symbol
IL2vsIL21_UP <- convertMouseGeneList(IL2vsIL21_UP)

IL2LDHivsIL21LDHi_UP <- read_excel(xls.pnas, sheet = 7, col_names = TRUE) %>%
  dplyr::select(Symbol, logFC, FDR) %>%
  dplyr::filter(logFC > 1) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::select(Symbol)
IL2LDHivsIL21LDHi_UP <- IL2LDHivsIL21LDHi_UP$Symbol
IL2LDHivsIL21LDHi_UP <- convertMouseGeneList(IL2LDHivsIL21LDHi_UP)

customPathways <- pathway.HALLMARK[c(1:18)]
var_names <- ls()
var_names <- var_names[grep("_UP", var_names)]
names(customPathways) <- var_names
customPathways$IL21LDHivsIL21_UP <- IL21LDHivsIL21_UP
customPathways$IL21vsNC_UP <- IL21vsNC_UP
customPathways$IL2LDHivsIL2_UP <- IL2LDHivsIL2_UP
customPathways$IL2LDHivsIL21LDHi_UP <- IL2LDHivsIL21LDHi_UP
customPathways$IL2vsIL21_UP <- IL2vsIL21_UP
customPathways$IL2vsNC_UP <- IL2vsNC_UP
customPathways$Texprog1VSTexprog2_UP <- Texprog1VSTexprog2_UP
customPathways$Texprog1VSTint_UP <- Texprog1VSTint_UP
customPathways$Texprog1VSTterm_UP <- Texprog1VSTterm_UP
customPathways$Texprog2VSTexprog1_UP <- Texprog2VSTexprog1_UP
customPathways$Texprog2VSTint_UP <- Texprog2VSTint_UP
customPathways$Texprog2VSTterm_UP <- Texprog2VSTterm_UP
customPathways$TintVSTexprog1_UP <- TintVSTexprog1_UP
customPathways$TintVSTexprog2_UP <- TintVSTexprog2_UP
customPathways$TintVSTterm_UP <- TintVSTterm_UP
customPathways$TtermVSTexprog1_UP <- TtermVSTexprog1_UP
customPathways$TtermVSTexprog2_UP <- TtermVSTexprog2_UP
customPathways$TtermVSTint_UP <- TtermVSTint_UP

all.pathways <- c(pathway.HALLMARK, customPathways)

setwd(workingDirectory)
save(all.pathways, file = "allPathways.rds")

filename <- "BM_Norm-BM_Hyp-SIG.txt"
dataDir <- "DESeqDay3"
path <- paste(PrimaryDirectory, dataDir, filename, sep = "/")
data <- read.delim(path)
rnk <-  data %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(res.rownames, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(res.rownames) %>%
  summarise(stat = mean(stat))
rev_rank <- rnk
rev_rank$stat <- rev_rank$stat * (-1)
ranks <- deframe(rev_rank)

save(ranks, file = "ranks.rds")

fgseaRes <- fgsea(all.pathways, ranks, maxSize = 500, minSize = 15, eps = 0)
plotEnrichment(all.pathways[["HALLMARK_HYPOXIA"]], ranks) + labs(title = "HALLMARK_HYPOXIA")
head(fgseaRes)
fgseaResSIG <- fgseaRes %>% dplyr::filter(padj < 0.05)
UP_sig <- fgseaResSIG %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(NES >1)
UP_sig
DOWN_sig <- fgseaResSIG %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(NES < -1)
DOWN_sig