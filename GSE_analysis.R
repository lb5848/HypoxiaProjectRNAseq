rm(list = ls())

install_packages <- FALSE

if(install_packages){
  BiocManager::install("DOSE")
  BiocManager::install("clusterProfiler")
  BiocManager::install("enrichplot")
  BiocManager::install("ggupset")
  BiocManager::install("gmt")
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
library(enrichplot)
library(org.Hs.eg.db)
library(ggupset)
library(gmt)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

ts <- timestamp(stamp = Sys.Date(), prefix = "", suffix = "")
workingDirectory <- paste(PrimaryDirectory, paste(ts, "workingDirectory", sep = "_"), sep = "/")
dir.create(workingDirectory)

#create ranks
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
ranks <- sort(ranks, decreasing = TRUE)

rnk_up <- (ranks > 0)
rnk_up <- names(rnk_up)
rnk_down <- (ranks < 0)
rnk_down <- names(rnk_down)

# compareCluster
rnk_up_entrez <- bitr(rnk_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
rnk_down_entrez <- bitr(rnk_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


ego <- enrichGO(names(ranks), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
goplot(ego)
barplot(ego, showCategory = 30)
dotplot(ego, showCategory = 30)

simple_ego <- simplify(ego)
cnetplot(simple_ego, foldChange = ranks)
cnetplot(simple_ego, foldChange = ranks, circular = TRUE, colorEdge = TRUE)

upsetplot(ego)

emapplot(simple_ego)



# CompClusEnrichedPath_up <- compareCluster(geneCluster = rnk_up_entrez$ENTREZID, fun = "enrichKEGG", pvalueCutoff = 0.05)
# CompClusEnrichedPath_down <- compareCluster(geneCluster = rnk_down_entrez$ENTREZID, fun = "enrichKEGG", pvalueCutoff = 0.05)


# pathways
all.pathways <- readRDS("allPathways.rds")

data("geneList")

fgseaRes <- fgsea(all.pathways, ranks, maxSize = 500, minSize = 15, eps = 0)

plotEnrichment(all.pathways[["HALLMARK_HYPOXIA"]], ranks) + labs(title = "HALLMARK_HYPOXIA")
head(fgseaRes)


writeGmtPathways(all.pathways, gmt.file = "allPathways.gmt")

all <- read.gmt(file.choose())

gseaALL <- GSEA(ranks, exponent = 1, minGSSize = 10, maxGSSize = 500, eps = 0, pvalueCutoff = 0.05,
                TERM2GENE = all, by = "fgsea")

emapplot(gseaALL)
dotplot(gseaALL, showCategory = 30)
dotplot(sig)

gseaplot2(gseaALL, geneSetID = "Texprog2VSTterm_UP")

simple_gseaALL <- simplify(gseaALL)
cnetplot(simple_gseaALL, foldChange = ranks)
setwd(workingDirectory)
