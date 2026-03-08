

# Transcriptomic Analysis of ASD-associated PWS Neurons
# Dataset: GSE178687

# 1. Load Required Libraries


library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)



# 2. Load Data

# Read RNA-seq count matrix
counts <- read.delim("data/GSE178687_scale_counts.csv", row.names = 1)

# Inspect dimensions and preview data
dim(counts)
head(counts[,1:5])


# 3. Construct Sample Metadata


metadata <- data.frame(
  row.names = colnames(counts),
  group = c(
    rep("Control", 4),
    rep("PWS_del", 4),
    rep("PWS_UPD_noASD", 4),
    rep("PWS_UPD_ASD", 4)
  )
)


# 4. Differential Expression Analysis


dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ group
)

dds <- DESeq(dds)

# Example contrast: PWS deletion vs Control
res_del_vs_control <- results(dds,
                              contrast = c("group","PWS_del","Control"))

sum(res_del_vs_control$padj < 0.05, na.rm = TRUE)



# 5. Visualization

## Volcano Plot

EnhancedVolcano(
  res_del_vs_control,
  lab = rownames(res_del_vs_control),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Volcano Plot: PWS_del vs Control"
)


## PCA Plot

vsd <- vst(dds)
plotPCA(vsd, intgroup = "group")

# 6. Mitochondrial Gene Exploration

res_df <- as.data.frame(res_del_vs_control)

mito_genes <- res_df[grep("^MT-", rownames(res_df)), ]

head(mito_genes)

# 7. ASD-Specific Differential Expression

res_asd_vs_noasd <- results(
  dds,
  contrast = c("group","PWS_UPD_ASD","PWS_UPD_noASD")
)

sum(res_asd_vs_noasd$padj < 0.05, na.rm = TRUE)


# 8. Prepare Ranked Gene List for GSEA


res_df <- as.data.frame(res_asd_vs_noasd)

res_df <- res_df[!is.na(res_df$log2FoldChange), ]

gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)

gene_list <- sort(gene_list, decreasing = TRUE)

# 9. Convert Gene Symbols to Entrez IDs

gene_conversion <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_list_filtered <- gene_list[gene_conversion$SYMBOL]
names(gene_list_filtered) <- gene_conversion$ENTREZID

gene_list_filtered <- sort(gene_list_filtered, decreasing = TRUE)

length(gene_list_filtered)

# 10. GO Biological Process Enrichment

gsea_go <- gseGO(
  geneList     = gene_list_filtered,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pvalueCutoff = 0.1,
  verbose      = FALSE
)

nrow(gsea_go@result)

dotplot(gsea_go, showCategory = 15)

# 11. KEGG Pathway Enrichment
gsea_kegg_asd <- gseKEGG(
  geneList = gene_list_filtered,
  organism = "hsa",
  pvalueCutoff = 0.1,
  verbose = FALSE
)

nrow(gsea_kegg_asd@result)

dotplot(gsea_kegg_asd, showCategory = 15)

# 12. Identify Significant KEGG Pathways


gsea_kegg_asd@result[,c("ID","Description","NES","p.adjust")]


# 13. Extract Core Enrichment Genes

core_string <- gsea_kegg_asd@result$core_enrichment[1]

core_entrez <- unlist(strsplit(core_string,"/"))

length(core_entrez)


core_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = core_entrez,
  columns = c("SYMBOL"),
  keytype = "ENTREZID"
)

core_symbols

# 14. Enrichment Visualizations

ridgeplot(gsea_kegg_asd, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")


gseaplot2(
  gsea_kegg_asd,
  geneSetID = 1,
  title = "NF-kB Signaling Pathway Enrichment"
)
