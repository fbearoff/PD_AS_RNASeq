# Inputs ------------------------------------------------------------------
WD <- "/home/frank/R_projects/FB38_June2021_30-506016007"
condition1 <- "10ug_ml_ASM"
condition2 <- "vehicle"

# Load Libraries ----------------------------------------------------------
library(BiocParallel)
library(tximport)
library(readr)
library(DESeq2)
library(data.table)
library(ggplot2)
library(gplots)
library(genefilter)
library(pheatmap)
library(ggrepel)
library(viridis)
register(MulticoreParam(10))

# TxImport ----------------------------------------------------------------
setwd(WD)
samples <-
    read.table(file.path(WD, "samples.txt"), header = TRUE)
samples$mapped_percent <-
    read.table(file.path(WD, "salmon_output", "mapped_percent.txt"))[, 2]
files <- file.path(WD, "salmon_output", "quants", samples$sample.id, "quant.sf")
names(files) <- samples$sample.id
all(file.exists(files))
tx2gene <- read.csv(file.path(WD, "tx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# DESeq2 ------------------------------------------------------------------
rownames(samples) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, samples, ~ condition)
dds <- DESeq(dds, parallel = TRUE)

# PCA Plot ------------------------------------------------------------
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "RIN", "DV200"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = "pca.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=2) + geom_text_repel(aes(label=name), show.legend=FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + labs(title = "Variance Stabilizing Transformation PCA Plot")
dev.off()

# Clustered Sample Distance Plot ------------------------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$condition, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- plasma(255)
pdf(file = "clustered_distance.pdf")
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
dev.off()

# Clustered Heatmap ---------------------------------------------------------
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
rownames(mat) <- tx2gene$gene_symbol[match(rownames(mat), tx2gene$gene_id)]
pdf(file = "clustered_heatmap.pdf")
anno <- as.data.frame(colData(vsd))
pheatmap (mat, annotation_col=anno, main = "Top 20 Most Variable Genes",
          color = plasma(255))
dev.off()

# DESeq2 Comparisons -----------------------------------------------------------
comparison <- paste(condition1, "_vs_", condition2, sep = "")
res <-
  results(dds,
          contrast = c("condition", condition1, condition2),
          parallel = TRUE, alpha = 0.05)
gene_synonym <- unique(tx2gene[,-1])
z <- data.frame(res)
z$gene_symbol <-
    gene_synonym$gene_symbol[match(rownames(res), gene_synonym$gene_id)]
z$chr <-
    gene_synonym$chr[match(rownames(res), gene_synonym$gene_id)]
z$start <-
    gene_synonym$start[match(rownames(res), gene_synonym$gene_id)]
z$end <-
    gene_synonym$end[match(rownames(res), gene_synonym$gene_id)]
z$description <-
    gene_synonym$description[match(rownames(res), gene_synonym$gene_id)]
z <- merge(z,txi$abundance, by=0)
setnames(z, "Row.names", "gene_id")
z <-
    z[, !(names(z) %in% row.names(subset(samples,
                                         condition != condition1
                                         & condition != condition2)))]
dir.create(file.path(WD, comparison))
fwrite(z, file = file.path(WD, comparison, paste(comparison, ".res.csv", sep="")))

# Expression Profile Plot -------------------------------------------------
pdf(file=file.path(WD, comparison, paste0(comparison, ".expression_profile.pdf")))
p <- ggplot(z,
            aes(
              x = baseMean + 0.01,
              y = log2FoldChange,
              col = paste(padj < 0.05, is.na(padj)),
              shape = paste(padj < 0.05, is.na(padj))
            ))
p + geom_point() +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(-3,3,1)) +
  scale_color_manual(
    labels = c("padj>0.05", "Below threshold", "padj<0.05"),
    values = c("lightblue", "salmon", "green"),
    guide_legend(title = "")) +
  scale_shape_manual(guide = "none", values = c(0, 1, 2)) +
  labs(title = paste(paste(condition1, "vs.", condition2), sep = ""),
       subtitle = "MA Plot", x = "baseMean") +
  guides(colour = guide_legend(override.aes = list(shape = c(0, 1, 2)))) +
  annotate(geom = "text", y = 3, x = 5,
           label = paste0(nrow(as.data.table(z)[padj <= 0.05 & padj != "NA"
                               & log2FoldChange > 0]), " Genes Upregulated in ",
                          condition1), hjust = "middle") +
  annotate(geom = "text", y = -3, x = 5,
           label = paste0(nrow(as.data.table(z)[padj <= 0.05 & padj != "NA"
                               & log2FoldChange < 0]), " Genes Downregulated in ",
                          condition1), hjust = "middle")
dev.off()

# Volcano Plot ------------------------------------------------------------
pdf(file=file.path(WD, comparison, paste0(comparison, ".volcano_plot.pdf")))
p <-
  ggplot(z, aes(
    x = log2FoldChange,
    y = 1 - padj,
    col = abs(log2FoldChange)))
p + geom_point(size = 1) + scale_y_log10()
dev.off()

# Clustered Top 50 Heatmap ------------------------------------------------
top.count <- 50
top.genes <- z$gene_id[order(z$padj)][1:top.count]
abund <- txi$abundance
abund <- abund[which(rownames(abund) %in% top.genes), ]
abund.scale <- t(scale(t(abund), center = TRUE, scale = TRUE))
rownames(abund.scale) <- tx2gene$gene_symbol[match(rownames(abund.scale),
                                                   tx2gene$gene_id)]
pdf(file = file.path(WD, comparison, paste0(comparison, ".heatmap.pdf")))
pheatmap(mat = abund.scale, main = paste0(comparison, "\n Top ", top.count,
                                          " Significant Genes"),
         annotation_col = anno, angle_col = 45, fontsize_row =6,
         annotation_names_col = FALSE, color = plasma(255))
dev.off()

# GO Enrichment -----------------------------------------------------------
library(gprofiler2)
library(stringr)
library(htmlwidgets)

sig_up <-
    as.data.table(z)[order(padj)][log2FoldChange > 1 & padj<= 0.05 ][, c("gene_id", "gene_symbol", "log2FoldChange", "padj","chr")]
sig_down <-
    as.data.table(z)[order(padj)][log2FoldChange < -1 & padj<= 0.05 ][, c("gene_id", "gene_symbol", "log2FoldChange", "padj","chr")]
sig_genes <- setNames(list(fifelse(as.character(sig_up$gene_symbol) != "",
                                   as.character(sig_up$gene_symbol),
                                   as.character(str_extract(sig_up$gene_id, "^.*(?=(\\.))"))),
                           fifelse(as.character(sig_down$gene_symbol) != "",
                                   as.character(sig_down$gene_symbol),
                                   as.character(str_extract(sig_down$gene_id, "^.*(?=(\\.))")))),
                      c(paste0(comparison, "_up"), paste0(comparison, "_down")))
gost_res <-
    gost(query = sig_genes,
         organism = "hsapiens",
         ordered_query = TRUE,
         multi_query = FALSE,
         evcodes = TRUE,
         exclude_iea = TRUE)
gost_link <-
    gost(query = sig_genes,
         organism = "hsapiens",
         ordered_query = TRUE,
         exclude_iea = TRUE,
         as_short_link = TRUE)
gp <-
    gostplot(gost_res,
             capped = FALSE,
             interactive = TRUE)
htmlwidgets::saveWidget(gp,
                        selfcontained = FALSE,
                        title = comparison,
                        file = file.path(WD, comparison,
                                         paste0(comparison, ".gProfiler_plot.html")))
fwrite(as.data.table(gost_res$result)[order(p_value)],
       file = file.path(WD, comparison, paste0(comparison, ".gProfiler.csv")))
cat(paste0("<html>\n<head>\n<meta http-equiv=\"refresh\" content=\"0; url=",
           gost_link, "\" />", "\n</head>\n<body>\n</body>\n</html>"),
    file=file.path(WD, comparison, paste0(comparison, ".gProfiler.html")))
browseURL(gost_link)


# Reporting Tool
library("ReportingTools")
library("ensembldb")
library("EnsDb.Hsapiens.v86")
htmlRep <- HTMLReport(shortName = paste0(comparison, "_report"),
                      title = comparison,
                      basePath = file.path(paste0(WD,"/", comparison)),
                      reportDirectory = "report")
publish(dds, htmlRep, pvalueCutoff = 0.05, annotation.db = "EnsDb.Hsapiens.v86", factor = colData(dds)$condition)
url <- finish(htmlRep)
browseURL(url)
