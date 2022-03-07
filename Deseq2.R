# Inputs ------------------------------------------------------------------
WD <- "/home/frank/R_projects/FB38_June2021_30-506016007"
condition1 <- "10ug_ml_ASPFF"
condition2 <- "10ug_ml_ASM"

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
library(dplyr)
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
tx2gene <- fread(file.path(WD, "tx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# DESeq2 ------------------------------------------------------------------
rownames(samples) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, samples, ~ condition)
dds <- DESeq(dds, parallel = TRUE)

# PCA Plot ------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition", "RIN", "DV200"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = "pca.pdf")
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 2) + geom_text_repel(aes(label = name), show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + labs(title = "Variance Stabilizing Transformation PCA Plot")
dev.off()

# Clustered Sample Distance Plot ------------------------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$condition, sep = ":")
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
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- tx2gene$gene_symbol[match(rownames(mat), tx2gene$gene_id)]
anno <- as.data.frame(colData(vsd)[, c("condition", "RIN", "DV200")])
pdf(file = "clustered_heatmap.pdf")
pheatmap(mat, annotation_col = anno, angle_col = 45,  main = "Top 20 Most Variable Genes",
          color = plasma(255))
dev.off()

# DESeq2 Comparisons -----------------------------------------------------------
comparison <- paste(condition1, "_vs_", condition2, sep = "")
res <-
  results(dds,
          contrast = c("condition", condition1, condition2),
          parallel = TRUE, alpha = 0.05)
gene_synonym <- unique(tx2gene[, -1])
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
z <- merge(z, txi$abundance, by = 0)
setnames(z, "Row.names", "gene_id")
z <-
    z[, !(names(z) %in% row.names(subset(samples,
                                         condition != condition1
                                         & condition != condition2)))]
z  <- as.data.table(z)
dir.create(file.path(WD, comparison))
fwrite(z, file = file.path(WD, comparison, paste(comparison, ".res.csv", sep = "")))

# Expression Profile Plot -------------------------------------------------
#highligts the top 5 transcripts by FC in each direction
top_10  <- z[padj <= 0.05][order(log2FoldChange, decreasing = TRUE)] %>% slice_head(n = 5)
top_10  <- rbind(top_10, z[padj <= 0.05][order(log2FoldChange, decreasing = TRUE)] %>% slice_tail(n = 5))
p <- ggplot(z,
            aes(
              x = baseMean + 0.01,
              y = log2FoldChange,
              color = paste(padj < 0.05, is.na(padj)),
              shape = paste(padj < 0.05, is.na(padj))
            ))
p + geom_point(na.rm = TRUE) +
    geom_text_repel(data = top_10,
                    aes(label = gene_symbol),
                    show.legend = FALSE,
                    color = "black") +
    scale_x_log10() +
    scale_y_continuous() +
    scale_color_manual(labels = c("padj>0.05", "Below threshold", "padj<0.05"),
    values = plasma(4),
    guide_legend(title = "")) +
    scale_shape_manual(guide = "none",
                       values = c(0, 1, 2)) +
    labs(title = paste(paste(condition1, "vs.", condition2), sep = ""),
         subtitle = "MA Plot", x = "baseMean") +
    guides(color = guide_legend(override.aes = list(shape = c(0, 1, 2)))) +
    annotate(geom = "text",
             y = max(z$log2FoldChange, na.rm = TRUE) / 2,
             x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
             label = paste0(nrow(z[padj <= 0.05 & padj != "NA" &
                                 log2FoldChange > 0]),
                            " Genes Upregulated in ",
                            condition1),
             hjust = "middle") +
    annotate(geom = "text",
             y = min(z$log2FoldChange, na.rm = TRUE) / 2,
             x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
             label = paste0(nrow(z[padj <= 0.05 & padj != "NA" &
                                 log2FoldChange < 0]),
                            " Genes Downregulated in ",
                            condition1),
             hjust = "middle") +
    theme_minimal() +
    theme(plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in"))
ggsave(last_plot(), file = file.path(WD, comparison, paste0(comparison, ".expression_profile.pdf")))

# Volcano Plot ------------------------------------------------------------
pdf(file = file.path(WD, comparison, paste0(comparison, ".volcano_plot.pdf")))
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
         annotation_col = anno, angle_col = 45, fontsize_row = 6,
         annotation_names_col = FALSE, color = plasma(255))
dev.off()

# GO Enrichment -----------------------------------------------------------
library(gprofiler2)
library(stringr)
library(htmlwidgets)
library(plyr)
library(ggtext)
library(forcats)


#MSigDB Canonical pathways v7.4 unique string
custom_gmt <- 'gp__IbKR_F2Rn_Tss'
sig_up <- z[order(padj)][log2FoldChange > 1 & padj <= 0.05 & chr %in% c(1:22, "X", "Y", "MT")][, c("gene_id", "gene_symbol", "log2FoldChange", "padj", "chr")]
sig_down <- z[order(padj)][log2FoldChange < -1 & padj <= 0.05 & chr %in% c(1:22, "X", "Y", "MT")][, c("gene_id", "gene_symbol", "log2FoldChange", "padj", "chr")]

sig_genes <- setNames(list(
                           fifelse(as.character(sig_up$gene_symbol) != "",
                                   as.character(sig_up$gene_symbol),
                                   as.character(str_extract(sig_up$gene_id, "^.*(?=(\\.))"))),
                           fifelse(as.character(sig_down$gene_symbol) != "",
                                   as.character(sig_down$gene_symbol),
                                   as.character(str_extract(sig_down$gene_id, "^.*(?=(\\.))")))),
                      c(paste0(comparison, "_up"),
                        paste0(comparison, "_down")))
gost_res <- gost(query = sig_genes,
                 organism = custom_gmt,
                 ordered_query = TRUE,
                 multi_query = FALSE,
                 evcodes = TRUE,
                 exclude_iea = TRUE)

gost_link <- gost(query = sig_genes,
                  organism = custom_gmt, #usually "hsapaiens"
                  ordered_query = TRUE,
                  exclude_iea = TRUE,
                  as_short_link = TRUE)
gp <- gostplot(gost_res,
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
           gost_link,
           "\" />", "\n</head>\n<body>\n</body>\n</html>"),
    file = file.path(WD,
                   comparison,
                   paste0(comparison, ".gProfiler.html")))
browseURL(gost_link)

gost_res$result$term_name <- str_replace_all(str_replace_all(gost_res$result$term_name,
                                                                  'http://www.gsea-msigdb.org/gsea/msigdb/cards/(.*?)_',
                                                                  ""),
                                                  "_", " ")
gost_res$result$source  <- str_extract(gost_res$result$term_id,
                                            "^[^_]*")
#collate lists of pathways for each gene
paths_down <- list()
for (gene in sig_down$gene_symbol) {
    pathway <- gost_res$result$term_id[grep(gene, gost_res$result$intersection)]
    paths_down[[gene]]  <-  pathway
}
fwrite(lapply(paths_down, `length<-`, max(lengths(paths_down))), file = paste(WD, comparison, "pathways_down.csv", sep = "/"))
paths_up <- list()
for (gene in sig_up$gene_symbol) {
    pathway <- gost_res$result$term_id[grep(gene, gost_res$result$intersection)]
    paths_up[[gene]]  <-  pathway
}
fwrite(lapply(paths_up, `length<-`, max(lengths(paths_up))), file = paste(WD, comparison, "pathways_up.csv", sep = "/"))

gem <- gost_res$result[, c("query", "term_id", "term_name", "p_value", "intersection", "source")]
colnames(gem) <- c("query", "GO.ID", "Description", "p.Val", "Genes", "source")
gem$FDR <- gem$p.Val
gem$Phenotype <- "+1"

gem %>% group_by(query) %>%
  group_walk(~
    write.table(data.frame(.x[, c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes", "source")]),
                file = paste0(comparison, "/", unique(.y$query), "_gem.txt"),
                sep = "\t", quote = F, row.names = F))


#return pathway genes as list per pathway filtered by direction and data source
query_gost_genes  <- function(direction, data_source) {
    query_list <- setNames(strsplit(as.data.table(gost_res$result)[query == paste0(comparison, "_", direction) & source == data_source, intersection], ","),
                           as.data.table(gost_res$result)[query == paste0(comparison, "_", direction) & source == data_source, term_name])
    print(query_list)
}

pathways <- as.data.table(gost_res$result)[, c("term_name", "p_value", "term_size", "intersection_size", "source", "query")][order(query, p_value)]
plot_order <- rownames(pathways[query == paste0(comparison, "_down")])
plot_order <- append(plot_order, rownames(pathways[query == paste0(comparison, "_up")]))
pathways$plot_order <- as.numeric(plot_order)
pathways[, `:=`(log_p_value = log(p_value))]
pathways[query == paste0(comparison, "_up"), `:=`(log_p_value = -1 * log_p_value)]

pathway_n <- 30

path_plot <- ggplot(pathways[plot_order <= pathway_n], aes(x = log_p_value, y = fct_rev(as.factor(plot_order)), fill = source)) +
    geom_col(position = position_stack(),
             color = "white") +
    geom_richtext(data = pathways[plot_order <= pathway_n & log_p_value > 1],
                  aes(label = paste0(term_name, " (**", intersection_size, "/", term_size, "**)"),
                      hjust = "left",
                      vjust = 0.5),
                  position = position_stack(vjust = 0),
                  size = 3,
                  fill = NA,
                  color = revalue(pathways[plot_order <= pathway_n & log_p_value > 1]$source,
                                c("REACTOME" = "black", "WP" = "black", "KEGG" = "white", "PID" = "black", "NABA" = "white", "BIOCARTA" = "black", "SIG" = "black")),
              label.color = NA) +
    geom_richtext(data = pathways[plot_order <= pathway_n & log_p_value < 1],
                  aes(label = paste0(term_name, " (**", intersection_size, "/", term_size, "**)"),
                      hjust = "right",
                      vjust = 0.5),
                  position = position_dodge(width = 1),
                  size = 3,
                  fill = NA,
                  color = "black",
                  label.color = NA) +
    scale_fill_manual(values = plasma(4)) +
    scale_x_continuous(name = "-log p-value", limits = c(-25, NA)) +
    scale_y_discrete(name = "Pathway Rank")  +
    theme_classic() +
    labs(title = paste0("Top ", pathway_n, " Enriched Pathways in ", comparison)) +
    theme(plot.margin = grid::unit(c(0.5, 0.5, 0, 0.5), "in"),
          plot.title = element_text(hjust = 0.5),
          legend.justification = "top")

# Most common genes in results
library(grid)
library(gridExtra)
library(ggplotify)

bar_up <- as.data.table(sort(table(unlist(strsplit(as.data.table(gost_res$result)[query == paste0(comparison, "_up"), intersection],
                                                   ","))),
                             decreasing = TRUE)[1:10])
bar_down <- as.data.table(sort(table(unlist(strsplit(as.data.table(gost_res$result)[query == paste0(comparison, "_down"), intersection],
                                                   ","))),
                             decreasing = TRUE)[1:10])
bar_up_grob <- as.grob(ggplot(bar_up, aes(x = reorder(V1, -N), y = N, fill = N)) +
                       theme_minimal() +
                       geom_col(show.legend = FALSE) +
                       scale_fill_viridis(option = "plasma",
                                          direction = -1) +
                       labs(title = paste0("Most Frequently Occuring Pathway Genes<br> *Upregulated*  in ", condition1),
                            y = "Pathway Occurrences") +
                       theme(aspect.ratio = 3 / 5,
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                             panel.grid = element_blank(),
                             plot.title = element_markdown(hjust = 0.5, size = 10)))
bar_down_grob <- as.grob(ggplot(bar_down, aes(x = reorder(V1, -N), y = N, fill = N)) +
                       theme_minimal() +
                       geom_col(show.legend = FALSE) +
                       scale_fill_viridis(option = "plasma",
                                          direction = -1) +
                       labs(title = paste0("Most Frequently Occuring Pathway Genes<br>*Downregulated*  in ", condition1),
                            y = "Pathway Occurrences") +
                       theme(aspect.ratio = 3 / 5,
                             axis.title.x = element_blank(),
                             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                             panel.grid = element_blank(),
                             plot.title = element_markdown(hjust = 0.5, size = 10)))
blankPlot <- ggplot() +
    geom_blank(aes(1, 1)) +
    theme_void()
path_plot + annotation_custom(as.grob(grid.arrange(blankPlot, arrangeGrob(bar_down_grob, blankPlot, blankPlot, blankPlot, bar_up_grob, ncol = 5), ncol = 1)), xmin = -25, xmax = Inf, ymin = 0)
ggsave(last_plot(), file = paste(WD, comparison, "pathways.pdf", sep = "/"))

# Plot Gene -----------------------------------------------------------
plot_gene <- function(gene_name) {
    gene_name  <- toupper(gene_name)
    gene_id  <- as.data.table(tx2gene)[gene_symbol == gene_name & chr %in% c(1:22, "X", "Y", "MT"), gene_id][1]
    plot_data <- plotCounts(dds, gene = as.character(gene_id), returnData = TRUE)
    ggplot(plot_data, aes(x = condition, y = count, color = condition)) +
           geom_point(show.legend = TRUE) +
           scale_color_manual(values = plasma(5)) +
           ggtitle(paste0(gene_name, "\n", gene_id)) +
           theme(plot.title = element_text(hjust = .5))
}

plot_gene("SLC18A2")

ggsave(last_plot(), file = "pathways.pdf")

for (gene in rownames(mat)) {
    pdf(file = paste0(WD, "/expression_plots/", gene, ".pdf"))
    final_plot <- plot_gene(gene)
    print(final_plot)
    dev.off()
}
