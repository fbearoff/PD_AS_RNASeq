# Inputs
WD <- "/home/frank/R_projects/PD_AS_RNASeq"
condition1 <- "10ug_ml_ASPFF"
condition2 <- "10ug_ml_ASM"

condition1_pretty <- "10μg/ml PFF"
condition2_pretty <- "10μg/ml ASM"

# Load Libraries
library(BiocParallel)
library(tximport)
library(readr)
library(DESeq2)
library(data.table)
library(ggplot2)
library(gplots)
library(genefilter)
library(ggrepel)
library(viridis)
library(dplyr)
library(Cairo)
library(ggtext)
library(ComplexHeatmap)
register(MulticoreParam(10))

# TxImport
setwd(WD)
samples <-
  read.table(file.path(WD, "samples.txt"), header = TRUE)
samples$mapped_percent <-
  read.table(file.path(WD, "salmon_output", "mapped_percent.txt"))[, 2]
files <- file.path(WD, "salmon_output", "quants", samples$sample.id, "quant.sf")
names(files) <- samples$sample.id
all(file.exists(files))
tx2gene <- fread(file.path(WD, "tx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2])
txi_abund <- as.data.frame(txi$abundance)

# DESeq2
rownames(samples) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
dds <- DESeq(dds, parallel = TRUE)

# PCA Plot
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd,
  intgroup = "condition",
  returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pdf(file = "pca.pdf")
ggplot(pca_data, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = name), show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  labs(title = "Variance Stabilizing Transformation PCA Plot")
dev.off()

# Clustered Sample Distance Plot
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(colnames(vsd), vsd$condition, sep = ":")
colnames(sample_dist_matrix) <- NULL
colors <- plasma(255)
pdf(file = "clustered_distance.pdf")
ComplexHeatmap::pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors
)
dev.off()

# Clustered Heatmap
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat <- assay(vsd)[top_var_genes, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- tx2gene$gene_symbol[match(rownames(mat), tx2gene$gene_id)]
anno <- as.data.frame(colData(vsd)[, c("condition", "RIN", "DV200")])
pdf(file = "clustered_heatmap.pdf")
ComplexHeatmap::pheatmap(mat,
  annotation_col = anno,
  main = "Top 20 Most Variable Genes",
  color = plasma(255)
)
dev.off()

## merge txi_abund
gene_synonym <- unique(tx2gene[, -1])

txi_abund <- merge(txi_abund,
  gene_synonym,
  by.x = "row.names",
  by.y = "gene_id"
)

txi_abund <- rename(txi_abund, gene_id = Row.names)
fwrite(txi_abund,
  file = file.path(WD, "txi.csv")
)

# DESeq2 Comparisons
comparison <- paste0(condition1, "_vs_", condition2)
res <-
  results(dds,
    contrast = c("condition", condition1, condition2),
    parallel = TRUE, alpha = 0.05
  )
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
  z[, !(names(z) %in% row.names(subset(
    samples,
    condition != condition1 &
      condition != condition2
  )))]
z <- as.data.table(z)
z <- z[chr %in% c(1:22, "MT", "X", "Y")][gene_symbol != ""]
dir.create(file.path(WD, comparison))
fwrite(z,
  file = file.path(
    WD,
    comparison,
    paste0(comparison, ".res.csv")
  )
)

## Expression Profile Plot
# highligts the top 12 transcripts by FC in each direction, baseMean of 10 or more
top_12 <- z[padj <= 0.05][baseMean > 10][order(log2FoldChange, decreasing = TRUE)] %>% slice_head(n = 12)
top_12 <- rbind(
  top_12,
  z[padj <= 0.05][baseMean > 10][order(log2FoldChange, decreasing = TRUE)] %>% slice_tail(n = 12)
)
p <- ggplot(
  z,
  aes(
    x = baseMean + 0.01,
    y = log2FoldChange,
    color = ifelse(is.na(padj),
      padj > 0.05,
      padj < 0.05 &
        abs(log2FoldChange) > 1 &
        baseMean > 10
    ),
    shape = ifelse(is.na(padj),
      padj > 0.05,
      padj < 0.05 &
        abs(log2FoldChange) > 1 &
        baseMean > 10
    )
  )
)
p + geom_point(na.rm = TRUE) +
  geom_text_repel(
    data = top_12,
    aes(label = paste0(
      gene_symbol,
      " (",
      round(log2FoldChange,
        digits = 1
      ),
      ")"
    )),
    show.legend = FALSE,
    color = "black",
    fontface = "bold"
  ) +
  scale_x_log10() +
  scale_y_continuous(trans = "pseudo_log") +
  scale_color_manual(
    labels = c("padj>0.05", "DE Transcripts", "padj=NA"),
    values = plasma(2, begin = .2, end = .7),
    guide_legend(title = "")
  ) +
  scale_shape_manual(
    guide = "none",
    na.value = 2,
    values = c(0, 1, 2)
  ) +
  labs(
    title = "α-Syn PFFs Alter Gene \nExpression in Monocytes",
    x = "baseMean",
    y = "log2FC"
  ) +
  guides(color = guide_legend(override.aes = list(
    shape = c(0, 1, 2),
    size = 4
  ))) +
  annotate(
    geom = "richtext",
    y = max(z$log2FoldChange, na.rm = TRUE) / 2,
    x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
    label = paste0(
      nrow(z[padj <= 0.05 & padj != "NA" &
        log2FoldChange > 1 &
        baseMean > 10]),
      " Genes <i>Upregulated</i> in ", "<b>",
      condition1_pretty, "</b>"
    ),
    hjust = "middle",
    size = 6
  ) +
  annotate(
    geom = "richtext",
    y = min(z$log2FoldChange, na.rm = TRUE) / 2,
    x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
    label = paste0(
      nrow(z[padj <= 0.05 & padj != "NA" &
        log2FoldChange < -1 &
        baseMean > 10]),
      " Genes <i>Downregulated</i> in <b>",
      condition1_pretty, "</b>"
    ),
    hjust = "middle",
    size = 6
  ) +
  theme_minimal() +
  theme(
    plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in"),
    legend.text = element_text(size = 15),
    legend.position = "bottom",
    plot.title = element_text(
      hjust = 0.5,
      size = 40
    )
  )
ggsave(last_plot(),
  file = file.path(
    WD,
    comparison,
    paste0(comparison, ".expression_profile.pdf")
  ),
  device = cairo_pdf,
  width = 26.5,
  height = 14.5,
  unit = "in"
)

## Volcano Plot
pdf(file = file.path(WD, comparison, paste0(comparison, ".volcano_plot.pdf")))
p <-
  ggplot(z, aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = abs(log2FoldChange)
  ))
p + geom_point(
  size = 1,
  na.rm = TRUE
) +
  scale_y_log10() +
  geom_hline(yintercept = 1) +
  scale_color_viridis(
    option = "plasma",
    begin = 0.1,
    end = 0.8,
    direction = -1
  )
dev.off()

## Clustered Top 50 Heatmap (comparison only)
top_count <- 50
top_genes <- z[abs(log2FoldChange) > 1][order(padj)]$gene_id[1:top_count]
abund <- txi$abundance[, !(colnames(txi$abundance) %in%
  row.names(subset(
    samples,
    condition != condition1 &
      condition != condition2
  )))]
abund <- abund[which(rownames(abund) %in% top_genes), ]
rownames(abund) <- tx2gene$gene_symbol[match(
  rownames(abund),
  tx2gene$gene_id
)]
abund_scale <- t(scale(t(abund),
  center = TRUE,
  scale = TRUE
))
anno_col <- samples[which(samples$condition %in% c(condition1, condition2)), ]
anno_col$condition[anno_col$condition == condition1] <- condition1_pretty
anno_col$condition[anno_col$condition == condition2] <- condition2_pretty

col_names <- list(condition = setNames(
  plasma(2, begin = 0.5, end = 0.7),
  c(condition1_pretty, condition2_pretty)
))

ha <- HeatmapAnnotation(
  df = anno_col["condition"],
  show_annotation_name = FALSE,
  border = FALSE,
  annotation_legend_param = list(
    nrow = 2,
    col = plasma(2),
    title_position = "topcenter",
    labels_gp = gpar(fontsize = 8)
  ),
  col = col_names,
  gp = gpar(col = "black")
)

ra <- rowAnnotation(order = anno_text(
  paste0(
    "(",
    rank(z[gene_symbol %in% rownames(abund)]$padj,
      ties.method = "first"
    ),
    ") "
  ),
  gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  just = "center",
  location = 0.5,
  show_name = FALSE
))

hm <- Heatmap(abund_scale,
  column_title = paste0(
    "Top ",
    top_count,
    " DE Genes\nAltered by α-Syn PFFs"
  ),
  column_title_gp = gpar(fontsize = 18),
  top_annotation = ha,
  row_split = 2,
  row_title = NULL,
  right_annotation = ra,
  row_names_gp = gpar(fontsize = 8),
  column_split = 2,
  show_column_names = FALSE,
  show_row_names = TRUE,
  show_parent_dend_line = FALSE,
  col = plasma(255, direction = -1),
  heatmap_legend_param = list(
    title = "Expression",
    legend_direction = "horizontal",
    title_position = "topcenter",
    labels_gp = gpar(fontsize = 8)
  )
)
Cairo(
  file = file.path(WD, comparison, paste0(comparison, ".heatmap.pdf")),
  type = "pdf",
  width = 8,
  height = 10,
  units = "in",
  dpi = "auto",
  family = "Arial"
)
draw(hm,
  merge_legend = TRUE,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom"
)
dev.off()

## GO Enrichment
library(gprofiler2)
library(stringr)
library(htmlwidgets)
library(plyr)
library(forcats)

# MSigDB C2 Canonical pathways
# GSEA c2.cp.v7.5.1 gp__Dg7I_dpHN_NPY
custom_gmt <- "gp__Dg7I_dpHN_NPY"
sig_up <- z[order(padj)][log2FoldChange > 1 & padj <= 0.05 & baseMean > 10][, c("gene_id", "gene_symbol", "log2FoldChange", "padj", "chr")]
sig_down <- z[order(padj)][log2FoldChange < -1 & padj <= 0.05 & baseMean > 10][, c("gene_id", "gene_symbol", "log2FoldChange", "padj", "chr")]
sig_genes <- setNames(
  list(
    fifelse(
      as.character(sig_up$gene_symbol) != "",
      as.character(sig_up$gene_symbol),
      as.character(str_extract(
        sig_up$gene_id,
        "^.*(?=(\\.))"
      ))
    ),
    fifelse(
      as.character(sig_down$gene_symbol) != "",
      as.character(sig_down$gene_symbol),
      as.character(str_extract(
        sig_down$gene_id,
        "^.*(?=(\\.))"
      ))
    )
  ),
  c(
    paste0(comparison, "_up"),
    paste0(comparison, "_down")
  )
)
gost_res <- gost(
  query = sig_genes,
  organism = custom_gmt,
  ordered_query = TRUE,
  multi_query = FALSE,
  evcodes = TRUE,
  exclude_iea = TRUE
)
gost_link <- gost(
  query = sig_genes,
  organism = custom_gmt,
  ordered_query = TRUE,
  exclude_iea = TRUE,
  as_short_link = TRUE
)
gp <- gostplot(gost_res,
  capped = FALSE,
  interactive = TRUE
)
htmlwidgets::saveWidget(gp,
  selfcontained = FALSE,
  title = comparison,
  file = file.path(
    WD, comparison,
    paste0(
      comparison,
      ".gProfiler_plot.html"
    )
  )
)

gdt <- as.data.table(gost_res$result)

# for GSEA GMTs
gdt$term_name <- str_replace_all(
  str_replace_all(
    gdt$term_name,
    "http://www.gsea-msigdb.org/gsea/msigdb/cards/(.*?)_",
    ""
  ),
  "_",
  " "
)
gdt$source <- str_extract(gdt$term_id, "^[^_]*")
fwrite(gdt[, !c("evidence_codes", "parents")][order(p_value)],
  file = file.path(WD, comparison, paste0(comparison, ".gProfiler.csv"))
)
cat(
  paste0(
    "<html>\n<head>\n<meta http-equiv=\"refresh\" content=\"0; url=",
    gost_link,
    "\" />", "\n</head>\n<body>\n</body>\n</html>"
  ),
  file = file.path(
    WD,
    comparison,
    paste0(comparison, ".gProfiler.html")
  )
)
browseURL(gost_link)

# collate lists of pathways for each gene
paths_down <- list()
for (gene in sig_down$gene_symbol) {
  pathway <- gdt[query == paste0(comparison, "_down")]$term_id[grep(gene, gdt[query == paste0(comparison, "_down")]$intersection)]
  paths_down[[gene]] <- pathway
}
fwrite(lapply(paths_down, `length<-`, max(lengths(paths_down))),
  file = paste(WD, comparison, "pathways_down.csv", sep = "/")
)
paths_up <- list()
for (gene in sig_up$gene_symbol) {
  pathway <- gdt[query == paste0(comparison, "_up")]$term_id[grep(gene, gdt[query == paste0(comparison, "_up")]$intersection)]
  paths_up[[gene]] <- pathway
}
fwrite(lapply(paths_up, `length<-`, max(lengths(paths_up))),
  file = paste(WD, comparison, "pathways_up.csv", sep = "/")
)

# gem creation
gem <- gdt[, c("query", "term_id", "term_name", "p_value", "intersection", "source")]
colnames(gem) <- c("query", "GO.ID", "Description", "p.Val", "Genes", "source")
gem$fdr <- gem$p.Val
gem$phenotype <- "+1"

gem %>%
  group_by(query) %>%
  group_walk(~
    write.table(data.frame(.x[, c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes", "source")]),
      file = paste0(comparison, "/", unique(.y$query), "_gem.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    ))

# return pathway genes as list per pathway filtered by direction and data source
query_gost_genes <- function(direction, data_source) {
  query_list <- setNames(
    strsplit(gdt[query == paste0(comparison, "_", direction) & source == data_source, intersection], ","),
    gdt[query == paste0(comparison, "_", direction) & source == data_source, term_name]
  )
  print(query_list)
}

## Pathways figure
pathways <- gdt[, c("term_name", "p_value", "term_size", "intersection_size", "source", "query")][order(query, p_value)]
plot_order <- rownames(pathways[query == paste0(comparison, "_down")])
plot_order <- append(
  plot_order,
  rownames(pathways[query == paste0(comparison, "_up")])
)
pathways$plot_order <- as.numeric(plot_order)
pathways[, `:=`(log_p_value = log(p_value))]
pathways[query == paste0(comparison, "_up"), `:=`(log_p_value = -1 * log_p_value)]

pathway_n <- 30

path_plot <- ggplot(
  pathways[plot_order <= pathway_n],
  aes(
    x = log_p_value,
    y = fct_rev(as.factor(plot_order)),
    fill = source
  )
) +
  geom_col(
    position = position_stack(),
    color = "white"
  ) +
  geom_richtext(
    data = pathways[plot_order <= pathway_n & log_p_value > 1],
    aes(
      label = paste0(
        term_name,
        " (**",
        intersection_size,
        "/",
        term_size,
        "**)"
      ),
      hjust = "left",
      vjust = 0.5
    ),
    position = position_stack(vjust = 0),
    size = 4,
    fill = NA,
    color = revalue(
      pathways[plot_order <= pathway_n & log_p_value > 1]$source,
      c(
        "REACTOME" = "black",
        "WP" = "black",
        "KEGG" = "white",
        "PID" = "black",
        "NABA" = "white",
        "BIOCARTA" = "black",
        "SIG" = "black"
      )
    ),
    label.color = NA
  ) +
  geom_richtext(
    data = pathways[plot_order <= pathway_n & log_p_value < 1],
    aes(
      label = paste0(
        term_name,
        " (**",
        intersection_size,
        "/",
        term_size,
        "**)"
      ),
      hjust = "right",
      vjust = 0.5
    ),
    position = position_dodge(width = 1),
    size = 4,
    fill = NA,
    color = "black",
    label.color = NA
  ) +
  scale_fill_manual(values = plasma(length(levels(as.factor(pathways$source[1:pathway_n]))))) +
  scale_x_continuous(name = "-log p-value", limits = c(-25, NA)) +
  scale_y_discrete(name = "Pathway Rank") +
  theme_classic() +
  labs(title = paste0(
    "Top ",
    pathway_n,
    " Pathways Enriched by α-Syn PFFs"
  )) +
  theme(
    plot.margin = grid::unit(c(0.5, 0.5, 0, 0.5), "in"),
    plot.title = element_text(
      hjust = 0.5,
      size = 50
    ),
    legend.justification = "top"
  )

# Most common genes in results
library(grid)
library(gridExtra)
library(ggplotify)

bar_up <- as.data.table(sort(
  table(unlist(strsplit(
    gdt[query == paste0(comparison, "_up"), intersection],
    ","
  ))),
  decreasing = TRUE
)[1:10])

bar_down <- as.data.table(sort(
  table(unlist(strsplit(
    gdt[query == paste0(comparison, "_down"), intersection],
    ","
  ))),
  decreasing = TRUE
)[1:10])

bar_up_grob <- as.grob(ggplot(bar_up, aes(
  x = reorder(V1, -N),
  y = N,
  fill = N
)) +
  theme_minimal() +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(
    option = "plasma",
    direction = -1
  ) +
  labs(
    title = paste0(
      "Most Frequently Occuring Pathway Genes<br> *Upregulated*  by ",
      condition1_pretty
    ),
    y = "Pathway Occurrences"
  ) +
  theme(
    aspect.ratio = 3 / 5,
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank(),
    plot.title = element_markdown(
      hjust = 0.5,
      size = 15,
      lineheight = 1.1
    )
  ))
bar_down_grob <- as.grob(ggplot(bar_down, aes(
  x = reorder(V1, -N),
  y = N,
  fill = N
)) +
  theme_minimal() +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis(
    option = "plasma",
    direction = -1
  ) +
  labs(
    title = paste0(
      "Most Frequently Occuring Pathway Genes<br>*Downregulated*  by ",
      condition1_pretty
    ),
    y = "Pathway Occurrences"
  ) +
  theme(
    aspect.ratio = 3 / 5,
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = element_blank(),
    plot.title = element_markdown(
      hjust = 0.5,
      size = 15,
      lineheight = 1.1
    )
  ))
blank_plot <- ggplot() +
  geom_blank(aes(1, 1)) +
  theme_void()
path_plot +
  annotation_custom(
    as.grob(grid.arrange(blank_plot,
      arrangeGrob(bar_down_grob,
        blank_plot,
        blank_plot,
        blank_plot,
        bar_up_grob,
        ncol = 5
      ),
      ncol = 1
    )),
    xmin = -25,
    xmax = Inf,
    ymin = 0
  )
ggsave(last_plot(),
  file = paste(WD,
    comparison,
    "pathways.pdf",
    sep = "/"
  ),
  device = cairo_pdf,
  width = 26,
  height = 14.5,
  unit = "in"
)
