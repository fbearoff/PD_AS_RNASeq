# inputs
WD <- "/home/frank/R_projects/PD_AS_RNASeq"
condition1 <- "10ug_ml_ASPFF"
condition2 <- "10ug_ml_ASM"

# libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(colorRamps)
library(ggtext)

# cibersort
comparison <- paste0(condition1, "_vs_", condition2)
cibersort <- fread(file.path(WD, comparison, "cibersort_import.csv"))
cibersort <- cibersort[PFF + ASM != 0]

cibersort_long <- as.data.table(gather(cibersort,
  Condition,
  Percent,
  ASM:PFF,
  factor_key = TRUE
))

cibersort_long <- cibersort_long[order(Condition, cell_cat, cell_type)]
cibersort <- cibersort[order(cell_cat, cell_type)]
cibersort_long$cell_type <- factor(cibersort_long$cell_type,
  levels = cibersort$cell_type
)

ggplot(
  data = cibersort_long,
  aes(
    x = Condition,
    y = Percent,
    fill = cell_type
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  scale_fill_manual(values = primary.colors(n = length(rownames(cibersort)))) +
  scale_y_continuous(expand = c(0, 0.1)) +
  geom_richtext(aes(label = if_else(Percent > 2,
    paste0(
      cell_type,
      "<b> (",
      round(Percent, 1),
      "%)</b>"
    ),
    NULL
  )),
  position = position_stack(vjust = 0.5),
  size = 5,
  label.color = NA,
  na.rm = TRUE,
  show.legend = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      size = 20,
      hjust = 0.5
    ),
    legend.box.spacing = unit(-0.5, "in"),
    legend.text = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 20),
    plot.margin = grid::unit(c(0.5, 0, 0, 0), "in")
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  ggtitle("PFFs Induce Macrophage Transcriptional Profiles") +
  geom_segment(
    data = cibersort,
    color = primary.colors(length(rownames(cibersort))),
    aes(
      x = 1 + 0.7 / 2,
      y = 100 - cumsum(ASM),
      xend = 2 - 0.7 / 2,
      yend = 100 - cumsum(PFF),
      color = cell_type
    )
  )

ggsave(last_plot(),
  file = file.path(
    WD,
    comparison,
    paste0(comparison, ".cibersort.pdf")
  ),
  device = cairo_pdf,
  width = 20,
  height = 20,
  unit = "in"
)
