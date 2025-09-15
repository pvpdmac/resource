# Setup chunk -------------------------------------------------------------

library(DESeq2)
library(patchwork)
library(UpSetR)
library(tidyverse)

source("Scripts/_pca_functions.R")

project_colours <- list(
  up_down = c("up" = "#f8766d", "down" = "#00bfc4", "not_sig" = "#c0c0c0"),
  vaccines = c(
    "VEH" = "#C0C0C0",
    "HBV" = "#0000FF",
    "BCG" = "#FF0000",
    "HBV+BCG" = "#B33FF8"
  ),
  DOL = c(
    "DOL0" = "#440154",
    "DOL1" = "#31688e",
    "DOL3" = "#35b779",
    "DOL7" = "#fde725"
  )
)

data(reactome_categories_HSA_L1_L2, package = "clusterPathways")

out_dir <- "Results_condition"
dir.create(out_dir)

# Setting the seed helps the networks and labeled PCA plot be consistent between
# runs
set.seed(123)


# Define helper functions -------------------------------------------------

volcano_plot <- function(
    x,
    main_title,
    subtitle = "default",
    prop_genes = 0.05,
    cutoff_p = 0.05,
    cutoff_fc = 1.5,
    xlims = NULL,
    ylims = NULL
) {

  y <- x %>%
    mutate(
      fc_dir = case_when(
        log2FoldChange < -log2(cutoff_fc) & padj < cutoff_p ~ "down",
        log2FoldChange > log2(cutoff_fc) & padj < cutoff_p ~ "up",
        TRUE ~ "not_sig"
      )
    )

  if (!is.null(prop_genes)) {
    genes_to_label <- y %>%
      filter(fc_dir %in% c("up", "down")) %>%
      group_by(fc_dir) %>%
      arrange(desc(abs(log2FoldChange)), padj) %>%
      slice_head(prop = prop_genes) %>%
      ungroup() %>%
      pull(hgnc_symbol)
  }

  z <- mutate(
    y,
    label_genes = if_else(
      hgnc_symbol %in% genes_to_label,
      hgnc_symbol,
      NA_character_
    )
  )

  if (subtitle == "default") {
    my_subtitle <- str_glue(
      "Down: {nrow(filter(z, fc_dir == 'down'))}; ",
      "Up: {nrow(filter(z, fc_dir == 'up'))}"
    )
  } else {
    my_subtitle <- subtitle
  }

  ggplot(
    z,
    aes(log2FoldChange, -log10(pvalue), fill = fc_dir, label = label_genes)
  ) +
    geom_point(size = 3, pch = 21, alpha = 0.5) +
    geom_hline(
      yintercept = -log10(cutoff_p),
      linetype = "dashed",
      alpha = 0.8
    ) +
    geom_vline(
      xintercept = c(-log2(cutoff_fc), log2(cutoff_fc)),
      linetype = "dashed",
      alpha = 0.8
    ) +
    ggrepel::geom_label_repel(fill = "white", size = 6) +
    scale_fill_manual(
      values = c("up" = "green", "down" = "yellow", "not_sig" = "#c0c0c0"),
      guide = NULL
    ) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    labs(
      title = main_title,
      x = "log2FC",
      y = "-log10(P Value)",
      subtitle = my_subtitle
    ) +
    theme_bw() +
    theme(
      title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 22, face = "bold"),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 22)
    )
}

rpa_wrapper <- function(x, gene_universe = NULL, species = "human") {
  x <- na.omit(x)

  message(str_glue(
    "Testing {length(x)} genes..."
  ), appendLF = FALSE)

  output <- ReactomePA::enrichPathway(
    gene = na.omit(x),
    organism = species,
    universe = gene_universe
  ) %>%
    pluck("result") %>%
    as_tibble() %>%
    janitor::clean_names()

  message("Done")
  return(output)
}

pathway_plot <- function(data, rows) {
  ggplot(
    mutate(data, direction = factor(direction, c("up", "down"))),
    aes(
      numerator,
      description,
      shape = direction,
      fill = rpa_gr,
      size = -log10(p_adjust)
    )
  ) +
    facet_grid(rows = vars({{rows}}), scales = "free_y", space = "free_y") +
    geom_point(colour = "black") +
    scale_shape_manual(values = c("up" = 24, "down" = 25)) +
    scale_fill_viridis_c() +
    scale_size_continuous(range = c(3, 7)) +
    theme_bw(base_size = 16) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Gene ratio",
      size = bquote(-log[10](italic(p[adjusted]))),
      shape = "Direction"
    ) +
    theme(
      text = element_text(family = "Helvetica", colour = "black"),
      axis.text.y = ggtext::element_markdown(size = 16),
      axis.text.x = element_text(size = 16, angle = 30, hjust = 1, colour = "black"),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0, size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    guides(
      size = guide_legend(order = 1, override.aes = list(shape = 24, fill = "white")),
      fill = guide_colorbar(order = 2),
      shape = guide_legend(order = 3, override.aes = list(size = 5))
    )
}

pathway_plot2 <- function(data, rows) {
  ggplot(
    mutate(data, direction = factor(direction, c("up", "down"))),
    aes(
      numerator,
      description,
      shape = direction,
      fill = rpa_gr,
      size = -log10(p_adjust)
    )
  ) +
    ggforce::facet_col(
      facets = vars({{rows}}),
      scales = "free_y",
      space = "free"
    ) +
    geom_point(colour = "black") +
    scale_shape_manual(values = c("up" = 24, "down" = 25)) +
    scale_fill_viridis_c() +
    scale_size_continuous(range = c(3, 7)) +
    theme_bw(base_size = 16) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Gene ratio",
      size = bquote(-log[10](italic(p[adjusted]))),
      shape = "Direction"
    ) +
    theme(
      text = element_text(family = "Helvetica", colour = "black"),
      axis.text.y = ggtext::element_markdown(size = 16, colour = "black"),
      axis.text.x = element_text(size = 16, angle = 30, hjust = 1, colour = "black"),
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16, vjust = 0),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16)
    ) +
    guides(
      size = guide_legend(order = 1, override.aes = list(shape = 24, fill = "white")),
      fill = guide_colorbar(order = 2),
      shape = guide_legend(order = 3, override.aes = list(size = 5))
    )
}


# Loading counts, sample, and supporting data -----------------------------

samples <-
  read_csv("Metadata/HIPC_P3_TissueConstruct_rnaseq_summary_response.csv") %>%
  filter(qc_notes == "Pass") %>%
  mutate(
    condition = factor(
      str_replace(condition, "Combo", "HBV+BCG"),
      levels = c("VEH", "HBV", "BCG", "HBV+BCG")
    ),
    date_submitted = as.factor(date_submitted),
    tissue_construct_id = as.factor(tissue_construct_id)
  )

count_matrix <- read.csv(
  "count_matrix_HIPC_tissue_construct.csv",
  row.names = 1,
  check.names = FALSE
) %>%
  dplyr::select(all_of(samples$library_name)) # Ensure consistent order

biomart_table <- read_csv("biomart_table_human.csv", col_types = cols())


# Run DESeq ---------------------------------------------------------------

cds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = as.data.frame(samples),
  design = ~tissue_construct_id + condition
)

dds <- DESeq(cds, parallel = TRUE)
vsd <- vst(dds, blind = FALSE)


# Create PCA plots --------------------------------------------------------

pca_data <- plotPCA(
  object = vsd,
  intgroup = c(
    "condition",
    "tissue_construct_id",
    "date_submitted",
    "location"
  ),
  returnData = TRUE
)

percent_var <- map2_chr(
  .x = c("PC1: ", "PC2: "),
  .y = paste0(round(100 * attr(pca_data, "percentVar")), "% Variance"),
  ~paste0(.x, .y)
)

# Condition and batch
ggplot(
  pca_data,
  aes(PC1, PC2, fill = condition, shape = date_submitted)
) +
  geom_point(colour = "black", size = 4) +
  scale_fill_manual(values = project_colours[["vaccines"]]) +
  scale_shape_manual(values = c(21, 24, 22)) +
  labs(
    x = percent_var[1],
    y = percent_var[2],
    fill = "Condition",
    shape = "Batch"
  ) +
  theme_bw(base_size = 18) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey"))
  )
ggsave(
  file.path(out_dir, "pca_plot_condition_batch.png"),
  width = 13,
  height = 6
)

# Condition and batch, with labels
ggplot(
  pca_data,
  aes(PC1, PC2, fill = condition, shape = date_submitted)
) +
  geom_point(colour = "black", size = 4) +
  ggrepel::geom_text_repel(aes(
    label = tissue_construct_id,
    colour = condition
  )) +
  scale_fill_manual(values = project_colours[["vaccines"]]) +
  scale_colour_manual(values = project_colours[["vaccines"]]) +
  scale_shape_manual(values = c(21, 24, 22)) +
  labs(
    x = percent_var[1],
    y = percent_var[2],
    fill = "Condition",
    shape = "Batch"
  ) +
  theme_bw(base_size = 18) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21)),
    shape = guide_legend(override.aes = list(fill = "grey")),
    colour = guide_none()
  )
ggsave(
  file.path(out_dir, "pca_plot_condition_batch_labeled.png"),
  width = 13,
  height = 6
)

# Location
ggplot(pca_data, aes(PC1, PC2, fill = location)) +
  geom_point(colour = "black", size = 4, pch = 21) +
  labs(
    x = percent_var[1],
    y = percent_var[2],
    fill  = "Location"
  ) +
  theme_bw(base_size = 18)
ggsave(
  file.path(out_dir, "pca_plot_location.png"),
  width = 13,
  height = 6
)


# Check variable and CIBSERSORT contributions -----------------------------

cibersort_results <-
  read_tsv("Cibersort/cibersort_local_results_transformed_counts.txt") %>%
  dplyr::rename("library_name" = Mixture) %>%
  janitor::clean_names() %>%
  dplyr::select(-c(p_value, rmse, correlation)) %>%
  column_to_rownames("library_name")

var_con_results <- model_var_contributions(
  deseq_object = vsd,
  cell_props = cibersort_results,
  vars_of_interest = c("condition", "date_submitted", "location"),
  n_genes = 500,
  n_PCs = 10
)

var_con_plot <- var_con_results$results %>%
  tRavis::tr_sort_alphanum(., "term") %>%
  mutate(
    PC = factor(PC, levels = rev(paste0("PC", 1:10))),
    term = fct_inorder(term)
  ) %>%
  ggplot(., aes(term, PC, colour = log_p_value)) +
  geom_point(size = 6) +
  scale_colour_gradientn(colours = viridisLite::plasma(256)) +
  theme_bw(base_size = 18) +
  labs(
    y = "Gene Expression\nPrinciple Components",
    x = NULL,
    colour = bquote(~-Log[10](~italic(p)[val])),
    caption = "Dots represent combinations of variables/PCs with p-value < 0.05"
  ) +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.text.x = element_text(
      family = "Helvetica",
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(family = "Helvetica", colour = "black"),
    legend.position = "top",
    plot.caption = element_text(hjust = 0)
  )

percent_vars_tbl <- var_con_results$percent_vars %>%
  enframe("PC", "percent_var") %>%
  mutate(PC = factor(PC, levels = rev(paste0("PC", 1:10))))

percent_var_plot <- ggplot(percent_vars_tbl, aes(percent_var, PC)) +
  geom_col() +
  scale_y_discrete(expand = expansion(add = 0.6), position = "right") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  theme_bw(base_size = 18) +
  labs(x = "Percent Variance") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.text.x = element_text(family = "Helvetica", colour = "black"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank()
  )

var_con_plot + percent_var_plot + plot_layout(widths = c(4, 1))
ggsave(
  file.path(out_dir, "variables_PC_correlations_de_analysis_condition.png"),
  width = 12,
  height = 7
)


# Check Cook's distance ---------------------------------------------------

cooks_data <- log10(assays(dds)[["cooks"]]) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    starts_with("RNA"),
    names_to = "library_name",
    values_to = "value"
  ) %>%
  left_join(., samples, by = "library_name") %>%
  mutate(library_name = fct_inorder(library_name))

ggplot(cooks_data, aes(value, library_name, fill = condition)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.5) +
  scale_fill_manual(values = project_colours[["vaccines"]]) +
  theme_classic(base_size = 16) +
  labs(x = NULL, y = NULL, fill = "Condition")
ggsave(
  file.path(out_dir, "cooks_distance.png"),
  width = 10,
  height = 22
)


# Get DE genes ------------------------------------------------------------

de_comparisons <- list(
  "condition_HBV_vs_VEH" = c("condition", "HBV", "VEH"),
  "condition_BCG_vs_VEH" = c("condition", "BCG", "VEH"),
  "condition_HBV+BCG_vs_VEH" = c("condition", "HBV+BCG", "VEH"),
  "condition_BCG_vs_HBV" = c("condition", "HBV", "BCG")
)

de_genes <-  de_comparisons %>% map(
  ~results(dds, contrast = .x) %>%
    tRavis::tr_clean_deseq2_result() %>%
    left_join(., biomart_table, by = c("gene" = "ensembl_gene_id"))
)

iwalk(
  de_genes,
  ~write_csv(
    x = .x,
    file = file.path(
      out_dir,
      paste0("deg_", str_replace(.y, "\\+", "_"), ".csv")
    )
  )
)

de_genes_plotdata <- de_genes %>% map_df(
  ~tibble(
    "total" = nrow(.x),
    "up" = filter(.x, log2FoldChange > 0) %>% nrow(),
    "down" = filter(.x, log2FoldChange < 0) %>% nrow()
  ), .id = "comp"
) %>%
  pivot_longer(total:down, names_to = "type", values_to = "num") %>%
  filter(type != "total") %>%
  mutate(
    numerator = fct_inorder(str_remove_all(comp, "^condition_|_vs_[A-Z]{3}$")),
    denominator = if_else(
      str_detect(comp, "vs_VEH"),
      "<i>cf.</i> Vehicle",
      "<i>cf.</i> HBV"
    ),
    denominator = fct_inorder(
      paste0("<i>cf.</i> ", str_extract(comp, "[A-Z]{3}$"))
    )
  )

ggplot(
  de_genes_plotdata,
  aes(numerator, num, fill = type, label = num)
) +
  facet_grid(cols = vars(denominator), scales = "free", space = "free", ) +
  geom_col(colour = "black", width = 0.8) +
  geom_text(
    size = 7,
    position = position_stack(vjust = 0.5),
    colour = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = project_colours[["up_down"]]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 22) +
  labs(
    x = "Test condition",
    y = "No. genes",
    fill = "Direction"
  ) +
  theme(
    text = element_text(family = "Helvetica"),
    strip.background.x = element_blank(),
    strip.text.x = ggtext::element_markdown(size = 22, face = "bold")
  )
ggsave(
  file.path(out_dir, "de_gene_summary_condition.png"),
  width = 8,
  height = 8
)


# Make volcano plots ------------------------------------------------------

de_genes_unfiltered <- de_comparisons %>% map(
  ~results(dds, contrast = .x, tidy = TRUE) %>%
    left_join(., biomart_table, by = c("row" = "ensembl_gene_id"))
)

volcano_plot(
  x = de_genes_unfiltered$condition_HBV_vs_VEH,
  main_title = "HBV vs VEH",
  xlims = c(-5, 10),
  ylims = c(0, 85)
)
ggsave(
  filename = file.path(out_dir, "volcano_plot_HBV.png"),
  width = 6,
  height = 9
)

volcano_plot(
  x = de_genes_unfiltered$condition_BCG_vs_VEH,
  main_title = "BCG vs VEH",
  xlims = c(-5, 10),
  ylims = c(0, 85)
)
ggsave(
  filename = file.path(out_dir, "volcano_plot_BCG.png"),
  width = 6,
  height = 9
)

volcano_plot(
  x = de_genes_unfiltered$`condition_HBV+BCG_vs_VEH`,
  main_title = "[HBV+BCG] vs VEH",
  xlims = c(-5, 10),
  ylims = c(0, 85)
)
ggsave(
  filename = file.path(out_dir, "volcano_plot_HBV_BCG.png"),
  width = 6,
  height = 9
)


# Upset plot of DE genes --------------------------------------------------

upset_input_list <- de_genes[str_detect(names(de_genes), "_vs_VEH")] %>%
  map(~pull(.x, "gene")) %>%
  set_names(~str_remove_all(.x, "condition_|_vs_VEH"))

svg(
  filename = file.path(out_dir, "upset_plot_condition_de_genes.svg"),
  width  = 10,
  height = 6
)
upset(
  data = fromList(upset_input_list),
  order.by = "freq",
  line.size = 2,
  point.size = 6,
  text.scale = c(3.2, 2.2, 2.2, 1.7, 2.2, 3.2),
  set_size.scale_max = 2500,
  sets.bar.color = c(
    project_colours[["vaccines"]]["HBV+BCG"],
    project_colours[["vaccines"]]["BCG"],
    project_colours[["vaccines"]]["HBV"]
  ),
  mb.ratio = c(0.6, 0.4)
)
dev.off()


# Run pathway enrichment --------------------------------------------------

my_gene_universe <- filter(
  biomart_table,
  ensembl_gene_id %in% rownames(count_matrix)
) %>%
  pull(entrez_gene_id) %>%
  as.character()

# Split DE genes into up and down regulated
de_genes_split <- de_genes %>%
  map(., function(x) {
    list(
      "up" = filter(x, log2FoldChange > 0),
      "down" = filter(x, log2FoldChange < 0)
    ) %>%
      discard(~nrow(.) == 0)
  })

reactomePA_results_raw <- map(de_genes_split, function(x)
  map(x, function(y)
    rpa_wrapper(
      x = as.character(y$entrez_gene_id),
      gene_universe = my_gene_universe,
      species = "human"
    )
  )
)

reactomePA_results_filtered <- map(reactomePA_results_raw, function(x)
  map_dfr(x, function(y)
    filter(y, p_adjust <= 0.05),
    .id = "direction"
  )
)

reactomePA_results_GR <- reactomePA_results_filtered %>% map(
  ~mutate(
    .x,
    rpa_gr = count / as.numeric(str_remove(bg_ratio, "/[0-9]{1,9}"))
  ) %>%
    dplyr::select(id, description, p_adjust, direction, rpa_gr) %>%
    left_join(., reactome_categories_HSA_L1_L2, by = c("id", "description")) %>%
    distinct(id, .keep_all = TRUE)
)

iwalk(
  reactomePA_results_GR,
  ~write_tsv(
    x = .x,
    file = file.path(
      out_dir,
      paste0("split_reactomePA_", str_replace(.y, "\\+", "_"), ".tsv")
    )
  )
)


# Plot pathway results ----------------------------------------------------

# |- Plot immune pathways -------------------------------------------------

pathways_condition_1 <- reactomePA_results_GR %>%
  bind_rows(.id = "condition") %>%
  filter(str_detect(condition, "VEH")) %>%
  mutate(
    numerator = fct_inorder(str_remove_all(condition, "^condition_|_vs_VEH$")),
    description = map_chr(description, ~tRavis::tr_trunc_neatly(.x, l = 50))
  ) %>%
  group_by(description) %>%
  mutate(
    is_different = if_else(length(unique(direction)) == 2, "yes", "no"),
    description  = if_else(
      is_different == "no",
      description,
      paste0("<span style='color:#4582ec;'><b>", description, "</b></span>")
    )
  ) %>%
  ungroup()

# Group any Level 2 terms into one "Other" category if they have < 5 pathways
# assigned to them, and highlight pathways with alternate directions
pathways_condition_2_immune <- pathways_condition_1 %>%
  filter(
    level_1 %in% c("Immune System", "Signal Transduction", "Disease")
  ) %>%
  group_by(level_2) %>%
  mutate(new_level_2 = case_when(n() < 5 ~ "Other", TRUE ~ level_2)) %>%
  ungroup() %>%
  mutate(new_level_2 = str_wrap(new_level_2, width = 25))

pathway_plot2(pathways_condition_2_immune, rows = new_level_2)
ggsave(
  filename = file.path(
    out_dir,
    "reactomePA_results_de_condition_split_immune.png"
  ),
  width = 11,
  height = 19,
  dpi = 100
)


# |- Plot remaining pathways ----------------------------------------------

pathways_condition_2_all <- pathways_condition_1 %>%
  filter(
    !level_1 %in% c("Immune System", "Signal Transduction", "Disease")
  ) %>%
  group_by(level_1) %>%
  mutate(new_level_1 = case_when(n() < 5 ~ "Other", TRUE ~ level_1)) %>%
  ungroup() %>%
  replace_na(list(new_level_1 = "Other")) %>%
  mutate(new_level_1 = str_wrap(new_level_1, width = 25))

pathway_plot2(pathways_condition_2_all, rows = new_level_1)
ggsave(
  filename = file.path(
    out_dir,
    "reactomePA_results_de_condition_split_other.png"
  ),
  width = 12,
  height = 19,
  dpi = 100
)


# Save session information ------------------------------------------------

write_lines(
  capture.output(sessionInfo()),
  file.path(out_dir, "session_info.txt")
)
