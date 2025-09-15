#' Remove columns with zero variance
#'
#' @param df
#'
#' @return
#' @export
#'
remove_zero_var_cols <- function(df) {
  cols_to_remove <- caret::nearZeroVar(as.data.frame(df), names = TRUE)
  df_removed <- df %>% dplyr::select(!one_of(cols_to_remove))
  return(df_removed)
}


#' Perform PCA
#'
#' @param object Object of class DESeqTransform, e.g. output of DESeq2::vst(dds)
#' @param ntop Number of genes to use in the PCA
#' @param nrank Number of principle components to calculate
#'
#' @return Data frame of PCA results and "intgroup" information
#' @export
#'
run_pca <- function(object, intgroup, ntop = 500, nrank = 2) {

  if (class(object) == "DESeqTransform") {
    message("DESeq object detected. Performing PCA...")

    if (!all(intgroup %in% names(colData(object)))) {
      stop("The argument 'intgroup' should specify columns of colData(dds)")
    }

    # Calculate the variance for each gene
    rv <- rowVars(assay(object))

    # Select the "ntop" genes by variance
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

    # Perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select, ]), rank. = nrank)

    # Contribution to the total variance for each component
    percent_var <- (pca$sdev ^ 2 / sum(pca$sdev ^ 2))[1:nrank]
    names(percent_var) <- paste0("PC", 1:nrank)

    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

    # Assemble the data for output
    d <- cbind(
      data.frame(library_name = colnames(object), intgroup.df),
      pca$x
    )

    attr(d, "percent_var") <- percent_var
    message("Done.\n")
    return(d)

  } else {
    # This section operates on the assumption that we aren't handling expression
    # data, instead looking at e.g. CIBERSORT results
    message("Data frame detected. Performing PCA...")

    pca <- object %>%
      as.data.frame() %>%
      remove_zero_var_cols() %>%
      prcomp(x = ., rank. = nrank, scale. = TRUE)

    percent_var <- (pca$sdev ^ 2 / sum(pca$sdev ^ 2))[1:nrank]
    names(percent_var) <- paste0("PC", 1:nrank)

    if (tail(cumsum(percent_var), 1) >= 0.9) {
      pc_percent_var_90 <-
        names(percent_var[seq(1, which(cumsum(percent_var) > 0.9)[1])])
    } else {
      pc_percent_var_90 <- names(percent_var)
    }

    x <- pca$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column("library_name") %>%
      select(library_name, all_of(pc_percent_var_90)) %>%
      rename_with(~gsub("PC", "CS", .x))

    attr(x, "percent_var") <- percent_var
    message("Done.\n")
    return(x)
  }
}


#' Find the contributions of metadata variables to principle components
#'
#' @param deseq_object DESeq dataset, such that `assay(deseq_object)` returns
#'   transformed counts
#' @param vars_of_interest Character vector; variables of interest to test.
#' @param n_genes Number of most-variable genes to use in the PCA
#' @param n_PCs Number of principle components to check
#'
#' @return Data frame containing significant results (p.value < 0.05)
#' @export
#'
model_var_contributions <- function(
  deseq_object,
  vars_of_interest,
  cell_props = NULL,
  n_genes = 500,
  n_PCs = 10
) {

  expr_data <- assay(deseq_object)

  # Remove low variance genes
  gene_vars <- expr_data %>% apply(1, function(x) {
    var(x)
  })
  remove <- c(which(gene_vars < 0.01))

  if (!is_empty(remove)) {
    expr_data <- expr_data[-remove, ]
  }

  pca_result <- run_pca(
    object   = deseq_object,
    intgroup = vars_of_interest,
    ntop     = n_genes,
    nrank    = n_PCs
  )

  if (!is.null(cell_props)) {
    cell_props_pca <- cell_props %>%
      run_pca(object = ., nrank = 10)

    pca_result <- left_join(
      pca_result,
      cell_props_pca,
      by = "library_name"
    )

    vars_of_interest <- c(vars_of_interest, str_subset(colnames(pca_result), "^CS"))
  }

  message("Run linear models...")
  pca_contribs_raw <- str_subset(colnames(pca_result), "^PC[0-9]{1,2}$") %>% map(
    ~lm(
      formula = paste0(.x, " ~ ", paste0(vars_of_interest, collapse = " + ")),
      data = pca_result
    )
  ) %>% set_names(paste0("PC", 1:n_PCs))

  pca_contribs <- str_subset(colnames(pca_result), "^PC[0-9]{1,2}$") %>% map(
    ~lm(
      formula = paste0(.x, " ~ ", paste0(vars_of_interest, collapse = " + ")),
      data = pca_result
    ) %>% broom::tidy()
  ) %>% set_names(paste0("PC", 1:n_PCs))

  pca_contribs_tidy <- pca_contribs %>%
    bind_rows(.id = "PC") %>%
    filter(
      term != "(Intercept)",
      p.value < 0.05
    ) %>%
    mutate(
      log_p_value = -log10(p.value),
      PC = factor(PC, levels = names(pca_contribs))
    ) %>%
    select(PC, term, estimate, statistic, p.value, log_p_value)

  message("Done.\n")
  return(list(
    "results" = pca_contribs_tidy,
    "percent_vars" = attr(pca_result, "percent_var"),
    "results_raw" = pca_contribs_raw
  ))
}
