# tabs/cohort_server.R
cohort_server <- function(input, output, session) {
  
  output$cohort_plot <- renderPlot({
    
    df <- summarized_out
    
    # Apply filters
    if (length(input$co_features) > 0) {
      df <- df[df$Feature %in% input$co_features, , drop = FALSE]
    }
    if (length(input$co_stages) > 0) {
      df <- df[df$Stage_code %in% input$co_stages, , drop = FALSE]
    }
    
    # Filter by genes include/exclude on BOTH Genotype string and Gene column
    inc <- input$co_genes_include
    exc <- input$co_genes_exclude
    
    if (length(inc) > 0) {
      inc_pat <- paste(inc, collapse = "|")
      df <- df[grepl(inc_pat, df$Genotype) & grepl(inc_pat, df$Gene), , drop = FALSE]
    }
    if (length(exc) > 0) {
      exc_pat <- paste(exc, collapse = "|")
      df <- df[!grepl(exc_pat, df$Genotype) & !grepl(exc_pat, df$Gene), , drop = FALSE]
    }
    
    # Must include logic (keep WT always)
    must <- input$co_must_include
    if (length(must) > 0) {
      must_pat <- paste(must, collapse = "|")
      keep <- grepl("WT", df$Genotype) | grepl(must_pat, df$Genotype)
      df <- df[keep, , drop = FALSE]
    }
    
    # recompute pct_express if you want threshold to be dynamic:
    # NOTE: currently pct_express was computed at Value>1 during preprocessing.
    # If you want this slider to truly change pct_express, we’d need to compute from `out`.
    # For now: keep med_express filter only, and use precomputed pct_express.
    df <- df[df$size_scale > input$co_size_cutoff, , drop = FALSE]
    
    # Factor ordering for nicer facets
    df$Stage_code <- factor(df$Stage_code, levels = c("D","CR","R"))
    if (length(inc) > 0) {
      df$Gene <- factor(df$Gene, levels = c("WT", setdiff(inc, "WT")))
    }
    
    p <- ggplot(df, aes(
      x = Genotype,
      y = pct_express,
      color = Patient,
      group = Patient,
      size = size_scale
    )) +
      geom_jitter(width = 0.05, alpha = 0.2) +
      geom_line(linewidth = 0.5) +
      theme_classic(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_grid(Feature ~ Stage_code)
    
    if (isTRUE(input$co_show_labels)) {
      p <- p + geom_label(aes(label = Patient), size = 2)
    }
    
    p
  }, res = 200)
  
  
  output$vio_plot <- renderPlot({
    
    grp <- input$vio_group
    
    # Cohort controls (shared)
    feats  <- input$co_features
    stages <- input$co_stages
    inc    <- input$co_genes_include
    exc    <- input$co_genes_exclude
    must   <- input$co_must_include
    
    md <- seurat_obj@meta.data
    md$cell_id <- rownames(md)
    
    # Stage filter
    if (!is.null(stages) && length(stages) > 0) {
      md <- md[md$Stage_code %in% stages, , drop = FALSE]
    }
    
    # Genotype include/exclude based on string matching (mirrors summarized_out logic)
    if (!is.null(inc) && length(inc) > 0) {
      inc_pat <- paste(inc, collapse = "|")
      md <- md[grepl(inc_pat, md$Genotype), , drop = FALSE]
    }
    if (!is.null(exc) && length(exc) > 0) {
      exc_pat <- paste(exc, collapse = "|")
      md <- md[!grepl(exc_pat, md$Genotype), , drop = FALSE]
    }
    
    # Must include (always keep WT)
    if (!is.null(must) && length(must) > 0) {
      must_pat <- paste(must, collapse = "|")
      keep <- grepl("WT", md$Genotype) | grepl(must_pat, md$Genotype)
      md <- md[keep, , drop = FALSE]
    }
    
    # Safety: if no features selected, nothing to plot
    feats <- as.character(feats)
    if (length(feats) == 0) return(NULL)
    
    DefaultAssay(seurat_obj) <- "Protein"
    
    # Pull protein expression for selected cells+features only
    mat <- GetAssayData(seurat_obj, assay = "Protein", layer = "data")[feats, md$cell_id, drop = FALSE]
    expr_df <- as.data.frame(t(as.matrix(mat)))
    expr_df$cell_id <- rownames(expr_df)
    
    df <- dplyr::left_join(md[, c("cell_id", "Patient", "Patient_code", "Genotype", "Stage_code")], expr_df, by = "cell_id")
    
    df_long <- tidyr::pivot_longer(
      df,
      cols = dplyr::all_of(feats),
      names_to = "Feature",
      values_to = "Value"
    )
    
    df_long[[grp]] <- factor(df_long[[grp]])
    
    ggplot(df_long, aes(x = .data[[grp]], y = Value)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.12, outlier.size = 0.2) +
      facet_wrap(~Feature, scales = "free_y", ncol = min(4, length(unique(df_long$Feature)))) +
      theme_classic(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = grp, y = "Protein expression (data slot)")
  }, res = 200)
  
  
  
}
