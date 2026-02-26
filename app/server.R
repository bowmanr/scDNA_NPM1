# server.R — Cohort (CSV) + Single Sample (Parquet)

server <- function(input, output, session) {

  
  get_clone_fill_map <- function(clones) {
    clones <- as.character(clones)
    
    # use stable ordering for consistent color assignment
    clones_sorted <- sort(unique(clones))
    
    cols <- rep(pals::tol(12), length.out = length(clones_sorted))
    setNames(cols, clones_sorted)
  }
  
  # --------------------
  # Clone tab (CSV)
  # --------------------
  cohort_df <- reactive({
    if (identical(input$cohort_clone_mode, "exact")) {
      df <- pseudo_bulk_exact_clone
      # Standardize to Clone
      if ("Genotype" %in% names(df)) df <- df %>% dplyr::rename(Clone = Genotype)
    } else {
      df <- pseudo_bulk_consolidated_clone
      if ("new_genotype" %in% names(df)) df <- df %>% dplyr::rename(Clone = new_genotype)
    }

    req(all(c("Protein", "Clone", "Expression") %in% names(df)))
    df
  })

  observeEvent(cohort_df(), {
    df <- cohort_df()

    proteins <- sort(unique(df$Protein))
    clones   <- sort(unique(df$Clone))

    default_proteins <- c("CD3","CD117","CD34","CD16","CD14","CD11b")
    default_clones <- c(
      "TET2-only","IDH-only",
      "TET2-NPM1","IDH-NPM1",
      "TET2-NPM1-MAPK","IDH-NPM1-MAPK",
      "WT"
    )

    protein_selected <- intersect(default_proteins, proteins)
    clone_selected   <- intersect(default_clones, clones)

    if (length(protein_selected) == 0) protein_selected <- head(proteins, 6)
    if (length(clone_selected) == 0)   clone_selected   <- head(clones, 6)

    updateSelectizeInput(session, "protein_select",
      choices  = proteins,
      selected = protein_selected
    )

    updateSelectizeInput(session, "clone_select",
      choices  = clones,
      selected = clone_selected
    )
  }, ignoreInit = FALSE)

  cohort_filtered <- reactive({
    df <- cohort_df()

    prot <- input$protein_select
    clon <- input$clone_select

    if (!is.null(prot) && length(prot) > 0) {
      df <- df %>% dplyr::filter(Protein %in% prot)
    }
    if (!is.null(clon) && length(clon) > 0) {
      df <- df %>% dplyr::filter(Clone %in% clon) %>%
        dplyr::mutate(Clone = factor(Clone, levels = clon))
    }

    df
  })

  output$cohort_plot <- renderPlot({
    df <- cohort_filtered()
    req(nrow(df) > 0)

    clones <- as.character(unique(df$Clone))
    clones <- clones[order(clones)]

    fill_map <- setNames(pals::tol(length(clones)), clones)
    
    facet_scales <- if (isTRUE(input$free_y)) "free_y" else "fixed"

    ylim <- NULL
    if (!isTRUE(input$free_y)) {
      y <- df$Expression
      y <- y[is.finite(y)]
      if (length(y) > 0) {
        rng <- range(y)
        pad <- diff(rng) * 0.05
        if (!is.finite(pad) || pad == 0) pad <- 1
        ylim <- c(rng[1] - pad, rng[2] + pad)
      }
    }

    ggplot(df, aes(x = Clone, y = Expression, fill = Clone)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.05, size = 0.4) +
      facet_wrap(~Protein, scales = facet_scales) +
      scale_fill_manual(values = fill_map) +
      { if (!is.null(ylim)) coord_cartesian(ylim = ylim) } +
      labs(x = NULL, y = "Expression") +
      theme_classic(base_size = 12) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey92", color = "grey50"),
        panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.4),
        panel.spacing = grid::unit(0.6, "lines")
      )
  })

  cohort_stats_filtered <- reactive({
    # Stats are consolidated-naming only; for exact mode we still show filtered proteins
    stats <- pseudo_bulk_consolidated_clone_statistics
    req(all(c("Protein","p_val","avg_log2FC","p_val_adj","clone1","clone2") %in% names(stats)))

    prot <- input$protein_select
    clon <- input$clone_select

    if (!is.null(prot) && length(prot) > 0) {
      stats <- stats %>% dplyr::filter(Protein %in% prot)
    }

    # Only meaningful when clone labels are consolidated
    if (identical(input$cohort_clone_mode, "consolidated") && !is.null(clon) && length(clon) > 0) {
      stats <- stats %>% dplyr::filter(clone1 %in% clon & clone2 %in% clon)
    }

    stats %>% dplyr::arrange(Protein, p_val_adj, p_val)
  })

  output$stats_table <- renderTable({
    stats <- cohort_stats_filtered()
    req(nrow(stats) > 0)

    stats %>%
      dplyr::mutate(
        avg_log2FC = round(avg_log2FC, 3),
        p_val      = formatC(p_val, format = "e", digits = 2),
        p_val_adj  = formatC(p_val_adj, format = "e", digits = 2)
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")


  # --------------------
  # Single sample tab (Parquet)
  # --------------------
  observeEvent(TRUE, {
    # patients
    if ("Patient" %in% names(cell_meta)) {
      pats <- sort(unique(as.character(cell_meta$Patient)))
      updateSelectizeInput(session, "ss_patient", choices = pats, selected = pats[1])
    }

    # proteins
    prots <- sort(unique(as.character(expr_long$Protein)))
    updateSelectizeInput(session, "ss_protein", choices = prots, selected = if ("CD3" %in% prots) "CD3" else prots[1])
  }, once = TRUE)

  observeEvent(input$ss_patient, {
    req(input$ss_patient)
    if (!("Stage" %in% names(cell_meta))) {
      updateSelectizeInput(session, "ss_stage", choices = character(0), selected = character(0))
      return()
    }

    stages <- cell_meta %>%
      dplyr::filter(Patient == input$ss_patient) %>%
      dplyr::pull(Stage) %>%
      unique() %>%
      as.character() %>%
      sort()

    updateSelectizeInput(session, "ss_stage", choices = stages, selected = stages)
  }, ignoreInit = FALSE)

  ss_clone_col <- reactive({
    if (identical(input$ss_clone_mode, "exact")) "Genotype" else "new_genotype"
  })

  output$ss_umap_clone <- renderPlot({
    req(cell_meta)
    req(!is.null(input$ss_patient), nzchar(input$ss_patient))
    
    # choose which clone label column to use
    clone_col <- if (identical(input$ss_clone_mode, "exact")) "Genotype" else "new_genotype"
    req(clone_col %in% names(cell_meta))
    
    # filter down to the selected cells
    df <- cell_meta %>%
      dplyr::filter(Patient == input$ss_patient)
    
    if (!is.null(input$ss_stage) && length(input$ss_stage) > 0) {
      df <- df %>% dplyr::filter(Stage %in% input$ss_stage)
    }
    
    req(nrow(df) > 0)
    
    # define Clone for plotting
    df <- df %>% dplyr::mutate(Clone = as.character(.data[[clone_col]]))
    
    # build a stable color map that matches your dotplot bar colors
    fill_map <- get_clone_fill_map(df$Clone)
    
    ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Clone)) +
      geom_point(size = 0.75) +
      scale_color_manual(values = fill_map) +
      theme_classic(base_size = 14) +
      labs(x = NULL, y = NULL, color = "Clone")
  })
  
  output$ss_umap_protein <- renderPlot({
    req(cell_meta, expr_long)
    req(!is.null(input$ss_patient), nzchar(input$ss_patient))
    req(!is.null(input$ss_protein), nzchar(input$ss_protein))
    
    meta <- ss_cells()
    req(nrow(meta) > 0)
    
    # now safe to select()
    meta_small <- meta %>% dplyr::select(Cell, UMAP_1, UMAP_2)
    
    ex <- expr_long %>%
      dplyr::filter(Protein == input$ss_protein, Cell %in% meta_small$Cell) %>%
      dplyr::select(Cell, Expression)
    
    dfp <- meta_small %>% dplyr::left_join(ex, by = "Cell")
    req(nrow(dfp) > 0)
    
    ggplot(dfp, aes(UMAP_1, UMAP_2, color = Expression)) +
      geom_point(size = 0.75) +
      theme_classic(base_size = 14) +
      scale_color_viridis_c()+
      labs(x = NULL, y = NULL, color = input$ss_protein)
  })
  
  ss_cells <- reactive({
    req(cell_meta)
    req(nrow(cell_meta) > 0)
    
    # don’t req input$ss_patient here; handle it in renderPlot
    df <- cell_meta
    
    if (!is.null(input$ss_patient) && nzchar(input$ss_patient)) {
      df <- df %>% dplyr::filter(Patient == input$ss_patient)
    }
    
    stg <- input$ss_stage
    if (!is.null(stg) && length(stg) > 0) {
      df <- df %>% dplyr::filter(Stage %in% stg)
    }
    
    df
  })
  
  output$ss_cluster_summary <- renderPlot({
    meta <- ss_cells()
    req(nrow(meta) > 0)
    
    # ---- Decide grouping (x-axis)
    group_var <- if (identical(input$ss_dot_x, "clone")) {
      # clone mode determines which column
      if (identical(input$ss_clone_mode, "exact")) "Genotype" else "new_genotype"
    } else {
      "Cluster"
    }
    req(group_var %in% names(meta))
    
    meta_g <- meta %>%
      dplyr::mutate(Group = as.character(.data[[group_var]]))
    
    is_clone_x <- identical(input$ss_dot_x, "clone")
    
    # ---- Abundance (bar)
    ab <- meta_g %>%
      dplyr::count(Group, name = "n") %>%
      dplyr::arrange(desc(n))
    
    group_levels <- ab$Group

    
    # ---- Dotplot data: ALL proteins (no user filtering)
    ex <- expr_long %>%
      dplyr::filter(Cell %in% meta_g$Cell) %>%
      dplyr::left_join(meta_g %>% dplyr::select(Cell, Group), by = "Cell") %>%
      dplyr::mutate(Group = factor(Group, levels = group_levels))
    
    dd <- ex %>%
      dplyr::group_by(Protein, Group) %>%
      dplyr::summarise(
        avg_expr = mean(Expression, na.rm = TRUE),
        pct_expr = mean(Expression > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      )
    
    req(nrow(dd) > 0)
    
    # Optional: show proteins ordered (alphabetical); reverse to match common dotplot style
    dd <- dd %>% dplyr::mutate(Protein = factor(Protein, levels = rev(sort(unique(Protein)))))
    
    # ---- Cluster proteins (rows) by their avg expression pattern across Group
    mat_wide <- dd %>%
      dplyr::select(Protein, Group, avg_expr) %>%
      tidyr::pivot_wider(names_from = Group, values_from = avg_expr, values_fill = 0)
    
    mat <- as.matrix(mat_wide[, -1, drop = FALSE])
    rownames(mat) <- mat_wide$Protein
    
    # optional: scale each protein to emphasize patterns rather than absolute levels
    # mat <- t(scale(t(mat)))
    
    hc <- hclust(dist(mat), method = "complete")
    prot_order <- rownames(mat)[hc$order]
    
    dd <- dd %>% dplyr::mutate(Protein = factor(Protein, levels = rev(prot_order)))
    
    if (is_clone_x) {
      clone_fill_map <- get_clone_fill_map(group_levels)
      
      p_bar <- ggplot(ab, aes(x = factor(Group, levels = group_levels), y = n, fill = Group)) +
        geom_col() +
        scale_fill_manual(values = clone_fill_map, guide = "none") +
        theme_classic(base_size = 14) +
        labs(x = NULL, y = "# cells") +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5)
        )
    } else {
      p_bar <- ggplot(ab, aes(x = factor(Group, levels = group_levels), y = n)) +
        geom_col() +
        theme_classic(base_size = 14) +
        labs(x = NULL, y = "# cells") +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5)
        )
    }
    
    
    p_dot <- ggplot(dd, aes(x = Group, y = Protein)) +
      geom_point(aes(size = pct_expr, color = avg_expr), alpha = 0.9) +
      scale_size_continuous(name = "% expressing", range = c(0.5, 8)) +
      scale_color_viridis_c(name = "Mean expr") +
      theme_classic(base_size = 14) +
      labs(x = if (identical(input$ss_dot_x, "clone")) "Clone" else "Cluster", y = NULL) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 5.5, 5.5, 5.5)
      )
    
    (p_bar / p_dot) + patchwork::plot_layout(heights = c(1, 5))
  })  
  
  make_boxplot <- function(df, group_var, y_var,
                           stage_mode = "exclude_CR",
                           palette = NULL) {
    
    # stage filter
    if (identical(stage_mode, "exclude_CR")) {
      df <- df %>% dplyr::filter(Stage != "CR")
    } else if (stage_mode %in% c("Diagnosis","CR","Relapse")) {
      df <- df %>% dplyr::filter(Stage == stage_mode)
    } else {
      # "all" -> no filtering
    }
    
    if (nrow(df) == 0) return(NULL)
    
    # Remove groups with fewer than 2 non-missing observations
    df_valid <- df %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::filter(sum(!is.na(.data[[y_var]])) >= 2) %>%
      dplyr::ungroup()
    
    if (length(unique(df_valid[[group_var]])) < 2) return(NULL)
    
    df_valid[[group_var]] <- droplevels(df_valid[[group_var]])
    
    stats_out <- df_valid %>%
      rstatix::pairwise_wilcox_test(
        formula = stats::as.formula(paste(y_var, "~", group_var)),
        p.adjust.method = "fdr"
      )
    
    max_y <- max(df_valid[[y_var]], na.rm = TRUE)
    if (!is.finite(max_y)) max_y <- 1
    
    # comparisons count drives headroom
    n_comp <- sum(stats_out$p.adj.signif != "ns", na.rm = TRUE)
    ylim_top <- max_y * (1 + max(1, n_comp) * 0.10)
    
    if (is.null(palette)) {
      # fallback if not provided
      pal <- pals::glasbey(length(unique(df_valid[[group_var]])))
      names(pal) <- sort(unique(as.character(df_valid[[group_var]])))
      palette <- pal
    }
    
    p <- ggplot2::ggplot(
      df_valid,
      ggplot2::aes(x = droplevels(.data[[group_var]]), y = .data[[y_var]], fill = .data[[group_var]])
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      ggplot2::geom_jitter(width = 0.25, size = 0.7, height = 0) +
      ggplot2::scale_fill_manual(values = palette) +
      ggpubr::theme_pubr(base_size = 14) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = ggplot2::element_blank()
      ) +
      ggplot2::coord_cartesian(clip = "off", ylim = c(0, ylim_top))
    
    # p-value annotations
    p <- p +
      ggpubr::stat_pvalue_manual(
        data = stats_out,
        y.position = max_y * 1.075,
        step.increase = 0.10,
        hide.ns = TRUE,
        vjust = 0.2,
        label = "p.adj.signif"
      )
    
    list(plot = p, stats = stats_out)
  }
  
  sample_stats_result <- reactive({
    req(sample_stats)
    req(input$ss_group_var, input$ss_metric, input$ss_stage_mode)
    
    grp <- input$ss_group_var
    met <- input$ss_metric
    
    pal <- color_scheme[[grp]]
    if (is.null(pal)) pal <- color_scheme[["Epi_group"]]
    
    make_boxplot(
      df = sample_stats,
      group_var = grp,
      y_var = met,
      stage_mode = input$ss_stage_mode,
      palette = pal
    )
  })
  
  output$sample_stats_plot <- renderPlot({
    res <- sample_stats_result()
    req(!is.null(res))
    res$plot
  })
  
  output$sample_stats_table <- renderTable({
    res <- sample_stats_result()
    req(!is.null(res))
    
    res$stats %>%
      dplyr::mutate(
        p = formatC(p, format = "e", digits = 2),
        p.adj = formatC(p.adj, format = "e", digits = 2)
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
}
