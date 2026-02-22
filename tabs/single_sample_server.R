# tabs/single_sample_server.R
single_sample_server <- function(input, output, session) {
  
  BASE_FONT    <- 20
  UMAP_PT_SIZE <- 1.6
  
  # -------------------------
  # Subset Seurat to sample
  # -------------------------
  ss_obj <- reactive({
    subset_to_sample(seurat_obj, input$ss_sample)
  })
  
  # -------------------------
  # UMAP (clusters)
  # -------------------------
  output$ss_umap_cluster <- renderPlot({
    obj <- ss_obj()
    
    cl_vec <- as.numeric(obj@meta.data[[CLUSTER_COL]])
    cl_levels <- sort(unique(cl_vec))
    cluster_colors <- make_cell_cluster_palette(cl_levels)
    
    DimPlot(
      obj,
      reduction = UMAP_REDUCTION,
      group.by = CLUSTER_COL,
      pt.size = UMAP_PT_SIZE
    ) +
      scale_color_manual(values = cluster_colors) +
      theme(
        text = element_text(size = BASE_FONT),
        axis.text = element_text(size = BASE_FONT - 2),
        legend.text = element_text(size = BASE_FONT - 2),
        legend.title = element_text(size = BASE_FONT)
      )
  })
  
  # -------------------------
  # UMAP (protein)
  # -------------------------
  output$ss_umap_protein <- renderPlot({
    obj <- ss_obj()
    DefaultAssay(obj) <- "Protein"
    
    FeaturePlot(
      obj,
      reduction = UMAP_REDUCTION,
      features = as.character(input$ss_protein)[1],
      pt.size = UMAP_PT_SIZE
    ) +
      scale_color_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 1
      ) +
      scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 1
      ) +
      theme(
        text = element_text(size = BASE_FONT),
        axis.text = element_text(size = BASE_FONT - 2),
        legend.text = element_text(size = BASE_FONT - 2),
        legend.title = element_text(size = BASE_FONT)
      )
  })
  
  # --------------------------------------------------
  # Cluster summary: bar + clustered dotplot (fancy)
  # --------------------------------------------------
  output$ss_cluster_summary <- renderPlot({
    
    obj <- ss_obj()
    DefaultAssay(obj) <- "Protein"
    
    # --------------------
    # Bar plot: cells per cluster
    # --------------------
    cl_vec <- as.character(obj@meta.data[[CLUSTER_COL]])
    df <- as.data.frame(table(cl_vec), stringsAsFactors = FALSE)
    colnames(df) <- c("cluster", "n_cells")
    
    suppressWarnings({
      cl_num <- as.numeric(df$cluster)
    })
    if (!any(is.na(cl_num))) {
      df <- df[order(cl_num), ]
      df$cluster <- factor(df$cluster, levels = df$cluster)
    } else {
      df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster)))
    }
    
    cluster_colors <- make_cell_cluster_palette(levels(df$cluster))
    
    p_bar <- ggplot(df, aes(x = cluster, y = n_cells, fill = cluster)) +
      geom_col() +
      scale_fill_manual(values = cluster_colors) +
      labs(x = "Cluster", y = "Cells") +
      theme_minimal(base_size = BASE_FONT) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.margin = margin(10, 10, 0, 10)
      )
    
    # --------------------
    # DotPlot: proteins × clusters
    # --------------------
    p_dot <- DotPlot(
      obj,
      features = rownames(GetAssayData(obj, assay = "Protein", layer = "data")),
      group.by = CLUSTER_COL
    ) +
      coord_flip() +
      scale_color_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 1
      ) +
      scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 1
      ) +
      labs(x = "Cluster", y = NULL, color = "Avg expr") +
      theme_minimal(base_size = BASE_FONT) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = BASE_FONT - 2),
        plot.margin = margin(0, 10, 10, 10)
      )
    
    cowplot::plot_grid(
      p_bar,
      p_dot,
      ncol = 1,
      rel_heights = c(1, 2.5),
      align = "v",
      axis = "lr"
    )
  })
  
  output$ss_umap_clone <- renderPlot({
    obj <- ss_obj()
    
    # set the metadata column name you want
    CLONE_COL <- "Genotype"
    
    cl_vec <- as.character(obj@meta.data[[CLONE_COL]])
    cl_levels <- sort(unique(cl_vec))
    clone_colors <- make_clone_cluster_palette(cl_levels)
    
    DimPlot(
      obj,
      reduction = UMAP_REDUCTION,
      group.by = CLONE_COL,
      pt.size = UMAP_PT_SIZE
    ) +
      scale_color_manual(values = clone_colors) +
      theme(
        text = element_text(size = BASE_FONT),
        axis.text = element_text(size = BASE_FONT - 2),
        legend.text = element_text(size = BASE_FONT - 2),
        legend.title = element_text(size = BASE_FONT)
      )
  })
  
  # -------------------------
  # Clonograph
  # -------------------------
  output$ss_clonograph <- renderPlot({
    
    
    p <- scDNA::clonograph(sce_list[[input$ss_sample]],complete_only = F)
    
    # Force predictable rendering in Shiny cards
    cowplot::ggdraw(p) +
      theme(
        plot.margin = margin(10, 10, 10, 10),
        text = element_text(size = BASE_FONT)
      )
  }, res = 250)
  
  
}
