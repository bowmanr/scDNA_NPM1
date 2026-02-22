# tabs/single_sample_ui.R
single_sample_ui <- function() {
  
  samples <- sort(unique(as.character(seurat_obj@meta.data[[SAMPLE_COL]])))
  proteins <- rownames(GetAssayData(seurat_obj, assay = "Protein", layer = "data"))
  
  tagList(
    # ------------------------
    # Controls
    # ------------------------
    bs4Card(
      title = "Single sample controls",
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      fluidRow(
        column(
          6,
          selectInput(
            "ss_sample",
            "Sample",
            choices = samples,
            selected = samples[1]
          )
        ),
        column(
          6,
          selectInput(
            "ss_protein",
            "Protein",
            choices = proteins,
            selected = proteins[1]
          )
        )
      )
    ),
    
    # ------------------------
    # Row 1: clonograph + UMAP clone coloring
    # ------------------------
    fluidRow(
      bs4Card(
        title = "Clonograph",
        width = 6,
        status = "warning",
        solidHeader = TRUE,
        plotOutput("ss_clonograph", height = "520px")
      ),
      bs4Card(
        title = "UMAP (clone_code)",
        width = 6,
        status = "info",
        solidHeader = TRUE,
        plotOutput("ss_umap_clone", height = "520px")
      )
    ),
    
    # ------------------------
    # Row 2: cluster UMAP + protein UMAP
    # ------------------------
    fluidRow(
      bs4Card(
        title = "UMAP (clusters)",
        width = 6,
        status = "info",
        solidHeader = TRUE,
        plotOutput("ss_umap_cluster", height = "520px")
      ),
      bs4Card(
        title = "UMAP (protein)",
        width = 6,
        status = "info",
        solidHeader = TRUE,
        plotOutput("ss_umap_protein", height = "520px")
      )
    ),
    
    # ------------------------
    # Row 3: dotplot + barplot (stacked in server)
    # ------------------------
    fluidRow(
      bs4Card(
        title = "Cluster summary (cells per cluster + protein dotplot)",
        width = 12,
        status = "secondary",
        solidHeader = TRUE,
        plotOutput("ss_cluster_summary", height = "950px")
      )
    )
  )
}
