# tools/queryTab.R
library(shiny)
library(bs4Dash)

tagList(
  fluidRow(
    bs4Card(
      title       = "Controls",
      width       = 4,
      status      = "primary",
      solidHeader = TRUE,
      
      # 1) Dataset selector
      selectInput(
        inputId  = "dataset",
        label    = "Select Dataset",
        choices  = names(dataset_paths),
        selected = names(dataset_paths)[1]
      ),
      
      # 2) Gene selector
      selectizeInput(
        inputId   = "genes",
        label     = "Select Genes",
        choices   = NULL,
        multiple  = TRUE
      ),
      
      # 3) Plot type (populated server-side)
      uiOutput("plot_type_ui"),
      
      # 4) Always-visible metadata selectors (choices updated server-side)
      selectInput(
        inputId  = "identity",
        label    = "Set Identity (Idents)",
        choices  = NULL
      ),
      selectInput(
        inputId  = "split_by",
        label    = "Split By",
        choices  = NULL
      ),
      selectInput(
        inputId  = "group_by",
        label    = "Group By",
        choices  = NULL
      )
    ),
    
    bs4Card(
      title       = "Plot",
      width       = 8,
      status      = "primary",
      solidHeader = TRUE,
      plotOutput("gene_plot", height = "700px")
    )
  )
)
