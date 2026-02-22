# tabs/cohort_ui.R
cohort_ui <- function() {
  
  tagList(
    bs4Card(
      title = "Cohort controls",
      width = 12,
      status = "primary",
      solidHeader = TRUE,
      fluidRow(
        column(
          3,
          selectizeInput("co_features", "Features",
                         choices = sort(unique(summarized_out$Feature)),
                         selected = c("CD14","CD117","CD34"),
                         multiple = TRUE)
        ),
        column(
          2,
          selectizeInput("co_stages", "Stages",
                         choices = sort(unique(summarized_out$Stage_code)),
                         selected = "D",
                         multiple = TRUE)
        ),
        column(
          3,
          selectizeInput("co_genes_include", "Genes include",
                         choices = c("WT","NRAS","FLT3","TET2","IDH2","NPM1","KRAS","DNMT3A","PTPN11","IDH1"),
                         selected = c("WT","NRAS","FLT3","TET2","IDH2","NPM1"),
                         multiple = TRUE)
        ),
        column(
          2,
          selectizeInput("co_genes_exclude", "Exclude",
                         choices = c("IDH1","DNMT3A","PTPN11","KRAS"),
                         selected = c("IDH1","DNMT3A","PTPN11"),
                         multiple = TRUE)
        ),
        column(
          2,
          numericInput("co_size_cutoff", "Min % cells", value = 1, min = 0, step = 0.5)
        )
      ),
      fluidRow(
        column(
          3,
          selectizeInput("co_must_include", "Must include in Genotype",
                         choices = c("NPM1","FLT3","NRAS","TET2","IDH2","KRAS"),
                         selected = "NPM1",
                         multiple = TRUE)
        ),
        column(
          3,
          sliderInput("co_expr_thresh", "Expression threshold (Value >)",
                      min = 0, max = 5, value = 1, step = 0.25)
        ),
        column(
          3,
          checkboxInput("co_show_labels", "Show patient labels (slow / busy)", value = FALSE)
        )
      )
    ),
    
    bs4Card(
      title = "Cohort protein expression by genotype",
      width = 12,
      status = "secondary",
      solidHeader = TRUE,
      plotOutput("cohort_plot", height = "850px")
    ),
  
  bs4Card(
    title = "Violins: protein distribution across cohort",
    width = 12,
    status = "info",
    solidHeader = TRUE,
    
    fluidRow(
      column(
        3,
        radioButtons(
          "vio_group",
          "Group by",
          choices = c("Genotype" = "Genotype", "Patient_code" = "Patient_code"),
          selected = "Genotype",
          inline = TRUE
        )
      )
    ),
    
    plotOutput("vio_plot", height = "800px")
  )
  )
}
