library(bs4Dash)
library(shiny)

ui <- bs4DashPage(
  
  title = "scDNA Dashboard",
  
  header = bs4DashNavbar(),
  
  sidebar = bs4DashSidebar(
    bs4SidebarMenu(
      bs4SidebarMenuItem(
        "Home",
        tabName = "home",
        icon = icon("home")
      ),
      
      bs4SidebarMenuItem(
        "Sample stats",
        tabName = "sample_stats",
        icon = icon("chart-line")
      ),
      
      bs4SidebarMenuItem(
        "Clone",
        tabName = "clone",
        icon = icon("chart-bar")
      ),
      
      bs4SidebarMenuItem(
        "Single Sample",
        tabName = "single_sample",
        icon = icon("project-diagram")
      )
      
    )
  ),
  
  body = bs4DashBody(
    
  
    bs4TabItems(
      # -----------------------
      # Home Tab
      # -----------------------
      bs4TabItem(
        tabName = "home",
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Welcome",
            status = "primary",
            solidHeader = TRUE,
            
            HTML("
        <h3>scDNA NPM1 Interactive Analysis Portal</h3>
        <p>
        This dashboard provides interactive visualization and statistical exploration 
        of single-cell DNA sequencing data from the NPM1 clonal evolution study.
        </p>

        <h4>Available Analysis Tabs</h4>

        <ul>
          <li><b>Cohort:</b> Protein expression comparisons across clones with statistical testing.</li>
          <li><b>Single Sample:</b> UMAP visualization, clone distributions, and cluster summaries for individual patients.</li>
          <li><b>Sample Statistics:</b> Group-level comparisons of diversity, mutation burden, and clonal structure.</li>
        </ul>

        <h4>Resources</h4>
        <p>
          Raw data and source code are available at:<br>
          <a href='https://github.com/bowmanr/scDNA_NPM1' target='_blank'>
            https://github.com/bowmanr/scDNA_NPM1
          </a>
        </p>

        <p>
          Preprint available at:<br>
          <a href='https://www.biorxiv.org/content/10.1101/2024.11.11.623033v1' target='_blank'>
            https://www.biorxiv.org/content/10.1101/2024.11.11.623033v1
          </a>
        </p>

        <h4>Contact</h4>
        <p>
          For questions about the dataset or analysis:
        </p>

        <ul>
          <li>Linde Miles: <span id='email_linde'></span></li>
          <li>Robert Bowman: <span id='email_robert'></span></li>
        </ul>

        <script>
          document.getElementById('email_linde').innerHTML =
            'linde' + '.' + 'miles' + '@' + 'cchmc' + '.' + 'org';

          document.getElementById('email_robert').innerHTML =
            'robert' + '.' + 'bowman' + '@' + 'pennmedicine' + '.' + 'upenn' + '.' + 'edu';
        </script>
      ")
          )
        )
      ) ,
      # -----------------------
      # Sample Stats Tab
      # -----------------------
      bs4TabItem(
        tabName = "sample_stats",
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Plot controls",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(4,
                     selectInput("ss_stage_mode", "Stage filter",
                                 choices = c("Diagnosis only"="Diagnosis",
                                             "CR only"="CR",
                                             "Relapse only"="Relapse",
                                             "All (exclude CR)"="exclude_CR",
                                             "All (include CR)"="all"),
                                 selected = "exclude_CR")
              ),
              column(4,
                     selectInput("ss_group_var", "Group by",
                                 choices = c("Epi_group","Sig_group","Stage"),
                                 selected = "Epi_group")
              ),
              column(4,
                     selectInput("ss_metric", "Metric",
                                 choices = sample_metrics,
                                 selected = "Shannon")
              )
            )
          )
        ) ,
        
        bs4Card(
          width = 12,
          title = "Sample statistics",
          status = "info",
          solidHeader = TRUE,
          
          fluidRow(
            column(
              width = 5,
              plotOutput("sample_stats_plot", height = "450px")
            ),
            column(
              width = 7,
              tableOutput("sample_stats_table")
            )
          )
        )
      ) ,
      # -----------------------
      # Clone Tab
      # -----------------------
      bs4TabItem(
        tabName = "clone",
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Clone Specific Immunophenotypes",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                4,
                radioButtons(
                  "cohort_clone_mode",
                  "Clone labels",
                  choices = c(
                    "Consolidated" = "consolidated",
                    "Exact Genes" = "exact"
                  ),
                  selected = "consolidated"
                )
              ),
              column(
                4,
                selectizeInput(
                  "protein_select",
                  "Proteins",
                  choices = NULL,
                  multiple = TRUE
                )
              ),
              column(
                4,
                selectizeInput(
                  "clone_select",
                  "Clones",
                  choices = NULL,
                  multiple = TRUE
                )
              )
            ),
            
            checkboxInput(
              "free_y",
              "Free y-scale per protein",
              value = FALSE
            )
          )
        ),
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Cohort Expression",
            status = "info",
            solidHeader = TRUE,
            plotOutput("cohort_plot", height = "650px")
          )
        ),
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Pairwise Statistics",
            status = "secondary",
            solidHeader = TRUE,
            tableOutput("stats_table")
          )
        )
        
      ),
      
      # -----------------------
      # Single Sample Tab
      # -----------------------
      bs4TabItem(
        tabName = "single_sample",
        
        fluidRow(
          bs4Card(
            width = 12,
            title = "Single Sample Controls",
            status = "primary",
            solidHeader = TRUE,
            
            fluidRow(
              column(
                3,
                selectInput(
                  "ss_patient",
                  "Patient",
                  choices = NULL
                )
              ),
              column(
                3,
                selectizeInput(
                  "ss_stage",
                  "Stage",
                  choices = NULL,
                  multiple = TRUE
                )
              ),
              column(
                3,
                radioButtons(
                  "ss_clone_mode",
                  "Clone labels",
                  choices = c(
                    "Consolidated (new_genotype)" = "consolidated",
                    "Exact (Genotype)" = "exact"
                  ),
                  selected = "consolidated"
                )
              ),
              column(
                3,
              radioButtons(
                "ss_dot_x",
                "DotPlot x-axis",
                choices = c("Clone" = "clone","Cluster" = "cluster" ),
                selected = "clone",
                inline = TRUE
              )
              ),
              column(
                3,
                selectInput(
                  "ss_protein",
                  "Protein",
                  choices = NULL
                )
              )
            )
          )
        ),
        
        fluidRow(
          bs4Card(
            width = 6,
            title = "UMAP colored by Clone",
            status = "info",
            solidHeader = TRUE,
            plotOutput("ss_umap_clone", height = "500px")
          ),
          
          bs4Card(
            width = 6,
            title = "UMAP colored by Protein Expression",
            status = "info",
            solidHeader = TRUE,
            plotOutput("ss_umap_protein", height = "500px")
          ),
          bs4Card(
            width = 12,
            title = "Cluster summary",
            status = "info",
            solidHeader = TRUE,
              plotOutput("ss_cluster_summary", height = "850px")
          )
        )
        
      )
      
    )
  ),
  
  controlbar = NULL,
  footer = NULL
)