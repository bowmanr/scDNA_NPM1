# tools/homeTab.R
# Home page content
tagList(
  div(id = "home",
      h2("Welcome to the Bowman Lab Datasets page!"),
      p("This web portal was created with the intent of making gene expression datasets more accessible to members of the Bowman lab. Here you will find tools to visualize gene queries across a breadth of datasets generated in the Bowman lab and by colleagues. For a comprehensive description of how each dataset was generated, check out the About tab. To get started with looking at a set of genes you are interested, take a look below."),
      h3("Getting Started:"),
      h4("Have a specific gene in mind?"),
      tags$ol(
        tags$li('Select the "Query" tab for individual gene analyses'),
        tags$li('Choose a dataset'),
        tags$li('Enter a Gene Symbol'),
        tags$li('Select one of the available plots')
      ),
      h4("Have a particular dataset that you want to look at?"),
      tags$ol(
        tags$li('Select the Explore tab'),
        tags$li('Choose a dataset'),
        tags$li('Select two sample sets to compare'),
        tags$li('Set fold change and FDR cutoffs')
      ),
      h4("To be completed:"),
      tags$ol(
        tags$li('Complete dataset description and methodology'),
        tags$li('Password protection and hosting')
      )
  )
)