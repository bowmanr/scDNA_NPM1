# ui.R
source("tabs/single_sample_ui.R")
source("tabs/cohort_ui.R")

ui <- bs4DashPage(
  title = "scDNA Explorer",
  header = bs4DashNavbar(),
  sidebar = bs4DashSidebar(
    skin = "light",
    status = "primary",
    title = "Explorer",
    bs4SidebarMenu(
      bs4SidebarMenuItem("Single sample", tabName = "single_sample", icon = icon("flask")),
      bs4SidebarMenuItem("Cohort", tabName = "cohort", icon = icon("users"))
    )
  ),
  body = bs4DashBody(
    bs4TabItems(
      bs4TabItem(tabName = "single_sample", single_sample_ui()),
      bs4TabItem(tabName = "cohort", cohort_ui())
    )
  )
)
