# server.R
source("tabs/single_sample_server.R")
source("tabs/cohort_server.R")

server <- function(input, output, session) {
  single_sample_server(input, output, session)
  cohort_server(input, output, session)
}