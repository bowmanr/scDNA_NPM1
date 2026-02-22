# tools/aboutTab.R
# About page UI
tagList(
  div(class = "about",
      includeMarkdown("tools/about.Rmd")
  )
)