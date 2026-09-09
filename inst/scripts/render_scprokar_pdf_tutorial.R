if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install 'rmarkdown' before rendering the PDF tutorial.", call. = FALSE)
}

if (!requireNamespace("tinytex", quietly = TRUE)) {
  message("If PDF rendering fails because LaTeX is missing, run: install.packages('tinytex'); tinytex::install_tinytex()")
}

rmarkdown::render(
  input = file.path("vignettes", "scprokar-pdf-tutorial.Rmd"),
  output_format = "pdf_document",
  output_file = "scprokar-comprehensive-tutorial.pdf",
  output_dir = file.path("vignettes"),
  clean = TRUE
)
