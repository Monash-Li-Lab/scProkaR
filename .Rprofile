## Only try to attach workflowr in interactive sessions, and do not fail the
## session if one of workflowr's dependencies is temporarily unavailable.
if (interactive() && requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  tryCatch(
    suppressPackageStartupMessages(library("workflowr")),
    error = function(e) {
      message("workflowr could not be attached: ", conditionMessage(e))
    }
  )
}
