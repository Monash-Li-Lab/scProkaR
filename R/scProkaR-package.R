#' SCProkaR: Microbial Single-Cell RNA-Seq Analysis for Bacterial Systems
#'
#' `SCProkaR` provides a `SingleCellExperiment`-centered toolkit for microbial
#' single-cell analysis in R. The package combines bacterial QC and integration
#' workflows with pseudobulk differential expression utilities and Time Aware
#' Trajectory Analysis (TATA) for time-informed trajectory abstraction.
#'
#' @section Example dataset:
#' The package includes `sce1`, a small `SingleCellExperiment` object derived
#' from a three-treatment bacterial single-cell experiment. The object contains
#' raw counts together with `sample`, `clusters`, `treatments`, and
#' `timepoints` metadata, and is used throughout the package workflow vignette.
#'
#' @keywords internal
"_PACKAGE"

#' Example bacterial single-cell dataset
#'
#' `sce1` is a packaged `SingleCellExperiment` example dataset for demonstrating
#' the core `SCProkaR` workflow. It contains one observation per bacterium and
#' includes three treatment groups (`control`, `PMB0.5`, and `PMB2`) measured
#' across multiple time points.
#'
#' @format A `SingleCellExperiment` with a counts assay and cell metadata
#' columns:
#' \describe{
#'   \item{sample}{original sample identifier}
#'   \item{clusters}{pre-existing cluster label}
#'   \item{treatments}{treatment group}
#'   \item{timepoints}{experimental sampling time}
#' }
#' @source Internal example dataset distributed with `SCProkaR`.
"sce1"
