## ---------------------------------------------------------------------------
## Provenance of data/sce1.rda
##
## Bioconductor requires that every object shipped in data/ can be traced back
## to its source. This script documents how `sce1` was produced.
##
## ---------------------------------------------------------------------------
## 1. SOURCE EXPERIMENT
## ---------------------------------------------------------------------------
##
## `sce1` is a downsampled subset of a bacterial single-cell RNA-seq time
## course in which a bacterial culture was profiled under polymyxin B (PMB)
## exposure at two concentrations against an untreated control.
##
## The design that is preserved in the packaged object is:
##
##   treatment  timepoint (h)   cells
##   control    0, 1, 4, 7      500 per timepoint  (2000 total)
##   PMB0.5     1, 4, 7         500 per timepoint  (1500 total)
##   PMB2       1, 4, 7         500 per timepoint  (1500 total)
##
## giving ten samples of exactly 500 cells each, i.e. 5000 cells in total.
##
## The following fields MUST be completed by the package authors before
## submission, because they cannot be recovered from the object itself:
##
##   ORGANISM AND STRAIN : <fill in, e.g. Acinetobacter baumannii strain ...>
##   REFERENCE GENOME    : <fill in accession, e.g. GCF_...>
##   ASSAY / PROTOCOL    : <fill in, e.g. BacDrop / MicroSPLiT / PETRI-seq>
##   RAW DATA ACCESSION  : <fill in, e.g. GEO GSE...... or SRA PRJNA......>
##   PROCESSING PIPELINE : <fill in aligner / counting tool and versions>
##   PUBLICATION         : <fill in DOI once available>
##
## ---------------------------------------------------------------------------
## 2. WHAT THE PACKAGED OBJECT CONTAINS
## ---------------------------------------------------------------------------
##
## class            : SingleCellExperiment
## dim              : 3722 features x 5000 cells
## assay            : "counts", a sparse dgCMatrix of raw integer UMI counts
##                    (97% zeros; per-cell library sizes 28-2828;
##                    20-509 features detected per cell)
## rownames         : 2868 locus tags of the form "ABUW-RSnnnnn" plus 854
##                    gene symbols such as "dnaA", "gyrB", "recF"
## colnames         : "<sample>_<cell barcode>", unique across samples
## colData          : sample     - factor, the ten sample identifiers
##                    clusters   - character, published cluster labels A-F
##                    treatments - character, "control", "PMB0.5", "PMB2"
##                    timepoints - character, "0", "1", "4", "7" (hours)
## reducedDims      : none; the vignettes compute these from the counts
##
## ---------------------------------------------------------------------------
## 3. HOW THE OBJECT WAS BUILT
## ---------------------------------------------------------------------------
##
## The code below reproduces `sce1` from the per-sample count matrices of the
## source experiment. Replace `raw_dir` with the directory holding those
## matrices. It is not executed when the package is built or checked.

if (FALSE) {

    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(Matrix)

    raw_dir <- "<path to the per-sample count matrices>"

    samples <- c(
        "control_0", "control_1", "control_4", "control_7",
        "PMB0.5_1", "PMB0.5_4", "PMB0.5_7",
        "PMB2_1", "PMB2_4", "PMB2_7"
    )

    ## Read one count matrix per sample and prefix the barcodes with the
    ## sample name so that cell identifiers stay unique after cbind().
    mats <- lapply(samples, function(s) {
        m <- readMM(file.path(raw_dir, s, "matrix.mtx.gz"))
        rownames(m) <- readLines(file.path(raw_dir, s, "features.tsv.gz"))
        colnames(m) <- paste(s, readLines(
            file.path(raw_dir, s, "barcodes.tsv.gz")
        ), sep = "_")
        as(m, "CsparseMatrix")
    })
    names(mats) <- samples

    ## Keep the genes shared by every sample, then downsample each sample to
    ## exactly 500 cells to keep the packaged object under the Bioconductor
    ## 5 MB data limit.
    common_genes <- Reduce(intersect, lapply(mats, rownames))

    set.seed(2026)
    mats <- lapply(mats, function(m) {
        m <- m[common_genes, , drop = FALSE]
        m[, sample(seq_len(ncol(m)), 500), drop = FALSE]
    })

    counts <- do.call(cbind, mats)

    sample_id <- sub("_[ACGT]+$", "", colnames(counts))
    sce1 <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            sample     = factor(sample_id, levels = samples),
            clusters   = NA_character_,  # published cluster labels A-F
            treatments = sub("_[0-9.]+$", "", sample_id),
            timepoints = sub("^.*_", "", sample_id),
            row.names  = colnames(counts)
        )
    )

    ## `clusters` holds the cluster labels reported for the source experiment;
    ## transfer them here from the published cell annotation.

    ## Compress with xz, which is what the shipped file uses.
    save(sce1, file = "data/sce1.rda", compress = "xz")
    tools::checkRdaFiles("data")
}
