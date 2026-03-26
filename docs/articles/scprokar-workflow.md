# SCProkaR Core Workflow

## Overview

This vignette shows three practical `SCProkaR` workflows:

- a three-batch preprocessing and integration workflow starting from
  Seurat `.rds` objects
- a pseudobulk differential-expression workflow on the same example
  dataset after integration
- a time-aware trajectory analysis workflow on the same example dataset

`SCProkaR` keeps the main analysis object as a `SingleCellExperiment`,
while still accepting Seurat input for convenient conversion.

## Load packages

``` r
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: 'sp'
#> The following object is masked from 'package:IRanges':
#> 
#>     %over%
#> 
#> Attaching package: 'SeuratObject'
#> The following object is masked from 'package:SummarizedExperiment':
#> 
#>     Assays
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     intersect
#> The following object is masked from 'package:Seqinfo':
#> 
#>     intersect
#> The following object is masked from 'package:IRanges':
#> 
#>     intersect
#> The following object is masked from 'package:S4Vectors':
#> 
#>     intersect
#> The following object is masked from 'package:BiocGenerics':
#> 
#>     intersect
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
#> 
#> Attaching package: 'Seurat'
#> The following object is masked from 'package:SummarizedExperiment':
#> 
#>     Assays
```

## Load Seurat `.rds` objects

The package includes a small example object, `sce1`, which we load
directly from `data/` in the same style as other Bioconductor workflow
packages.

``` r
data("sce1", package = "SCProkaR")
sce <- sce1

seu <- CreateSeuratObject(counts = assay(sce, "counts"), assay = "RNA", meta.data = as.data.frame(colData(sce)))

seu1 <- seu[, seu$treatments == "PMB0.5"]
seu2 <- seu[, seu$treatments == "PMB2"]
seu3 <- seu[, seu$treatments == "control"]
```

For clarity, add sample origin columns if they do not already exist:

``` r

RunDimRed <- function(sce){
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(sce)
  sce <- ScaleData(sce, features = all.genes)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:10)
  sce <- RunUMAP(sce, dims = 1:10)
  return(sce)
}

seu1 <- RunDimRed(seu1)
#> Normalizing layer: counts
#> Finding variable features for layer counts
#> Centering and scaling data matrix
#> PC_ 1 
#> Positive:  basJ, ABUW-RS14490, ABUW-RS02205, ABUW-RS09975, tcuB, bauA, rplD, ABUW-RS18285, rpsJ, ABUW-RS14475 
#>     hutU, ABUW-RS18155, ABUW-RS02150, ABUW-RS13085, ABUW-RS00175, ABUW-RS12220, ABUW-RS11275, ABUW-RS06670, ABUW-RS18015, ABUW-RS15565 
#>     ABUW-RS14160, coaBC, ABUW-RS15985, ABUW-RS00195, ABUW-RS09630, ABUW-RS00740, maf, ABUW-RS06595, ABUW-RS13780, cysE 
#> Negative:  katE, ABUW-RS17440, ABUW-RS15165, ABUW-RS14025, ABUW-RS14775, ABUW-RS12405, ABUW-RS11840, ABUW-RS12725, ABUW-RS11875, ABUW-RS11850 
#>     ABUW-RS18980, carO, ABUW-RS09995, ABUW-RS07145, ABUW-RS13220, ABUW-RS07170, ABUW-RS07960, ABUW-RS15330, ABUW-RS07500, ABUW-RS05995 
#>     ABUW-RS03195, ABUW-RS15280, ABUW-RS20195, ABUW-RS14240, ABUW-RS13010, acnA, ABUW-RS00260, ABUW-RS05940, ABUW-RS20155, ABUW-RS11870 
#> PC_ 2 
#> Positive:  ABUW-RS07500, ABUW-RS00825, ABUW-RS15335, ABUW-RS20195, ABUW-RS15330, ABUW-RS03380, ABUW-RS18450, ABUW-RS15165, ABUW-RS03195, ABUW-RS02175 
#>     ABUW-RS16155, ABUW-RS14490, acnA, ABUW-RS07170, ABUW-RS11840, ABUW-RS09210, carO, clpB, ABUW-RS08700, ABUW-RS13120 
#>     ABUW-RS12725, ABUW-RS16690, ABUW-RS15280, ABUW-RS18740, ABUW-RS18980, ABUW-RS11850, ABUW-RS12805, ABUW-RS20475, ABUW-RS03180, ABUW-RS20155 
#> Negative:  rplB, rpsC, secY, rplV, rpsB, rplD, rpsD, rplE, rplF, rpoA 
#>     rpsK, rplX, rpsM, rplN, rpsJ, rpmC, rpsQ, rplM, rpsN, rpsS 
#>     rplW, rplP, rpsE, ABUW-RS04480, rplQ, atpE, ABUW-RS16090, rpsH, ABUW-RS05155, rpsF 
#> PC_ 3 
#> Positive:  ABUW-RS14490, basJ, tcuB, hutU, ABUW-RS02205, ABUW-RS09975, ABUW-RS05480, ABUW-RS12220, ABUW-RS18155, ABUW-RS01895 
#>     hutH, ABUW-RS11930, ABUW-RS05310, ABUW-RS07640, rho, macA, ABUW-RS18750, macB, ABUW-RS02965, ABUW-RS17650 
#>     ABUW-RS10205, ABUW-RS15565, coaBC, ABUW-RS05420, ABUW-RS15700, ABUW-RS01690, ABUW-RS06595, ABUW-RS18705, ABUW-RS02490, ABUW-RS00165 
#> Negative:  ABUW-RS11870, rpsC, ABUW-RS16090, katE, rpsQ, secY, rplP, rpsK, rpsB, cydB 
#>     nuoG, ABUW-RS06040, rplV, rpmC, ABUW-RS11840, atpG, rplF, tuf, rpsM, rplR 
#>     nuoF, ABUW-RS08740, rpoA, ABUW-RS11875, rplB, atpD, rplQ, ABUW-RS13120, ABUW-RS17440, rpsH 
#> PC_ 4 
#> Positive:  ABUW-RS18750, ABUW-RS02965, ABUW-RS05420, ABUW-RS02985, macA, dnaK, macB, ABUW-RS18450, ABUW-RS18745, ABUW-RS12910 
#>     ABUW-RS07640, ABUW-RS15890, rplU, dkgB, baeS, ABUW-RS08975, uvrB, ABUW-RS04480, rpsF, sthA 
#>     ABUW-RS05105, sucC, groL, ABUW-RS03195, lolA, ABUW-RS02175, dapD, ABUW-RS10315, ABUW-RS14475, bfr.1 
#> Negative:  ABUW-RS11920, ABUW-RS11935, ABUW-RS11930, ABUW-RS11945, ABUW-RS11960, ABUW-RS11940, hppD, ABUW-RS10180, ABUW-RS00045, ABUW-RS12975 
#>     ABUW-RS00350, ABUW-RS01710, ABUW-RS11950, hutH, ABUW-RS00335, ABUW-RS11340, hutU, ABUW-RS04045, ABUW-RS11345, ABUW-RS05480 
#>     ABUW-RS07925, ABUW-RS11335, ABUW-RS10165, ABUW-RS11925, ABUW-RS08650, ABUW-RS01720, ABUW-RS00110, ABUW-RS10190, pta, fahA 
#> PC_ 5 
#> Positive:  ABUW-RS18750, macA, ABUW-RS02965, lolA, ABUW-RS05420, ABUW-RS18155, macB, ABUW-RS01895, ABUW-RS02985, ABUW-RS07640 
#>     ABUW-RS13935, fahA, ABUW-RS11930, lolB, hppD, ABUW-RS11870, tilS, baeS, ABUW-RS06025, ABUW-RS11635 
#>     ABUW-RS11495, ABUW-RS17125, ABUW-RS05105, ispE, ABUW-RS18745, rimO, putA, coaBC, ABUW-RS17440, ABUW-RS04855 
#> Negative:  ABUW-RS02205, ABUW-RS14490, ABUW-RS03160, atpG, ABUW-RS09275, rplV, ABUW-RS09285, ABUW-RS20235, ABUW-RS19140, rpsC 
#>     rlmN, ABUW-RS14475, ABUW-RS15335, rplP, ABUW-RS06145, ABUW-RS05240, ABUW-RS02900, ABUW-RS05875, ABUW-RS11620, ABUW-RS02765 
#>     adeB, bauA, ABUW-RS00775, ABUW-RS16690, rpmJ, ABUW-RS14460, ABUW-RS08700, ABUW-RS04040, rpsS, ABUW-RS04705
#> Computing nearest neighbor graph
#> Computing SNN
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
#> 13:30:15 UMAP embedding parameters a = 0.9922 b = 1.112
#> 13:30:15 Read 1500 rows and found 10 numeric columns
#> 13:30:15 Using Annoy for neighbor search, n_neighbors = 30
#> 13:30:15 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> **************************************************|
#> 13:30:15 Writing NN index file to temp file /var/folders/rm/kkt8lwv515qdxks_k5nc00mh0000gp/T//RtmpENHApB/file141dc1d55b835
#> 13:30:15 Searching Annoy index using 1 thread, search_k = 3000
#> 13:30:15 Annoy recall = 100%
#> 13:30:15 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#> 13:30:16 Initializing from normalized Laplacian + noise (using RSpectra)
#> 13:30:16 Commencing optimization for 500 epochs, with 60244 positive edges
#> 13:30:16 Using rng type: pcg
#> 13:30:17 Optimization finished
seu1$batch <- "batch1"

seu2 <- RunDimRed(seu2)
#> Normalizing layer: counts
#> Finding variable features for layer counts
#> Centering and scaling data matrix
#> PC_ 1 
#> Positive:  ABUW-RS10310, ABUW-RS04625, ABUW-RS01045, ABUW-RS18025, ABUW-RS01770, ABUW-RS02800, ABUW-RS19590, benA, ABUW-RS19140, dmeF 
#>     ABUW-RS08235, ABUW-RS16150, ABUW-RS12715, ABUW-RS09580, ABUW-RS14180, ABUW-RS15320, ABUW-RS12810, ABUW-RS01150, entE, ABUW-RS01115 
#>     barA, ABUW-RS02740, garD, ABUW-RS03365, ABUW-RS15750, atzF, add, ABUW-RS17640, ABUW-RS00350, ABUW-RS03420 
#> Negative:  secY, ABUW-RS03195, ABUW-RS13935, rplB, rplC, fusA, rpsC, rpsF, rpsD, ABUW-RS00825 
#>     ABUW-RS04280, rplE, rpsM, rpsK, rplN, rpsA, ABUW-RS19280, ABUW-RS14775, rplD, rplX 
#>     rpoA, macA, rpsB, rpsH, rplF, tuf, rpsL, rpsE, rpsJ, cyoB 
#> PC_ 2 
#> Positive:  rpsB, rplB, rpsC, rpsS, secY, rplW, rpmC, rplF, rpsL, rpsD 
#>     rpsK, rpsJ, rplC, ABUW-RS19280, rpoA, rplR, ABUW-RS03195, rpsF, rplU, rplO 
#>     rplP, rpmJ, ABUW-RS10685, infB, rpsA, mnmG, rplD, rplV, rplK, cydB 
#> Negative:  ABUW-RS04435, ABUW-RS15345, ABUW-RS18115, ABUW-RS05885, ABUW-RS13745, tolA, gspD, alaS, rdgB, ABUW-RS17780 
#>     ABUW-RS17335, ABUW-RS00405, ABUW-RS15315, ABUW-RS01035, traN, ABUW-RS17300, lnt, mdcH, ABUW-RS05910, ABUW-RS16970 
#>     holA, ABUW-RS00855, ABUW-RS19165, ABUW-RS00170, lolD, ABUW-RS01250, hemH, ABUW-RS16850, ABUW-RS07415, ABUW-RS19140 
#> PC_ 3 
#> Positive:  ssuC, crp, ABUW-RS01035, ABUW-RS15595, ABUW-RS11755, ABUW-RS05885, ABUW-RS13135, thiM, traW, ABUW-RS00995 
#>     ABUW-RS14160, ABUW-RS17335, puuE, nusA, rluB, leuC, ABUW-RS10960, ABUW-RS16970, ABUW-RS17835, ABUW-RS06305 
#>     ABUW-RS15200, ABUW-RS16920, ABUW-RS17425, iscU, ABUW-RS11265, ABUW-RS12670, ggt, dksA, ABUW-RS00060, ABUW-RS05325 
#> Negative:  pncA, ABUW-RS00045, ABUW-RS15555, ABUW-RS17325, rdgB, baeS, ABUW-RS11485, ABUW-RS03790, ABUW-RS08770, ABUW-RS18950 
#>     alaS, pobR, ABUW-RS16450, dnaK, ABUW-RS08835, folD, ABUW-RS13370, ABUW-RS14520, purF, ABUW-RS19805 
#>     ABUW-RS13330, ABUW-RS07810, rseP, ABUW-RS05625, ABUW-RS16790, ABUW-RS19140, ABUW-RS10705, hcaR, ABUW-RS06665, ABUW-RS05400 
#> PC_ 4 
#> Positive:  ABUW-RS14540, holA, ABUW-RS05910, ABUW-RS10325, katE, ABUW-RS19805, bauF, lspA, ABUW-RS13830, tolA 
#>     ABUW-RS16970, ABUW-RS13245, benA, mrdA, ABUW-RS10830, dkgB, ABUW-RS15435, ABUW-RS05885, ABUW-RS04435, ABUW-RS03325 
#>     ABUW-RS16975, ABUW-RS14990, ABUW-RS07465, map, rbtA, ABUW-RS03485, ppc, ABUW-RS04790, ABUW-RS10200, ABUW-RS15820 
#> Negative:  ABUW-RS17325, ABUW-RS14695, ABUW-RS00405, ABUW-RS14410, ABUW-RS11945, ABUW-RS05625, ABUW-RS14750, alr.1, ABUW-RS18310, entE 
#>     pdxA, ABUW-RS15345, ABUW-RS12180, ABUW-RS04740, hchA, ABUW-RS15750, ABUW-RS15595, ABUW-RS16480, ABUW-RS05375, csuC 
#>     ABUW-RS17835, ABUW-RS10315, ABUW-RS10705, ABUW-RS04800, secB, ABUW-RS11755, ABUW-RS17275, ABUW-RS09505, ABUW-RS05965, ABUW-RS16460 
#> PC_ 5 
#> Positive:  ABUW-RS14040, lnt, ABUW-RS12010, aroC, eutC, ABUW-RS20015, crp, ABUW-RS13610, ABUW-RS13160, ABUW-RS07810 
#>     ABUW-RS01170, ABUW-RS03345, ABUW-RS05885, mrdA, ABUW-RS15280, ABUW-RS16790, ABUW-RS13345, ABUW-RS18030, rnc, ABUW-RS14650 
#>     ABUW-RS03140, ABUW-RS03485, ttcA, ABUW-RS15760, ABUW-RS01070, ABUW-RS05400, ABUW-RS13460, ABUW-RS01305, ABUW-RS19140, eat 
#> Negative:  ABUW-RS05185, ABUW-RS01255, ABUW-RS18750, secF, ABUW-RS04215, ABUW-RS16805, bauC, ABUW-RS14895, ABUW-RS17700, ABUW-RS13325 
#>     ABUW-RS05255, ABUW-RS02130, add, ABUW-RS02965, ABUW-RS11920, ABUW-RS06665, ABUW-RS16905, lon, ABUW-RS02475, ABUW-RS04665 
#>     ABUW-RS17810, glnD, ABUW-RS10830, ABUW-RS02255, ABUW-RS16605, ABUW-RS18780, hutH, ABUW-RS17275, ABUW-RS02165, ABUW-RS00610 
#> Computing nearest neighbor graph
#> Computing SNN
#> 13:30:18 UMAP embedding parameters a = 0.9922 b = 1.112
#> 13:30:18 Read 1500 rows and found 10 numeric columns
#> 13:30:18 Using Annoy for neighbor search, n_neighbors = 30
#> 13:30:18 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> **************************************************|
#> 13:30:19 Writing NN index file to temp file /var/folders/rm/kkt8lwv515qdxks_k5nc00mh0000gp/T//RtmpENHApB/file141dc73902c6b
#> 13:30:19 Searching Annoy index using 1 thread, search_k = 3000
#> 13:30:19 Annoy recall = 100%
#> 13:30:19 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#> 13:30:20 Initializing from normalized Laplacian + noise (using RSpectra)
#> 13:30:20 Commencing optimization for 500 epochs, with 62430 positive edges
#> 13:30:20 Using rng type: pcg
#> 13:30:21 Optimization finished
seu2$batch <- "batch2"

seu3 <- RunDimRed(seu3)
#> Normalizing layer: counts
#> Finding variable features for layer counts
#> Centering and scaling data matrix
#> PC_ 1 
#> Positive:  ABUW-RS11870, ABUW-RS19720, katE, ABUW-RS02205, ABUW-RS14490, ABUW-RS15165, ABUW-RS11875, ABUW-RS10640, ABUW-RS11840, ABUW-RS10635 
#>     entE, basD, ABUW-RS10585, basB, ABUW-RS15280, ABUW-RS11850, ABUW-RS06610, ABUW-RS05995, ABUW-RS08700, ABUW-RS07960 
#>     ABUW-RS11860, ABUW-RS13010, ABUW-RS00895, ABUW-RS10625, ABUW-RS13120, ABUW-RS11960, ABUW-RS10620, ABUW-RS13250, raiA, ABUW-RS14025 
#> Negative:  rplJ, rpsC, rplB, rpsM, rpsD, rpsB, rpoA, rplC, rpsK, secY 
#>     rplF, rplD, rplX, fusA, omp33-36, atpD, rplE, rplK, rpsE, rplM 
#>     ABUW-RS18175, rplA, rplP, rpsJ, rplV, atpG, rplN, ABUW-RS17440, rplL, infB 
#> PC_ 2 
#> Positive:  ABUW-RS02205, ABUW-RS10640, ABUW-RS14490, entE, basD, basB, ABUW-RS10635, ABUW-RS19720, ABUW-RS10625, ABUW-RS10585 
#>     basJ, ABUW-RS10620, ABUW-RS10220, ABUW-RS06610, basG, ABUW-RS10520, ABUW-RS15460, basC, basF, ABUW-RS10630 
#>     ABUW-RS01870, ABUW-RS15470, ABUW-RS17685, barB, raiA, ABUW-RS14495, ABUW-RS18090, kdpA, ABUW-RS14485, ABUW-RS10705 
#> Negative:  katE, ABUW-RS11840, ABUW-RS11875, ABUW-RS15165, ABUW-RS03090, ABUW-RS05995, ABUW-RS14025, ABUW-RS07960, ABUW-RS11850, ABUW-RS09995 
#>     ABUW-RS11860, ABUW-RS11870, ABUW-RS12725, ABUW-RS18980, ABUW-RS12605, ABUW-RS17440, ABUW-RS17145, ABUW-RS13010, ABUW-RS14240, ABUW-RS07965 
#>     ABUW-RS07155, ABUW-RS13220, ABUW-RS17420, ABUW-RS15520, ABUW-RS15280, lysM, ABUW-RS07150, ABUW-RS16155, ABUW-RS01020, otsB 
#> PC_ 3 
#> Positive:  ABUW-RS02205, ABUW-RS14490, ABUW-RS19720, ABUW-RS00660, ABUW-RS03350, raiA, ABUW-RS13815, ABUW-RS13135, ABUW-RS03355, ABUW-RS01495 
#>     ABUW-RS15470, pilB, ABUW-RS15460, ABUW-RS18175, sucD, ABUW-RS17440, ABUW-RS06610, ABUW-RS01715, ABUW-RS14715, ABUW-RS10640 
#>     ABUW-RS13140, ABUW-RS01445, nuoF, ABUW-RS05265, tuf, ABUW-RS03345, ABUW-RS14485, cas6f, ABUW-RS08585, ABUW-RS02765 
#> Negative:  ABUW-RS01870, ABUW-RS16100, ABUW-RS11870, ABUW-RS15035, rplC, ABUW-RS13090, trmD, ABUW-RS16090, mnmG, rho 
#>     ABUW-RS10700, rimM, ABUW-RS15170, ABUW-RS09970, infB, ABUW-RS16655, ABUW-RS05155, ABUW-RS09975, rplD, pth 
#>     ABUW-RS05405, tssC, ABUW-RS15030, ABUW-RS10220, rpsJ, ABUW-RS13330, gltP, ABUW-RS07745, ABUW-RS12760, rpsF 
#> PC_ 4 
#> Positive:  rpsS, secY, rpsE, rpsC, atpD, rpmJ, atpG, adeB, mqo, rpsH 
#>     rplJ, rplB, nuoF, ABUW-RS06040, ABUW-RS03360, ABUW-RS15920, rplU, rplR, rpsB, rpsK 
#>     rplF, rplP, rpsR, rpsN, ABUW-RS18175, ABUW-RS02930, nuoM, rpmB, ABUW-RS05910, ptsP.1 
#> Negative:  rph, grxC, ABUW-RS10220, ABUW-RS09975, basD, basB, ABUW-RS05405, entE, ABUW-RS10640, ABUW-RS15460 
#>     ABUW-RS18705, yidD, ABUW-RS01250, hppD, ubiE, basJ, benA, ABUW-RS17170, ABUW-RS10520, ABUW-RS01720 
#>     ABUW-RS10700, ABUW-RS18710, fumC, basF, basG, ABUW-RS15280, ABUW-RS15250, phoR, ABUW-RS01480, ABUW-RS10635 
#> PC_ 5 
#> Positive:  pilB, groL, ABUW-RS12860, nuoF, ABUW-RS01250, ABUW-RS13815, pepN, ABUW-RS03350, hscA, ABUW-RS13135 
#>     ABUW-RS11870, ABUW-RS03355, ABUW-RS13140, ychF, ABUW-RS05120, ABUW-RS04500, ABUW-RS10220, ABUW-RS03345, ABUW-RS06740, ABUW-RS15350 
#>     grpE, eno, ABUW-RS03075, ptsP.1, ABUW-RS07595, cueR, ABUW-RS17090, grxC, ABUW-RS07925, aroE 
#> Negative:  ABUW-RS18625, betA, ABUW-RS09965, ABUW-RS02205, ABUW-RS04495, ABUW-RS14490, ABUW-RS15470, trmD, ABUW-RS09975, ABUW-RS13250 
#>     infB, cmlA, ABUW-RS19720, mqo, ABUW-RS17200, ABUW-RS19305, ABUW-RS13680, ABUW-RS02935, nusA, ABUW-RS18615 
#>     ABUW-RS18705, ABUW-RS15280, rplO, betT, rplI, ABUW-RS07515, ABUW-RS00895, gdhA, ABUW-RS01315, ABUW-RS04830 
#> Computing nearest neighbor graph
#> Computing SNN
#> 13:30:22 UMAP embedding parameters a = 0.9922 b = 1.112
#> 13:30:22 Read 2000 rows and found 10 numeric columns
#> 13:30:22 Using Annoy for neighbor search, n_neighbors = 30
#> 13:30:22 Building Annoy index with metric = cosine, n_trees = 50
#> 0%   10   20   30   40   50   60   70   80   90   100%
#> [----|----|----|----|----|----|----|----|----|----|
#> **************************************************|
#> 13:30:22 Writing NN index file to temp file /var/folders/rm/kkt8lwv515qdxks_k5nc00mh0000gp/T//RtmpENHApB/file141dc62afe986
#> 13:30:22 Searching Annoy index using 1 thread, search_k = 3000
#> 13:30:23 Annoy recall = 100%
#> 13:30:23 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
#> 13:30:24 Initializing from normalized Laplacian + noise (using RSpectra)
#> 13:30:24 Commencing optimization for 500 epochs, with 85436 positive edges
#> 13:30:24 Using rng type: pcg
#> 13:30:25 Optimization finished
seu3$batch <- "batch3"
```

## Create SCProkaR objects

[`CreateBacObject()`](https://monash-li-lab.github.io/scProka/reference/CreateBacObject.md)
can extract counts, metadata, and existing reductions such as `pca` and
`umap` directly from each Seurat object.

``` r
sce1 <- CreateBacObject(
  seu1,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "sample",
  batch_col = "batch"
)

sce2 <- CreateBacObject(
  seu2,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "sample",
  batch_col = "batch"
)

sce3 <- CreateBacObject(
  seu3,
  seurat_assay = "RNA",
  seurat_layer = "counts",
  sample_col = "sample",
  batch_col = "batch"
)
```

## Quality control

Compute QC metrics on each sample separately:

``` r
sce1 <- RunBacQC(sce1)
sce2 <- RunBacQC(sce2)
sce3 <- RunBacQC(sce3)
```

The QC results are stored in `colData(sce)`:

``` r
head(as.data.frame(colData(sce1))[, c(
  "total_counts",
  "detected_features",
  "rrna_fraction",
  "pct_rrna"
)])
#>                           total_counts detected_features rrna_fraction pct_rrna
#> P1_1_GGAGCTGAGAGGTCCGTCTA          504               167             0        0
#> P1_1_TATTGTCATCGGAGCGCTGC          583               176             0        0
#> P1_1_GTATCTTCATCGACCGTAAT          377                98             0        0
#> P1_1_GCAGGCAATCGGAGCGATAT          428               138             0        0
#> P1_1_GCAGGTAAGGCGCCTGGAGC          538               168             0        0
#> P1_1_GGTCCAGTTGGCGAAGATGA         1544               359             0        0
```

Filter cells using your chosen thresholds:

``` r
sce1 <- FilterBacCells(
  sce1,
  min_counts = 50,
  min_features = 50,
  max_rrna_fraction = 0.20
)

sce2 <- FilterBacCells(
  sce2,
  min_counts = 50,
  min_features = 50,
  max_rrna_fraction = 0.20
)

sce3 <- FilterBacCells(
  sce3,
  min_counts = 50,
  min_features = 50,
  max_rrna_fraction = 0.20
)
```

## Merge samples

Use
[`MergeBacObjects()`](https://monash-li-lab.github.io/scProka/reference/MergeBacObjects.md)
instead of [`cbind()`](https://rdrr.io/r/base/cbind.html) to avoid
duplicate-barcode issues and to align genes safely across samples.

``` r
sce12 <- MergeBacObjects(sce1, sce2, gene_mode = "intersect")
sce <- MergeBacObjects(sce12, sce3, gene_mode = "intersect")
```

If cell barcodes overlap between samples,
[`MergeBacObjects()`](https://monash-li-lab.github.io/scProka/reference/MergeBacObjects.md)
automatically makes them unique and stores the original barcode in
`colData(sce)$original_cell_id`.

If both input objects contain compatible reductions such as `pca` and
`umap`,
[`MergeBacObjects()`](https://monash-li-lab.github.io/scProka/reference/MergeBacObjects.md)
carries them forward into the merged object.

## Visualize unintegrated structure

If your Seurat objects already contained a UMAP, it is typically
available after conversion:

``` r
PlotReduction(sce, reduction = "umap", colour_by = "batch")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-10-1.png)

``` r
PlotReduction(sce, reduction = "umap", colour_by = "clusters")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-10-2.png)

Interpretation:

- colouring by `batch` shows whether samples are separated before
  integration
- colouring by `seurat_clusters` or another biological label shows
  whether expected structure is present

`PlotReduction(..., facet_by = "batch")` is the `SCProkaR` equivalent of
Seurat’s `split.by` visualization.

## Run integration

### Option 1: MNN

``` r
sce <- IntegrateBacData(
  sce,
  batch_col = "batch",
  method = "mnn",
  dims = 1:30
)
```

### Option 2: Harmony

``` r
sce <- IntegrateBacData(
  sce,
  batch_col = "batch",
  method = "harmony",
  dims = 1:30
)
```

The output embedding is stored as `integrated_mnn` or
`integrated_harmony`.

## Visualize the integrated latent space

``` r
PlotReduction(sce, reduction = "integrated_mnn", colour_by = "batch")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-13-1.png)

``` r
PlotReduction(sce, reduction = "integrated_harmony", colour_by = "batch")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-13-2.png)

This gives a quick direct view of the integrated latent space before
running a new UMAP.

## Run UMAP on the integrated embedding

[`RunIntegratedUMAP()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedUMAP.md)
computes a UMAP from any stored embedding:

``` r
sce <- RunIntegratedUMAP(sce, reduction = "integrated_mnn")
sce <- RunIntegratedUMAP(sce, reduction = "integrated_harmony")
```

This stores a new reduced dimension named `umap_integrated_mnn`.

``` r
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "batch")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-15-1.png)

``` r
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "clusters")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-15-2.png)

``` r
PlotReduction(sce, reduction = "umap_integrated_harmony", colour_by = "batch")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-15-3.png)

``` r
PlotReduction(sce, reduction = "umap_integrated_harmony", colour_by = "clusters")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-15-4.png)

## Cluster cells after integration

Future downstream analyses should use clusters identified from the
integrated embedding rather than the pre-integration batch structure.

[`RunIntegratedClustering()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedClustering.md)
builds a kNN graph on the corrected embedding and stores the resulting
labels in `colData(sce)`.

``` r
sce <- RunIntegratedClustering(
  sce,
  reduction = "integrated_mnn",
  cluster_col = "integrated_clusters",
  k = 20,
  resolution = 1
)
```

Here, `k` controls the neighbour graph size, while `resolution` controls
how finely that graph is partitioned into clusters. This mirrors the
Seurat logic of `FindNeighbors(k.param = ...)` followed by
`FindClusters(resolution = ...)`.

Visualize the integrated clusters on the corrected UMAP:

``` r
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "integrated_clusters")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-17-1.png)

``` r
PlotReduction(
  sce,
  reduction = "umap_integrated_mnn",
  colour_by = "integrated_clusters",
  facet_by = "batch"
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-17-2.png)

## Seurat-style integration summary panels

The Seurat integration tutorial emphasizes three visual checks after
alignment:

- compare unintegrated and integrated embeddings coloured by batch
- compare unintegrated and integrated embeddings coloured by cluster or
  label
- split the integrated embedding by batch/condition while keeping
  cluster colours fixed

`SCProkaR` provides
[`PlotIntegrationOverview()`](https://monash-li-lab.github.io/scProka/reference/PlotIntegrationOverview.md)
for this pattern.

### MNN example

``` r
mnn_plots <- PlotIntegrationOverview(
  sce,
  unintegrated_reduction = "umap",
  integrated_reduction = "umap_integrated_mnn",
  batch_col = "batch",
  label_col = "integrated_clusters",
  split_by = "batch"
)

mnn_plots$unintegrated_by_batch
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-18-1.png)

``` r
mnn_plots$integrated_by_batch
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-18-2.png)

``` r
mnn_plots$unintegrated_by_label
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-18-3.png)

``` r
mnn_plots$integrated_by_label
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-18-4.png)

``` r
mnn_plots$integrated_split
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-18-5.png)

### Harmony example

``` r
harmony_plots <- PlotIntegrationOverview(
  sce,
  unintegrated_reduction = "umap",
  integrated_reduction = "umap_integrated_harmony",
  batch_col = "batch",
  label_col = "integrated_clusters",
  split_by = "batch"
)

harmony_plots$unintegrated_by_batch
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-19-1.png)

``` r
harmony_plots$integrated_by_batch
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-19-2.png)

``` r
harmony_plots$unintegrated_by_label
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-19-3.png)

``` r
harmony_plots$integrated_by_label
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-19-4.png)

``` r
harmony_plots$integrated_split
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-19-5.png)

Interpretation:

- `unintegrated_by_batch` should show stronger sample separation than
  `integrated_by_batch`
- `integrated_by_label` should preserve biologically meaningful group
  separation
- `integrated_split` is the closest equivalent to Seurat’s
  `DimPlot(..., split.by = "stim")`

After successful integration:

- batches should mix better than in the original UMAP
- biological clusters should remain interpretable

## Benchmark integration quality

If you have a biological label such as `seurat_clusters`, you can
benchmark the integration with the package’s scIB-inspired metric panel:

``` r
bench <- BenchmarkIntegration(
  sce,
  batch_col = "batch",
  label_col = "clusters",
  methods = NULL
)

bench$ranking
#>    method     score rank
#> 1 harmony 0.5683985    1
#> 2     mnn 0.5504337    2
bench$plots$overall
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-1.png)

``` r
bench$plots$key_metrics
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-2.png)

``` r
bench$plots$batch_removal
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-3.png)

``` r
bench$plots$bio_conservation
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-4.png)

``` r
bench$plots$tradeoff
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-5.png)

``` r
bench$plots$heatmap
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-20-6.png)

The key comparison figures are:

- `key_metrics`: per-metric score comparison across methods
- `batch_removal`: average of the batch-removal metrics
- `bio_conservation`: average of the biological-conservation metrics
- `tradeoff`: batch-removal versus bio-conservation, which is often the
  most useful summary of integration quality

## Pseudobulk differential expression

`SCProkaR` also includes pseudobulk differential-expression utilities
built around edgeR. Here we continue with the same example dataset used
above rather than switching to a separate simulated object.

The merged `sce` object currently contains one sample per
treatment-timepoint combination, with:

- `treatments` describing the perturbation group
- `timepoints` describing the experimental sampling time
- `sample` identifying each original sample

For this reason, we create pseudo-replicates inside each sample so the
edgeR workflow has replicate-like pseudobulk columns to operate on. This
is useful for demonstrating the mechanics of the workflow, but in a real
study you should aggregate by the true biological replicate, donor, or
independent sample.

### Visualize the dataset before pseudobulk aggregation

``` r
sce$timepoints_factor <- factor(
  sce$timepoints,
  levels = sort(unique(sce$timepoints))
)

PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "treatments")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-21-1.png)

``` r
PlotReduction(sce, reduction = "umap_integrated_mnn", colour_by = "timepoints_factor")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-21-2.png)

### Create pseudo-replicate sample labels

We split cells inside each original sample into pseudo-replicates. This
keeps the DE setup aligned with the real dataset already used in the
vignette.

``` r
assign_pseudoreplicates <- function(
    sce,
    group_cols = c("sample", "treatments", "timepoints"),
    n_reps = 3,
    seed = 101
) {
  set.seed(seed)

  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  strata <- interaction(meta[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  rep_id <- integer(ncol(sce))

  for (level_name in levels(strata)) {
    idx <- which(strata == level_name)
    rep_id[idx] <- sample(rep(seq_len(n_reps), length.out = length(idx)))
  }

  sce$pb_replicate <- factor(rep_id, levels = seq_len(n_reps))
  sce$pb_sample_id <- paste(
    as.character(sce$sample),
    paste0("r", sce$pb_replicate),
    sep = "_"
  )

  sce
}

sce_de <- assign_pseudoreplicates(sce, n_reps = 3, seed = 101)

with(
  as.data.frame(SummarizedExperiment::colData(sce_de)),
  table(sample, pb_replicate)
)
#>            pb_replicate
#> sample        1   2   3
#>   PMB0.5_1  167 167 166
#>   PMB0.5_4  167 167 166
#>   PMB0.5_7  167 167 166
#>   PMB2_1    167 167 166
#>   PMB2_4     55  54  54
#>   PMB2_7    125 125 124
#>   control_0 167 167 166
#>   control_1 136 135 135
#>   control_4  88  88  87
#>   control_7  65  65  65
```

### Aggregate pseudobulk samples

[`aggregate_pseudobulk()`](https://monash-li-lab.github.io/scProka/reference/aggregate_pseudobulk.md)
sums counts across cells that share the same sample definition and keeps
the sample-level metadata alongside the aggregated counts. Here we
aggregate the pseudo-replicated object, which preserves all three
treatments while creating replicate-like pseudobulk columns inside each
sample-timepoint stratum.

``` r
pb <- aggregate_pseudobulk(
  sce_de,
  sample_cols = c("pb_sample_id", "sample", "pb_replicate", "clusters", "treatments", "timepoints"),
  aggregation = "sum",
  min_cells = 20
)

pb <- normalize_pseudobulk(pb, method = "TMM")
pb$timepoints_factor <- factor(pb$timepoints, levels = sort(unique(pb$timepoints)))

pb_pca <- stats::prcomp(
  t(as.matrix(SummarizedExperiment::assay(pb, "logcounts"))),
  scale. = TRUE
)
SingleCellExperiment::reducedDim(pb, "pb_pca") <- pb_pca$x[, 1:10, drop = FALSE]
```

``` r
PlotReduction(pb, reduction = "pb_pca", colour_by = "treatments", point_size = 5)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-24-1.png)

``` r
PlotReduction(pb, reduction = "pb_pca", colour_by = "timepoints_factor", point_size = 5)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-24-2.png)

### Pairwise differential expression

The first DE mode compares one group against another. A common default
is to compare each timepoint against baseline while adjusting for
treatment as a covariate. In this dataset, the `0` timepoint is only
present in the control group, so these contrasts should be interpreted
as overall progression away from the untreated baseline rather than a
perfectly balanced treatment-adjusted design. With all three treatments
present in the model, the treatment term absorbs broad
perturbation-specific shifts while the timepoint contrast captures the
shared movement away from baseline.

``` r
pairwise_time <- run_edger_pairwise_de(
  pb,
  group_col = "timepoints",
  covariates = "treatments",
  contrast_type = "reference",
  reference_level = "0",
  normalization = "TMM"
)

names(pairwise_time$tables)
#> [1] "1_vs_0" "4_vs_0" "7_vs_0"
```

``` r
pairwise_sig_counts <- aggregate(
  significant ~ contrast,
  data = pairwise_time$combined_table,
  FUN = sum
)

ggplot2::ggplot(pairwise_sig_counts, ggplot2::aes(x = contrast, y = significant, fill = contrast)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title = "Significant genes in timepoint-versus-baseline contrasts",
    x = "contrast",
    y = "number of DE genes"
  )
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-26-1.png)

Inspect one contrast in more detail:

``` r
head(
  pairwise_time$tables[["1_vs_0"]][
    order(pairwise_time$tables[["1_vs_0"]]$FDR),
    c("gene", "logFC", "FDR")
  ],
  10
)
#>                      gene     logFC          FDR
#> ABUW-RS11870 ABUW-RS11870  9.429892 4.085509e-26
#> ABUW-RS11875 ABUW-RS11875  9.812129 2.952245e-18
#> ABUW-RS15330 ABUW-RS15330  5.767058 5.246736e-18
#> omp33-36         omp33-36 -4.232055 1.812547e-17
#> ABUW-RS11840 ABUW-RS11840  9.395086 9.193228e-17
#> carO                 carO  4.200643 9.231608e-17
#> ABUW-RS20195 ABUW-RS20195  7.086212 1.868548e-16
#> ABUW-RS18175 ABUW-RS18175 -3.297693 4.409101e-16
#> rplJ                 rplJ -3.923307 5.592220e-16
#> katE                 katE  6.341821 9.201106e-16

plot_pairwise_de_volcano(
  pairwise_time,
  contrast = "1_vs_0",
  fdr_cutoff = 0.05,
  lfc_cutoff = 0.5,
  top_n_labels = 8
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-27-1.png)

``` r

plot_pairwise_de_ma(
  pairwise_time,
  contrast = "1_vs_0",
  fdr_cutoff = 0.05
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-27-2.png)

`SCProkaR` also supports composite treatment-timepoint groups when you
want an explicit within-treatment comparison such as `PMB0.5 at day 7`
versus `PMB0.5 at day 1`, while keeping all treatment groups in the same
pseudobulk object.

``` r
pairwise_composite <- run_edger_pairwise_de(
  pb,
  group_col = "treatments",
  combine_group_cols = c("treatments", "timepoints"),
  contrast_type = "manual",
  manual_pairs = data.frame(
    group1 = c("PMB0.5__7", "PMB2__7", "control__7", "PMB0.5__4"),
    group2 = c("PMB0.5__1", "PMB2__1", "control__1", "PMB0.5__1"),
    stringsAsFactors = FALSE
  ),
  normalization = "TMM"
)

plot_pairwise_de_volcano(
  pairwise_composite,
  contrast = "PMB0.5__7_vs_PMB0.5__1",
  fdr_cutoff = 0.05,
  lfc_cutoff = 0.5,
  top_n_labels = 8
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-28-1.png)

### Spline-based time-series differential expression

The second DE mode treats time as continuous and uses a spline basis to
detect genes that change over time within each treatment group. Because
the example object now contains all three treatments, the result returns
a separate spline fit for each perturbation trajectory.

``` r
spline_de <- run_edger_spline_de(
  pb,
  time_col = "timepoints",
  condition_col = "treatments",
  df = 3,
  normalization = "TMM",
  return_curve_for = "significant"
)

names(spline_de$by_condition)
#> [1] "control" "PMB0.5"  "PMB2"
```

``` r
spline_sig_counts <- aggregate(
  significant ~ treatments,
  data = spline_de$combined_table,
  FUN = sum
)

ggplot2::ggplot(spline_sig_counts, ggplot2::aes(x = treatments, y = significant, fill = treatments)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title = "Genes with significant temporal change within each treatment",
    x = "treatment",
    y = "number of spline DE genes"
  )
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-30-1.png)

``` r
top_by_condition <- lapply(
  split(spline_de$combined_table, spline_de$combined_table$treatments),
  function(df) {
    keep_cols <- intersect(c("gene", "F", "PValue", "FDR"), colnames(df))
    utils::head(df[order(df$FDR), keep_cols, drop = FALSE], 100)
  }
)

top_by_condition
#> $control
#>                              gene         F       PValue          FDR
#> control.ABUW-RS11870 ABUW-RS11870 198.52145 8.643393e-66 2.800459e-63
#> control.katE                 katE 178.60583 3.470688e-61 5.622515e-59
#> control.ABUW-RS02205 ABUW-RS02205 175.00746 1.573835e-60 1.699742e-58
#> control.ABUW-RS14490 ABUW-RS14490 138.62055 7.913372e-54 6.409831e-52
#> control.ABUW-RS10640 ABUW-RS10640 106.34761 1.039887e-44 6.738470e-43
#> control.ABUW-RS11840 ABUW-RS11840 108.56957 1.340614e-44 7.239316e-43
#> control.ABUW-RS15165 ABUW-RS15165 102.59970 3.689458e-43 1.707692e-41
#> control.ABUW-RS11875 ABUW-RS11875 102.98549 9.282119e-43 3.759258e-41
#> control.ABUW-RS15330 ABUW-RS15330  99.03895 6.053414e-42 2.179229e-40
#> control.entE                 entE  99.06560 1.158112e-41 3.752283e-40
#> control.rplJ                 rplJ  81.33042 1.431876e-36 4.217527e-35
#> control.ABUW-RS19720 ABUW-RS19720  79.65804 6.549129e-36 1.768265e-34
#> control.basD                 basD  77.69868 1.959697e-35 4.884168e-34
#> control.omp33-36         omp33-36  75.69790 3.964654e-35 9.175343e-34
#> control.ABUW-RS07500 ABUW-RS07500  75.67349 7.201328e-35 1.555487e-33
#> control.rpsC                 rpsC  75.86185 1.045011e-34 2.116148e-33
#> control.carO                 carO  73.09841 1.012801e-33 1.930280e-32
#> control.ABUW-RS10585 ABUW-RS10585  71.51752 3.871232e-33 6.968218e-32
#> control.infB                 infB  66.54358 1.704688e-31 2.906941e-30
#> control.ABUW-RS11850 ABUW-RS11850  67.11078 1.974124e-31 3.016775e-30
#> control.ABUW-RS05995 ABUW-RS05995  66.32347 2.027674e-31 3.016775e-30
#> control.secY                 secY  66.49198 2.048428e-31 3.016775e-30
#> control.basJ                 basJ  66.55031 2.276951e-31 3.207530e-30
#> control.basB                 basB  64.96492 4.878911e-31 6.586530e-30
#> control.ABUW-RS18175 ABUW-RS18175  64.72693 8.444214e-31 1.094370e-29
#> control.ABUW-RS07960 ABUW-RS07960  63.26662 2.502163e-30 3.118080e-29
#> control.ABUW-RS20195 ABUW-RS20195  62.91565 2.974370e-30 3.569244e-29
#> control.rpoA                 rpoA  62.59870 5.078360e-30 5.876388e-29
#> control.ABUW-RS10220 ABUW-RS10220  62.18032 6.294309e-30 7.032263e-29
#> control.ABUW-RS10635 ABUW-RS10635  60.44335 1.148935e-29 1.240850e-28
#> control.rpsD                 rpsD  61.15144 1.714798e-29 1.792240e-28
#> control.ABUW-RS17440 ABUW-RS17440  61.12846 2.161549e-29 2.188568e-28
#> control.ABUW-RS14025 ABUW-RS14025  60.16284 2.938995e-29 2.885559e-28
#> control.fusA                 fusA  60.28005 3.685596e-29 3.512157e-28
#> control.rplD                 rplD  59.93110 4.828295e-29 4.469622e-28
#> control.rplC                 rplC  59.79105 5.646817e-29 5.082135e-28
#> control.atpD                 atpD  58.20976 2.111979e-28 1.849409e-27
#> control.rplE                 rplE  55.58106 2.172529e-27 1.838104e-26
#> control.rplF                 rplF  55.45018 2.212533e-27 1.838104e-26
#> control.ABUW-RS08740 ABUW-RS08740  54.58762 6.286635e-27 5.092174e-26
#> control.ABUW-RS03350 ABUW-RS03350  53.73827 7.573972e-27 5.985285e-26
#> control.rplK                 rplK  53.71237 1.143061e-26 8.817900e-26
#> control.rpsQ                 rpsQ  52.66864 1.790945e-26 1.349456e-25
#> control.rpsB                 rpsB  53.19575 1.846483e-26 1.359683e-25
#> control.atpG                 atpG  52.46238 2.780345e-26 2.001849e-25
#> control.rplP                 rplP  52.31324 3.018355e-26 2.125972e-25
#> control.rpsK                 rpsK  52.47197 3.519489e-26 2.426201e-25
#> control.rplV                 rplV  52.05979 3.961111e-26 2.673750e-25
#> control.rplL                 rplL  51.27920 7.748201e-26 5.123300e-25
#> control.rplA                 rplA  51.10659 1.167644e-25 7.566334e-25
#> control.rpsE                 rpsE  50.64644 1.552398e-25 9.862296e-25
#> control.rpsM                 rpsM  50.81434 1.604785e-25 9.999043e-25
#> control.ABUW-RS14740 ABUW-RS14740  49.41221 4.063051e-25 2.483828e-24
#> control.ABUW-RS11860 ABUW-RS11860  49.00988 5.744116e-25 3.446469e-24
#> control.rplM                 rplM  49.24450 6.113371e-25 3.577388e-24
#> control.ABUW-RS10625 ABUW-RS10625  48.63181 6.183140e-25 3.577388e-24
#> control.rplX                 rplX  48.64796 1.160484e-24 6.596433e-24
#> control.rpsN                 rpsN  47.50723 2.244516e-24 1.253833e-23
#> control.ABUW-RS03355 ABUW-RS03355  47.51308 2.590857e-24 1.422776e-23
#> control.rplN                 rplN  47.74264 2.745773e-24 1.461708e-23
#> control.ABUW-RS06610 ABUW-RS06610  47.04668 2.751981e-24 1.461708e-23
#> control.ABUW-RS03090 ABUW-RS03090  46.73139 5.682611e-24 2.969622e-23
#> control.ABUW-RS07145 ABUW-RS07145  46.91183 6.123260e-24 3.149105e-23
#> control.wecB.2             wecB.2  46.70216 8.951106e-24 4.531497e-23
#> control.nusA                 nusA  45.84752 1.695415e-23 8.450990e-23
#> control.ABUW-RS01495 ABUW-RS01495  44.92179 3.775012e-23 1.853188e-22
#> control.rpsJ                 rpsJ  44.54935 5.664806e-23 2.739399e-22
#> control.ABUW-RS15280 ABUW-RS15280  44.11982 7.539499e-23 3.588728e-22
#> control.rpsF                 rpsF  44.21082 7.680455e-23 3.588728e-22
#> control.rplR                 rplR  43.67258 7.753426e-23 3.588728e-22
#> control.rpsS                 rpsS  43.30005 1.308533e-22 5.971333e-22
#> control.ABUW-RS10620 ABUW-RS10620  42.93177 1.513821e-22 6.812192e-22
#> control.ABUW-RS09995 ABUW-RS09995  42.88024 3.066486e-22 1.361016e-21
#> control.nuoG                 nuoG  41.54322 1.014261e-21 4.440817e-21
#> control.ABUW-RS13010 ABUW-RS13010  40.98123 1.844609e-21 7.968713e-21
#> control.rplU                 rplU  40.48387 2.819221e-21 1.201878e-20
#> control.tuf                   tuf  39.80388 5.924668e-21 2.492977e-20
#> control.rpmB                 rpmB  39.26863 9.938497e-21 4.128299e-20
#> control.rpsH                 rpsH  38.90954 1.129330e-20 4.631681e-20
#> control.rpmC                 rpmC  38.87562 1.320040e-20 5.346162e-20
#> control.atpE                 atpE  38.16530 2.534334e-20 1.013734e-19
#> control.ABUW-RS14775 ABUW-RS14775  37.86517 4.085723e-20 1.614359e-19
#> control.ABUW-RS15035 ABUW-RS15035  37.75573 4.569429e-20 1.783729e-19
#> control.ABUW-RS18615 ABUW-RS18615  37.30612 4.746678e-20 1.830861e-19
#> control.rplO                 rplO  37.49738 5.124970e-20 1.953518e-19
#> control.ABUW-RS11960 ABUW-RS11960  37.28209 7.320731e-20 2.758043e-19
#> control.ABUW-RS18980 ABUW-RS18980  36.90786 7.643796e-20 2.846655e-19
#> control.ABUW-RS05570 ABUW-RS05570  36.61621 1.159260e-19 4.268185e-19
#> control.ABUW-RS16100 ABUW-RS16100  36.59110 1.418381e-19 5.163546e-19
#> control.basF                 basF  36.48706 1.800869e-19 6.483128e-19
#> control.ABUW-RS04495 ABUW-RS04495  36.19087 1.895900e-19 6.750237e-19
#> control.cydB                 cydB  36.19671 1.969467e-19 6.935947e-19
#> control.ABUW-RS18180 ABUW-RS18180  35.72233 2.600763e-19 9.060722e-19
#> control.basG                 basG  35.45592 5.428132e-19 1.870973e-18
#> control.ABUW-RS04045 ABUW-RS04045  34.98868 5.718252e-19 1.950225e-18
#> control.ABUW-RS12405 ABUW-RS12405  35.03018 6.249059e-19 2.109057e-18
#> control.ABUW-RS03380 ABUW-RS03380  34.96297 6.826430e-19 2.280168e-18
#> control.rplW                 rplW  34.96838 7.098554e-19 2.346869e-18
#> control.wecC                 wecC  34.89193 7.665626e-19 2.508750e-18
#> control.ABUW-RS00825 ABUW-RS00825  34.71408 1.046219e-18 3.389750e-18
#> 
#> $PMB0.5
#>                             gene         F       PValue          FDR
#> PMB0.5.ABUW-RS14490 ABUW-RS14490 207.09631 1.506805e-18 1.158733e-15
#> PMB0.5.ABUW-RS12405 ABUW-RS12405 128.70044 5.119385e-18 1.968404e-15
#> PMB0.5.katE                 katE 170.65568 2.683124e-17 6.877742e-15
#> PMB0.5.ABUW-RS09995 ABUW-RS09995  95.00191 1.889356e-16 3.632288e-14
#> PMB0.5.ABUW-RS12725 ABUW-RS12725  91.09269 3.889249e-16 5.968662e-14
#> PMB0.5.ABUW-RS14025 ABUW-RS14025 128.00562 4.656953e-16 5.968662e-14
#> PMB0.5.basJ                 basJ 150.01205 6.682333e-16 7.341020e-14
#> PMB0.5.ABUW-RS02205 ABUW-RS02205 112.87996 7.271422e-15 6.989654e-13
#> PMB0.5.basG                 basG 119.52304 1.334176e-14 1.139980e-12
#> PMB0.5.ABUW-RS11840 ABUW-RS11840  85.43877 1.575019e-14 1.211190e-12
#> PMB0.5.basC                 basC 116.98718 1.753829e-14 1.226086e-12
#> PMB0.5.basD                 basD 104.34243 2.103654e-14 1.348092e-12
#> PMB0.5.tcuA                 tcuA  88.26247 4.073965e-14 2.409907e-12
#> PMB0.5.ABUW-RS18980 ABUW-RS18980  80.35495 4.515103e-14 2.480081e-12
#> PMB0.5.ABUW-RS17440 ABUW-RS17440 123.98479 5.086234e-14 2.607543e-12
#> PMB0.5.tcuB                 tcuB  91.13126 8.454573e-14 3.921108e-12
#> PMB0.5.basF                 basF 102.85174 8.668249e-14 3.921108e-12
#> PMB0.5.ABUW-RS11875 ABUW-RS11875  86.24074 1.341602e-13 5.731621e-12
#> PMB0.5.ABUW-RS15165 ABUW-RS15165  89.24135 1.668712e-13 6.753890e-12
#> PMB0.5.ABUW-RS07960 ABUW-RS07960  68.72588 2.588485e-13 9.952724e-12
#> PMB0.5.ABUW-RS10640 ABUW-RS10640  85.60945 2.842755e-13 1.040990e-11
#> PMB0.5.basB                 basB  84.85434 3.200648e-13 1.118772e-11
#> PMB0.5.ABUW-RS14495 ABUW-RS14495  81.19286 4.308270e-13 1.440461e-11
#> PMB0.5.ABUW-RS10585 ABUW-RS10585  89.18445 4.974837e-13 1.594021e-11
#> PMB0.5.ABUW-RS07170 ABUW-RS07170  61.19434 7.300919e-13 2.245763e-11
#> PMB0.5.ABUW-RS13010 ABUW-RS13010  64.53336 8.794548e-13 2.601157e-11
#> PMB0.5.ABUW-RS07965 ABUW-RS07965  61.57093 1.392901e-12 3.967189e-11
#> PMB0.5.bauF                 bauF  81.50603 1.556705e-12 4.275378e-11
#> PMB0.5.basA                 basA  79.39009 1.869130e-12 4.956417e-11
#> PMB0.5.ABUW-RS14775 ABUW-RS14775  71.52750 2.888729e-12 7.404775e-11
#> PMB0.5.ABUW-RS10520 ABUW-RS10520  71.38293 3.013230e-12 7.474753e-11
#> PMB0.5.ABUW-RS11870 ABUW-RS11870  68.57500 4.823634e-12 1.159179e-10
#> PMB0.5.ABUW-RS06610 ABUW-RS06610  68.34097 5.091317e-12 1.172075e-10
#> PMB0.5.ABUW-RS05995 ABUW-RS05995  63.72752 5.182125e-12 1.172075e-10
#> PMB0.5.ABUW-RS10635 ABUW-RS10635  67.02443 6.486936e-12 1.425273e-10
#> PMB0.5.cydB                 cydB  61.80285 7.530302e-12 1.608556e-10
#> PMB0.5.ABUW-RS17680 ABUW-RS17680  67.93491 1.266994e-11 2.633293e-10
#> PMB0.5.bauE                 bauE  67.35122 1.357154e-11 2.711910e-10
#> PMB0.5.ABUW-RS10840 ABUW-RS10840  63.20971 1.375351e-11 2.711910e-10
#> PMB0.5.ABUW-RS20155 ABUW-RS20155  51.19627 1.530687e-11 2.893715e-10
#> PMB0.5.ABUW-RS10080 ABUW-RS10080  66.81576 1.542813e-11 2.893715e-10
#> PMB0.5.ABUW-RS10630 ABUW-RS10630  60.86984 2.199063e-11 4.026380e-10
#> PMB0.5.ABUW-RS07155 ABUW-RS07155  49.15796 2.539188e-11 4.381616e-10
#> PMB0.5.ABUW-RS07805 ABUW-RS07805  61.73398 2.543634e-11 4.381616e-10
#> PMB0.5.ABUW-RS10625 ABUW-RS10625  59.96293 2.564015e-11 4.381616e-10
#> PMB0.5.ABUW-RS17685 ABUW-RS17685  59.60557 2.776767e-11 4.642030e-10
#> PMB0.5.wecB.2             wecB.2  63.07824 3.466590e-11 5.671932e-10
#> PMB0.5.ABUW-RS14240 ABUW-RS14240  46.68910 4.425675e-11 6.982658e-10
#> PMB0.5.ABUW-RS18090 ABUW-RS18090  61.59561 4.449288e-11 6.982658e-10
#> PMB0.5.ABUW-RS07145 ABUW-RS07145  53.64626 5.506892e-11 8.297509e-10
#> PMB0.5.ABUW-RS09980 ABUW-RS09980  59.17171 5.542783e-11 8.297509e-10
#> PMB0.5.ABUW-RS18105 ABUW-RS18105  59.16288 5.610799e-11 8.297509e-10
#> PMB0.5.otsB                 otsB  57.41299 5.989823e-11 8.690894e-10
#> PMB0.5.ABUW-RS16850 ABUW-RS16850  56.82622 7.229196e-11 1.029491e-09
#> PMB0.5.ABUW-RS14485 ABUW-RS14485  53.56186 7.489752e-11 1.047204e-09
#> PMB0.5.bfr                   bfr  42.95712 8.359506e-11 1.147939e-09
#> PMB0.5.ABUW-RS08740 ABUW-RS08740  47.23525 9.360537e-11 1.262851e-09
#> PMB0.5.ABUW-RS05940 ABUW-RS05940  51.46687 9.758606e-11 1.293857e-09
#> PMB0.5.ABUW-RS00260 ABUW-RS00260  45.57234 1.279479e-10 1.667660e-09
#> PMB0.5.entE                 entE  51.85415 1.462736e-10 1.874740e-09
#> PMB0.5.ABUW-RS09970 ABUW-RS09970  55.05737 1.599485e-10 2.016400e-09
#> PMB0.5.bauB                 bauB  53.65730 1.875965e-10 2.310929e-09
#> PMB0.5.basH                 basH  51.69996 1.893219e-10 2.310929e-09
#> PMB0.5.tal                   tal  42.48045 2.201289e-10 2.644986e-09
#> PMB0.5.ABUW-RS10075 ABUW-RS10075  52.15178 2.884199e-10 3.412230e-09
#> PMB0.5.ABUW-RS11850 ABUW-RS11850  42.67572 3.194074e-10 3.721581e-09
#> PMB0.5.ABUW-RS20230 ABUW-RS20230  46.18665 3.319669e-10 3.810187e-09
#> PMB0.5.ABUW-RS15280 ABUW-RS15280  48.19657 3.477838e-10 3.933025e-09
#> PMB0.5.ABUW-RS16155 ABUW-RS16155  44.27340 3.671009e-10 4.091312e-09
#> PMB0.5.ABUW-RS09210 ABUW-RS09210  43.45243 4.055182e-10 4.454907e-09
#> PMB0.5.ABUW-RS10530 ABUW-RS10530  47.31161 6.513967e-10 7.055268e-09
#> PMB0.5.ABUW-RS12225 ABUW-RS12225  43.12068 7.046584e-10 7.526143e-09
#> PMB0.5.alr.1               alr.1  46.54509 7.184492e-10 7.568321e-09
#> PMB0.5.ABUW-RS16650 ABUW-RS16650  45.00855 7.749597e-10 8.053297e-09
#> PMB0.5.ABUW-RS05470 ABUW-RS05470  46.96099 8.574444e-10 8.791663e-09
#> PMB0.5.ABUW-RS15460 ABUW-RS15460  43.54298 1.103479e-09 1.116547e-08
#> PMB0.5.ABUW-RS18420 ABUW-RS18420  45.51632 1.253135e-09 1.251508e-08
#> PMB0.5.ABUW-RS07150 ABUW-RS07150  34.49763 1.466290e-09 1.443834e-08
#> PMB0.5.ABUW-RS10620 ABUW-RS10620  42.45268 1.493918e-09 1.443834e-08
#> PMB0.5.ABUW-RS03090 ABUW-RS03090  38.39377 1.502038e-09 1.443834e-08
#> PMB0.5.ABUW-RS13220 ABUW-RS13220  41.64028 1.622878e-09 1.540733e-08
#> PMB0.5.bfr.1               bfr.1  34.28916 1.865755e-09 1.749715e-08
#> PMB0.5.barA                 barA  41.56491 1.934985e-09 1.792775e-08
#> PMB0.5.ABUW-RS18410 ABUW-RS18410  37.51603 2.253200e-09 2.061600e-08
#> PMB0.5.ABUW-RS15755 ABUW-RS15755  42.18152 2.278752e-09 2.061600e-08
#> PMB0.5.ABUW-RS00905 ABUW-RS00905  42.41793 2.635887e-09 2.356974e-08
#> PMB0.5.ABUW-RS14475 ABUW-RS14475  37.49237 3.172291e-09 2.804014e-08
#> PMB0.5.ABUW-RS11945 ABUW-RS11945  37.49018 3.414240e-09 2.983580e-08
#> PMB0.5.bauC                 bauC  40.94877 3.660799e-09 3.163095e-08
#> PMB0.5.thiC                 thiC  38.91221 3.823057e-09 3.266590e-08
#> PMB0.5.lysM                 lysM  38.10559 4.402949e-09 3.720734e-08
#> PMB0.5.ABUW-RS03195 ABUW-RS03195  38.00515 5.006966e-09 4.185171e-08
#> PMB0.5.ABUW-RS10220 ABUW-RS10220  37.62022 5.601077e-09 4.631429e-08
#> PMB0.5.hemF                 hemF  35.34159 6.844039e-09 5.576799e-08
#> PMB0.5.ABUW-RS09965 ABUW-RS09965  36.94130 6.889414e-09 5.576799e-08
#> PMB0.5.carO                 carO  35.80543 9.567218e-09 7.663740e-08
#> PMB0.5.barB                 barB  34.94179 1.300064e-08 1.030669e-07
#> PMB0.5.ABUW-RS17035 ABUW-RS17035  30.75094 1.843651e-08 1.446702e-07
#> PMB0.5.ABUW-RS09185 ABUW-RS09185  30.88163 2.036715e-08 1.582055e-07
#> PMB0.5.ABUW-RS12760 ABUW-RS12760  33.97172 2.219101e-08 1.706489e-07
#> 
#> $PMB2
#>                           gene         F       PValue          FDR
#> PMB2.ABUW-RS13935 ABUW-RS13935 44.775158 3.758922e-09 6.120790e-07
#> PMB2.secY                 secY 42.376787 6.375823e-09 6.120790e-07
#> PMB2.ABUW-RS14775 ABUW-RS14775 39.858476 1.045649e-08 6.692154e-07
#> PMB2.ssrS                 ssrS 36.802239 5.793375e-08 2.780820e-06
#> PMB2.ABUW-RS04280 ABUW-RS04280 32.520384 8.452640e-08 3.245814e-06
#> PMB2.ABUW-RS05935 ABUW-RS05935 30.729018 1.403990e-07 3.960423e-06
#> PMB2.rplE                 rplE 29.992163 1.443904e-07 3.960423e-06
#> PMB2.rplN                 rplN 25.470144 4.612756e-07 1.107061e-05
#> PMB2.ABUW-RS05420 ABUW-RS05420 23.557994 5.711332e-07 1.218418e-05
#> PMB2.ABUW-RS03195 ABUW-RS03195 25.689886 6.769298e-07 1.299705e-05
#> PMB2.macA                 macA 23.539949 1.216195e-06 1.833229e-05
#> PMB2.cydB                 cydB 22.758573 1.231013e-06 1.833229e-05
#> PMB2.ABUW-RS10840 ABUW-RS10840 22.655581 1.280648e-06 1.833229e-05
#> PMB2.ABUW-RS01745 ABUW-RS01745 22.056205 1.336729e-06 1.833229e-05
#> PMB2.ABUW-RS17530 ABUW-RS17530 24.759050 1.632430e-06 2.089511e-05
#> PMB2.rplB                 rplB 22.816092 1.867232e-06 2.223042e-05
#> PMB2.nuoL                 nuoL 20.144743 1.968318e-06 2.223042e-05
#> PMB2.fusA                 fusA 21.832622 2.677396e-06 2.581224e-05
#> PMB2.ABUW-RS00825 ABUW-RS00825 21.825250 2.683939e-06 2.581224e-05
#> PMB2.rpsH                 rpsH 21.678642 2.688775e-06 2.581224e-05
#> PMB2.ABUW-RS19635 ABUW-RS19635 21.597586 4.230755e-06 3.868119e-05
#> PMB2.lolA                 lolA 19.278918 5.331765e-06 4.569590e-05
#> PMB2.ABUW-RS08740 ABUW-RS08740 19.467568 5.473988e-06 4.569590e-05
#> PMB2.ABUW-RS17440 ABUW-RS17440 18.919186 5.893135e-06 4.714508e-05
#> PMB2.ABUW-RS04445 ABUW-RS04445 20.067725 8.240414e-06 6.328638e-05
#> PMB2.ABUW-RS17525 ABUW-RS17525 18.429489 1.032125e-05 7.621850e-05
#> PMB2.rplX                 rplX 17.716379 1.132983e-05 8.056767e-05
#> PMB2.tuf                   tuf 17.744886 1.333121e-05 9.141402e-05
#> PMB2.atpB                 atpB 17.022134 1.537553e-05 1.017966e-04
#> PMB2.ffs                   ffs 17.590841 1.825198e-05 1.168126e-04
#> PMB2.ABUW-RS13430 ABUW-RS13430 17.748080 2.002369e-05 1.240177e-04
#> PMB2.ABUW-RS03410 ABUW-RS03410 17.586157 2.149346e-05 1.289607e-04
#> PMB2.ABUW-RS17255 ABUW-RS17255 16.453440 2.372388e-05 1.374609e-04
#> PMB2.rpsP                 rpsP 15.848215 2.434203e-05 1.374609e-04
#> PMB2.rplW                 rplW 16.073813 2.515216e-05 1.379775e-04
#> PMB2.adeI                 adeI 15.176351 2.681463e-05 1.430113e-04
#> PMB2.carO                 carO 15.064832 2.762329e-05 1.433425e-04
#> PMB2.ABUW-RS15630 ABUW-RS15630 16.362697 2.869037e-05 1.449619e-04
#> PMB2.ABUW-RS17510 ABUW-RS17510 16.001052 2.977880e-05 1.466033e-04
#> PMB2.rpsK                 rpsK 15.845099 3.139754e-05 1.472267e-04
#> PMB2.ABUW-RS15625 ABUW-RS15625 16.656360 3.143902e-05 1.472267e-04
#> PMB2.rpsF                 rpsF 15.206136 4.232574e-05 1.934891e-04
#> PMB2.rpsE                 rpsE 14.910846 4.571980e-05 2.041442e-04
#> PMB2.rpsM                 rpsM 14.650771 5.491779e-05 2.396412e-04
#> PMB2.odhB                 odhB 14.195971 5.900712e-05 2.517637e-04
#> PMB2.ABUW-RS02985 ABUW-RS02985 13.069594 6.650788e-05 2.775981e-04
#> PMB2.cyoA                 cyoA 14.092466 7.164601e-05 2.926816e-04
#> PMB2.ABUW-RS18340 ABUW-RS18340 14.447573 7.828661e-05 3.131464e-04
#> PMB2.sucC                 sucC 13.577948 8.462825e-05 3.316046e-04
#> PMB2.ABUW-RS05820 ABUW-RS05820 13.057907 9.839722e-05 3.778453e-04
#> PMB2.lpdA                 lpdA 13.003882 1.019032e-04 3.836354e-04
#> PMB2.rplK                 rplK 13.069533 1.180552e-04 4.358960e-04
#> PMB2.ABUW-RS18750 ABUW-RS18750 12.339233 1.261773e-04 4.501474e-04
#> PMB2.rpsD                 rpsD 12.920317 1.266040e-04 4.501474e-04
#> PMB2.ABUW-RS01895 ABUW-RS01895 12.300452 1.442983e-04 5.037323e-04
#> PMB2.rpoA                 rpoA 12.513677 1.545097e-04 5.297477e-04
#> PMB2.ABUW-RS05265 ABUW-RS05265 12.458803 1.610215e-04 5.404339e-04
#> PMB2.rpsU                 rpsU 11.427817 1.632561e-04 5.404339e-04
#> PMB2.rplC                 rplC 12.220206 1.798463e-04 5.852625e-04
#> PMB2.ABUW-RS06040 ABUW-RS06040 11.867717 1.915597e-04 6.129911e-04
#> PMB2.rplF                 rplF 11.949993 2.088967e-04 6.575109e-04
#> PMB2.rpmC                 rpmC 11.796905 2.260551e-04 7.000416e-04
#> PMB2.ABUW-RS02640 ABUW-RS02640 12.044884 2.526893e-04 7.624344e-04
#> PMB2.ABUW-RS15035 ABUW-RS15035 11.564867 2.541448e-04 7.624344e-04
#> PMB2.gltA                 gltA 11.467439 2.657105e-04 7.848679e-04
#> PMB2.atpE                 atpE 11.003649 2.720342e-04 7.913722e-04
#> PMB2.macB                 macB 11.401186 2.799413e-04 8.022198e-04
#> PMB2.ABUW-RS00170 ABUW-RS00170 10.817572 3.260436e-04 9.205937e-04
#> PMB2.minD                 minD 10.315456 3.673101e-04 1.022080e-03
#> PMB2.ABUW-RS04410 ABUW-RS04410 10.721148 3.841943e-04 1.053790e-03
#> PMB2.rpsL                 rpsL 10.639840 4.145587e-04 1.121060e-03
#> PMB2.ABUW-RS19280 ABUW-RS19280 10.453616 4.616099e-04 1.230960e-03
#> PMB2.ABUW-RS04655 ABUW-RS04655  9.725592 5.827476e-04 1.532706e-03
#> PMB2.rpoB                 rpoB  9.844415 6.558314e-04 1.701617e-03
#> PMB2.ABUW-RS19720 ABUW-RS19720 10.152212 6.680798e-04 1.710284e-03
#> PMB2.trmD                 trmD  9.495893 6.955336e-04 1.757137e-03
#> PMB2.rplS                 rplS  9.563148 7.322857e-04 1.825959e-03
#> PMB2.rplD                 rplD  9.400980 8.454739e-04 2.081166e-03
#> PMB2.ABUW-RS07780 ABUW-RS07780  9.064400 1.024584e-03 2.490127e-03
#> PMB2.rpsC                 rpsC  9.010341 1.057524e-03 2.538058e-03
#> PMB2.omp33-36         omp33-36  8.982888 1.080048e-03 2.560114e-03
#> PMB2.ABUW-RS19190 ABUW-RS19190  9.097290 1.104279e-03 2.585629e-03
#> PMB2.ABUW-RS02185 ABUW-RS02185  9.147778 1.151417e-03 2.656883e-03
#> PMB2.rpmB                 rpmB  8.847770 1.162386e-03 2.656883e-03
#> PMB2.bfmR                 bfmR  8.515281 1.189835e-03 2.687626e-03
#> PMB2.ABUW-RS01250 ABUW-RS01250  8.600730 1.349849e-03 3.008719e-03
#> PMB2.adeJ                 adeJ  8.596272 1.363326e-03 3.008719e-03
#> PMB2.ABUW-RS04045 ABUW-RS04045  8.308487 1.384655e-03 3.021065e-03
#> PMB2.ahpC                 ahpC  8.533211 1.418885e-03 3.060966e-03
#> PMB2.rpsJ                 rpsJ  8.456379 1.483430e-03 3.135358e-03
#> PMB2.ftsZ                 ftsZ  8.307926 1.486029e-03 3.135358e-03
#> PMB2.rplR                 rplR  8.225129 1.667455e-03 3.479906e-03
#> PMB2.rpsG                 rpsG  8.166122 1.770269e-03 3.654749e-03
#> PMB2.thrS                 thrS  7.865336 1.954906e-03 3.959214e-03
#> PMB2.putA                 putA  7.925730 1.958986e-03 3.959214e-03
#> PMB2.nuoG                 nuoG  7.937747 2.039388e-03 4.078776e-03
#> PMB2.nuoM                 nuoM  7.769664 2.214093e-03 4.382535e-03
#> PMB2.ABUW-RS02255 ABUW-RS02255  7.412943 2.253090e-03 4.414217e-03
#> PMB2.ABUW-RS05940 ABUW-RS05940  7.439262 2.518358e-03 4.884089e-03
#> PMB2.atpG                 atpG  7.487426 2.665131e-03 5.117051e-03
```

Plot fitted curves for genes with strong time-dependent changes:

``` r
curve_genes <- unique(unlist(lapply(top_by_condition, function(df) df$gene[seq_len(min(6, nrow(df)))])))

plot_time_series_deg_curves(
  spline_de,
  genes = curve_genes,
  ncol = 4
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-32-1.png)

For users who prefer a one-step interface,
[`run_scproka_de()`](https://monash-li-lab.github.io/scProka/reference/run_scproka_de.md)
wraps the pseudobulk and edgeR steps into one function call.

## Time-aware trajectory analysis with TATA

`SCProkaR` also includes Time Aware Trajectory Analysis (TATA), a
PAGA-inspired trajectory abstraction method that uses the observed
experimental time labels when scoring and directing cluster-level
transitions.

For the practical workflow here, we continue with the same integrated
object and run TATA on the batch-corrected latent space. This gives a
compact cluster-level trajectory graph and a pseudotime estimate that
can be interpreted on the integrated UMAP.

### Run TATA on the integrated embedding

We use the MNN-integrated embedding for graph construction and the
`timepoints` column as the real experimental time variable.

``` r
tata <- run_tata(
  sce = sce,
  dimred = "integrated_mnn",
  time_col = "timepoints",
  cluster_col = "tata_cluster",
  k = 20,
  cluster_method = "louvain",
  alpha = 1,
  beta = 1,
  direction_threshold = 0.10,
  prune_threshold = 0.001,
  time_weight_mode = "directional_confidence",
  do_pseudotime = TRUE,
  seed = 101
)

sce_tata <- tata$sce
sce_tata$tata_cluster <- factor(sce_tata$tata_cluster)
```

The result stores the updated `SingleCellExperiment`, the cell-cell
graph, the abstract cluster graph, the cluster edge table, and the
cell-level pseudotime.

``` r
names(tata)
#>  [1] "sce"                    "cell_graph"             "adjacency_matrix"      
#>  [4] "cluster_graph"          "cluster_edge_table"     "cluster_vertex_table"  
#>  [7] "cluster_pseudotime"     "cell_pseudotime"        "cell_pseudotime_scaled"
#> [10] "parameters"             "topology_table"         "time_table"

head(
  tata$cluster_edge_table[
    order(-tata$cluster_edge_table$tata_weight),
    c(
      "cluster_a",
      "cluster_b",
      "topology_weight",
      "flow_score",
      "time_weight",
      "tata_weight",
      "direction_label",
      "kept"
    )
  ],
  10
)
#>    cluster_a cluster_b topology_weight  flow_score time_weight tata_weight
#> 10        C4        C5     0.717792713 -0.45394737   0.2730263 0.521816413
#> 7         C2        C5     0.664168441  0.47860963   0.7393048 0.491022925
#> 8         C3        C4     0.578898891  0.05062082   0.5253104 0.304101614
#> 9         C3        C5     0.468419780 -0.29692557   0.3515372 0.303752794
#> 5         C2        C3     0.185994078  0.36193772   0.6809689 0.126656175
#> 2         C1        C3     0.090666514  0.39957265   0.6997863 0.063447187
#> 6         C2        C4     0.045895442  0.51075269   0.7553763 0.034668331
#> 1         C1        C2     0.056543071 -0.17889908   0.4105505 0.033329287
#> 4         C1        C5     0.045438643  0.12605042   0.5630252 0.025583101
#> 3         C1        C4     0.005385108  0.20689655   0.6034483 0.003249634
#>    direction_label kept
#> 10        C5 -> C4 TRUE
#> 7         C2 -> C5 TRUE
#> 8         C3 -- C4 TRUE
#> 9         C5 -> C3 TRUE
#> 5         C2 -> C3 TRUE
#> 2         C1 -> C3 TRUE
#> 6         C2 -> C4 TRUE
#> 1         C2 -> C1 TRUE
#> 4         C1 -> C5 TRUE
#> 3         C1 -> C4 TRUE
```

### Visualize TATA clusters and pseudotime

``` r
PlotReduction(sce_tata, reduction = "umap_integrated_mnn", colour_by = "tata_cluster")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-35-1.png)

``` r
PlotReduction(sce_tata, reduction = "umap_integrated_mnn", colour_by = "tata_pseudotime_scaled")
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-35-2.png)

### Overlay the abstract TATA graph on the embedding

[`plot_tata_cluster_graph()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_cluster_graph.md)
places the cluster-level graph directly on top of the embedding so the
abstract transitions can be interpreted in the same coordinate system as
the cells.

``` r
plot_tata_cluster_graph(
  tata,
  dimred = "umap_integrated_mnn",
  layout_mode = "embedding",
  show_cells = TRUE,
  colour_by = "timepoints",
  node_colour_by = "cluster",
  main = "TATA abstract graph aligned to the integrated UMAP"
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-36-1.png)

### Plot the inferred trajectories

For a cleaner trajectory-style summary,
[`plot_tata_trajectory_embedding()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_trajectory_embedding.md)
draws fitted paths on top of the embedding. The `combined` mode is
especially useful in a workflow vignette because it compresses the graph
into a smaller set of representative trajectories.

``` r
plot_tata_trajectory_embedding(
  tata,
  dimred = "umap_integrated_mnn",
  colour_by = "timepoints",
  mode = "combined",
  max_paths = 3,
  smooth_paths = TRUE
)
```

![](scprokar-workflow_files/figure-html/unnamed-chunk-37-1.png)

This practical example is meant to show how TATA fits into the main
`SCProkaR` workflow on a real dataset. For the full step-by-step TATA
method walkthrough, simulated multi-branch examples, and benchmarking
against other trajectory tools, see the separate TATA vignette.

## Notes

- [`CreateBacObject()`](https://monash-li-lab.github.io/scProka/reference/CreateBacObject.md)
  accepts Seurat input directly, but the downstream analysis object
  remains a `SingleCellExperiment`.
- [`MergeBacObjects()`](https://monash-li-lab.github.io/scProka/reference/MergeBacObjects.md)
  is the recommended merge path for multiple samples in `SCProkaR`.
- [`PlotReduction()`](https://monash-li-lab.github.io/scProka/reference/PlotReduction.md)
  is a convenient visualization helper for both original and integrated
  reductions.
- [`PlotIntegrationOverview()`](https://monash-li-lab.github.io/scProka/reference/PlotIntegrationOverview.md)
  provides the standard Seurat-style integration comparison panels.
- [`RunIntegratedClustering()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedClustering.md)
  provides the post-integration cluster labels for future downstream
  analyses.
- [`RunIntegratedUMAP()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedUMAP.md)
  works on integrated embeddings as well as any custom latent space
  registered with
  [`RegisterIntegrationEmbedding()`](https://monash-li-lab.github.io/scProka/reference/RegisterIntegrationEmbedding.md).
- Pseudobulk differential expression is covered in the same package
  workflow via
  [`aggregate_pseudobulk()`](https://monash-li-lab.github.io/scProka/reference/aggregate_pseudobulk.md),
  [`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md),
  and
  [`run_edger_spline_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_spline_de.md).
- TATA can be run directly on the integrated object with
  [`run_tata()`](https://monash-li-lab.github.io/scProka/reference/run_tata.md)
  and visualized with
  [`plot_tata_cluster_graph()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_cluster_graph.md)
  and
  [`plot_tata_trajectory_embedding()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_trajectory_embedding.md).
- The separate TATA vignette focuses on the step-wise method
  explanation, simulated data, and benchmarking.
