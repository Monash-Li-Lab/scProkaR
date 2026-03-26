# Package index

## Data and object creation

- [`sce1`](https://monash-li-lab.github.io/scProka/reference/sce1.md) :
  Example bacterial single-cell dataset
- [`CreateBacObject()`](https://monash-li-lab.github.io/scProka/reference/CreateBacObject.md)
  : Create a standardized bacterial single-cell object
- [`MergeBacObjects()`](https://monash-li-lab.github.io/scProka/reference/MergeBacObjects.md)
  : Merge multiple SCProkaR or SingleCellExperiment objects

## QC, integration, and visualization

- [`RunBacQC()`](https://monash-li-lab.github.io/scProka/reference/RunBacQC.md)
  : Compute bacterial QC metrics
- [`FilterBacCells()`](https://monash-li-lab.github.io/scProka/reference/FilterBacCells.md)
  : Filter cells using bacterial QC thresholds
- [`IntegrateBacData()`](https://monash-li-lab.github.io/scProka/reference/IntegrateBacData.md)
  : Integrate microbial single-cell batches
- [`RegisterIntegrationEmbedding()`](https://monash-li-lab.github.io/scProka/reference/RegisterIntegrationEmbedding.md)
  : Register an external embedding for benchmarking and downstream
  analysis
- [`RunIntegratedUMAP()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedUMAP.md)
  : Run UMAP on an existing embedding
- [`RunIntegratedClustering()`](https://monash-li-lab.github.io/scProka/reference/RunIntegratedClustering.md)
  : Cluster cells on an integrated embedding
- [`PlotReduction()`](https://monash-li-lab.github.io/scProka/reference/PlotReduction.md)
  : Plot any stored reduced dimension
- [`PlotIntegrationOverview()`](https://monash-li-lab.github.io/scProka/reference/PlotIntegrationOverview.md)
  : Plot Seurat-style integration result panels

## Differential expression

- [`aggregate_pseudobulk()`](https://monash-li-lab.github.io/scProka/reference/aggregate_pseudobulk.md)
  : Aggregate single-cell data into pseudobulk samples
- [`filter_pseudobulk_samples()`](https://monash-li-lab.github.io/scProka/reference/filter_pseudobulk_samples.md)
  : Filter pseudobulk samples by size
- [`normalize_pseudobulk()`](https://monash-li-lab.github.io/scProka/reference/normalize_pseudobulk.md)
  : Normalize pseudobulk samples for exploratory analysis
- [`run_edger_pairwise_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_pairwise_de.md)
  : Run pairwise pseudobulk differential expression with edgeR
- [`run_edger_spline_de()`](https://monash-li-lab.github.io/scProka/reference/run_edger_spline_de.md)
  : Run pseudobulk spline or time-series differential expression with
  edgeR
- [`run_scproka_de()`](https://monash-li-lab.github.io/scProka/reference/run_scproka_de.md)
  : Run a full SCProkaR differential expression workflow
- [`plot_pairwise_de_volcano()`](https://monash-li-lab.github.io/scProka/reference/plot_pairwise_de_volcano.md)
  : Plot a volcano plot for a pairwise edgeR contrast
- [`plot_pairwise_de_ma()`](https://monash-li-lab.github.io/scProka/reference/plot_pairwise_de_ma.md)
  : Plot an MA plot for a pairwise edgeR contrast
- [`plot_time_series_deg_curves()`](https://monash-li-lab.github.io/scProka/reference/plot_time_series_deg_curves.md)
  : Plot fitted spline expression curves from time-series DE

## TATA trajectory analysis

- [`simulate_tata_sce()`](https://monash-li-lab.github.io/scProka/reference/simulate_tata_sce.md)
  [`simulate_tata_divergent_sce()`](https://monash-li-lab.github.io/scProka/reference/simulate_tata_sce.md)
  : Simulate a divergent TATA dataset
- [`simulate_tata_multidrug_sce()`](https://monash-li-lab.github.io/scProka/reference/simulate_tata_multidrug_sce.md)
  : Simulate a multi-drug branching TATA dataset
- [`simulate_tata_spiral_sce()`](https://monash-li-lab.github.io/scProka/reference/simulate_tata_spiral_sce.md)
  : Simulate an extreme spiral-plus-divergent TATA dataset
- [`run_tata()`](https://monash-li-lab.github.io/scProka/reference/run_tata.md)
  : Run the full TATA workflow
- [`build_knn_graph()`](https://monash-li-lab.github.io/scProka/reference/build_knn_graph.md)
  : Build a cell-cell kNN graph
- [`cluster_graph_states()`](https://monash-li-lab.github.io/scProka/reference/cluster_graph_states.md)
  : Cluster graph-connected cell states
- [`compute_topology_weights()`](https://monash-li-lab.github.io/scProka/reference/compute_topology_weights.md)
  : Compute topology-based cluster connectivity weights
- [`compute_time_weights()`](https://monash-li-lab.github.io/scProka/reference/compute_time_weights.md)
  : Compute time consistency for cluster-cluster edges
- [`build_tata_graph()`](https://monash-li-lab.github.io/scProka/reference/build_tata_graph.md)
  : Build the TATA cluster graph
- [`compute_tata_pseudotime()`](https://monash-li-lab.github.io/scProka/reference/compute_tata_pseudotime.md)
  : Compute TATA pseudotime
- [`plot_tata_cluster_graph()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_cluster_graph.md)
  : Plot the cluster-level TATA graph
- [`plot_tata_trajectory_embedding()`](https://monash-li-lab.github.io/scProka/reference/plot_tata_trajectory_embedding.md)
  : Plot TATA trajectories over an embedding

## Benchmarking

- [`BenchmarkIntegration()`](https://monash-li-lab.github.io/scProka/reference/BenchmarkIntegration.md)
  : Benchmark integration quality with scIB-inspired metrics
- [`orient_pseudotime_to_time()`](https://monash-li-lab.github.io/scProka/reference/orient_pseudotime_to_time.md)
  : Orient pseudotime to match increasing real time
- [`local_temporal_order_metrics()`](https://monash-li-lab.github.io/scProka/reference/local_temporal_order_metrics.md)
  : Compute local temporal order metrics on a cell graph
- [`adjacent_timepoint_auc()`](https://monash-li-lab.github.io/scProka/reference/adjacent_timepoint_auc.md)
  : Compute adjacent-timepoint AUC for a pseudotime ordering
- [`branch_onset_metrics()`](https://monash-li-lab.github.io/scProka/reference/branch_onset_metrics.md)
  : Measure branch-emergence timing from pseudotime
- [`directed_reachability_concordance()`](https://monash-li-lab.github.io/scProka/reference/directed_reachability_concordance.md)
  : Quantify how well a directed cluster graph follows real time
- [`anticausal_edge_mass()`](https://monash-li-lab.github.io/scProka/reference/anticausal_edge_mass.md)
  : Measure the weight of anti-causal edges in a directed graph
