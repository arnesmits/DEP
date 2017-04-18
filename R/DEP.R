#' DEP: A package for Differential Enrichment analysis of Proteomics data.
#'
#' The DEP package provides all functionalities to analyze label free protein quantification (LFQ) and tandem-mass-tags (TMT) labeled proteomics data
#'
#' @section Shiny apps:
#' \itemize{
#'   \item \code{\link{run_app}}: Shiny apps for interactive analysis
#' }
#'
#' @section Workflows:
#' \itemize{
#'   \item \code{\link{LFQ}}: Label-free quantification (LFQ) workflow
#'   \item \code{\link{TMT}}: Tandem-mass-tags (TMT) workflow
#'   \item \code{\link{report}}: Create a rmarkdown report
#'   \item \code{\link{iBAQ}}: iBAQ-based stoichiometry workflow
#' }
#'
#' @section Functions:
#' \itemize{
#'   \item \code{\link{unique_names}}: Generate unique names
#'   \item \code{\link{make_se_parse}}: Turn data.frame into SummarizedExperiment by parsing column names
#'   \item \code{\link{make_se}}: Turn data.frame into SummarizedExperiment using an experimental design
#'   \item \code{\link{filter_missval}}: Filter on missing values
#'   \item \code{\link{norm_vsn}}: Normalize data using vsn
#'   \item \code{\link{imputation}}: Impute missing values
#'   \item \code{\link{linear_model}}: Differential expression analysis
#'   \item \code{\link{cutoffs}}: Denote significant proteins
#'   \item \code{\link{results}}: Generate a results table
#' }
#'
#' @section Visualization functions:
#' \itemize{
#'   \item \code{\link{plot_single}}: Barplot for a protein of interest
#'   \item \code{\link{plot_volcano}}: Volcano plot for a specified contrast
#'   \item \code{\link{plot_heatmap}}: Heatmap of all significant proteins
#'   \item \code{\link{plot_ibaq}}: Scatter plot of iBAQ intensities versus log fold changes for a specified contrast
#'   \item \code{\link{plot_norm}}: Boxplots to inspect normalization
#'   \item \code{\link{plot_detect}}: Density and CumSum plots to inspect detection limit
#'   \item \code{\link{plot_imp}}: Density plots to inspect imputation
#'   \item \code{\link{plot_missval}}: Heatmap to inspect missing values
#'   \item \code{\link{plot_numbers}}: Barplot of proteins identified
#'   \item \code{\link{plot_frequency}}: Barplot of protein identification overlap between conditions
#'   \item \code{\link{plot_coverage}}: Barplot of the protein coverage in conditions
#' }
#'
#' @section iBAQ functions:
#' \itemize{
#'   \item \code{\link{ibaq_merge}}: Merge iBAQ intensities of protein groups based on shared peptides
#'   \item \code{\link{stoichiometry}}: Calculate relative stoichiometry of all significant proteins
#'   \item \code{\link{plot_stoi}}: Barplot of the relative stoichiometries
#' }
#'
#' @docType package
#' @name DEP
#'
#'
#' @import ComplexHeatmap limma magrittr tidyverse grid SummarizedExperiment assertthat
#' @import shinydashboard shiny
#'
NULL
