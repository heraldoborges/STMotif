#' STMotif: Discovery and Ranking of Motifs in Spatial-Time Series
#'
#' Identifies motifs in spatial-time series using the Combined Series
#' Approach (CSA). A motif is a previously unknown subsequence of a spatial
#' time series with a relevant number of occurrences. The package provides
#' functions for SAX encoding, motif discovery, ranking, and visualization.
#'
#' The main workflow consists of:
#' \enumerate{
#'   \item \code{\link{NormSAX}}: Normalize and encode the dataset with SAX.
#'   \item \code{\link{SearchSTMotifs}}: Discover motifs in spatio-temporal blocks.
#'   \item \code{\link{RankSTMotifs}}: Rank motifs by quality metrics.
#' }
#'
#' These three steps are wrapped in the convenience function
#' \code{\link{CSAMiningProcess}}.
#'
#' Visualization functions:
#' \itemize{
#'   \item \code{\link{display_motifsDataset}}: Heatmap with highlighted motifs.
#'   \item \code{\link{display_motifsSTSeries}}: Time series plots with highlighted motifs.
#' }
#'
#' @seealso \code{vignette("STMotif")} for a comprehensive introduction.
#'
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTE for NSE variable used in facet formula
utils::globalVariables("variable")
