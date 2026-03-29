#' Convert wide data frame to long format
#'
#' Replacement for \code{reshape2::melt(df, id.vars = 1)}.
#' Assumes the first column is the id variable.
#'
#' @param df Data frame with id column in position 1.
#' @return Data frame with columns \code{x}, \code{variable}, \code{value}.
#' @keywords internal
#' @noRd
.wide_to_long <- function(df) {
  id_col <- df[[1]]
  measure_names <- colnames(df)[-1]
  n_measures <- length(measure_names)
  n_rows <- nrow(df)
  
  data.frame(
    x        = rep(id_col, times = n_measures),
    variable = factor(rep(measure_names, each = n_rows), levels = measure_names),
    value    = unlist(df[, -1, drop = FALSE], use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


#' Convert matrix to long format
#'
#' Replacement for \code{reshape2::melt(mat)} applied to a matrix.
#' Produces columns Var1 (row index), Var2 (column index), value.
#'
#' @param mat A numeric matrix.
#' @return Data frame with columns \code{Var1}, \code{Var2}, \code{value}.
#' @keywords internal
#' @noRd
.melt_matrix <- function(mat) {
  data.frame(
    Var1  = as.vector(row(mat)),
    Var2  = as.vector(col(mat)),
    value = as.vector(mat)
  )
}


#' Generate a safe color palette
#'
#' Wrapper around \code{RColorBrewer::brewer.pal()} that handles
#' \code{n < 3} (the minimum required by brewer.pal).
#'
#' @param n Number of colors needed.
#' @param palette Palette name (default: "Spectral").
#' @return Character vector with \code{n} hex color strings.
#' @keywords internal
#' @noRd
.safe_palette <- function(n, palette = "Spectral") {
  if (n < 1) return(character(0))
  n_pal <- max(n, 3)
  cores <- RColorBrewer::brewer.pal(n_pal, palette)
  cores[seq_len(n)]
}


#' Plot spatial-time series with colors
#'
#' Internal rendering function. Accepts a long-format data frame with
#' columns: x, variable, value, color.
#'
#' @param series Data frame in long format.
#' @param label_x X axis label.
#' @param label_y Y axis label.
#' @return A ggplot object.
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes
#' @keywords internal
#' @noRd
.plot_series <- function(series, label_x = "", label_y = "") {
  ggplot2::ggplot(
    data = series,
    ggplot2::aes(x = .data$x, y = .data$value,
                 colour = .data$color, group = 1)
  ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_grid(variable ~ .) +
    ggplot2::xlab(label_x) +
    ggplot2::ylab(label_y) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_blank()
    )
}


# --- Exported functions -----------------------------------------------------

#' Plot Spatial-Time Series with Highlighted Motifs
#'
#' Displays the selected spatial-time series and highlights the segments
#' corresponding to the discovered motifs using distinct colors.
#'
#' @param dataset Data frame or matrix containing numeric values.
#'   Each column represents a spatial series, each row a time point.
#' @param rstmotifs List of ranked motifs, as returned by
#'   \code{\link{RankSTMotifs}} or \code{\link{CSAMiningProcess}}.
#' @param space Integer vector specifying which columns (spatial series)
#'   to display. Defaults to all columns.
#' @return A \code{\link[ggplot2]{ggplot}} object showing the time series
#'   with motif occurrences highlighted in color.
#'
#' @importFrom rlang .data
#'
#' @examples
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' display_motifsSTSeries(
#'   dataset = STMotif::example_dataset,
#'   rstmotifs[c(1:4)],
#'   space = c(1:4, 10:12)
#' )
#'
#' @export
display_motifsSTSeries <- function(dataset, rstmotifs,
                                   space = seq_len(ncol(dataset))) {
  dataset <- as.data.frame(dataset)
  n_series <- ncol(dataset)
  colnames(dataset) <- as.character(seq_len(n_series))
  
  # Series labels: "ST1", "ST2", ...
  series_labels <- paste0("ST", colnames(dataset))
  
  # Select columns and convert to long format
  selected <- as.data.frame(dataset[, space, drop = FALSE])
  colnames(selected) <- series_labels[space]
  selected <- data.frame(x = seq_len(nrow(selected)), selected,
                         check.names = FALSE)
  long_data <- .wide_to_long(selected)
  long_data$color <- "black"
  
  # Color motif segments
  if (!is.null(rstmotifs) && length(rstmotifs) > 0) {
    palette <- .safe_palette(length(rstmotifs))
    motif_size <- nchar(rstmotifs[[1]]$isaxcod)
    
    for (pos in seq_along(rstmotifs)) {
      motif <- rstmotifs[[pos]]
      for (j in seq_along(motif$vecst$s)) {
        s_idx <- motif$vecst$s[j]
        t_idx <- motif$vecst$t[j]
        
        if (s_idx %in% space) {
          target_label <- series_labels[s_idx]
          time_range <- t_idx:(t_idx + motif_size - 1)
          
          mask <- long_data$variable == target_label &
            long_data$x %in% time_range
          long_data$color[mask] <- palette[pos]
        }
      }
    }
  }
  
  .plot_series(long_data)
}


#' Plot Heatmap with Highlighted Motifs
#'
#' Displays the dataset as a heatmap (encoded via SAX binning) and overlays
#' colored markers at the positions where the selected motifs occur.
#'
#' @param dataset Data frame or matrix containing numeric values.
#'   Each column represents a spatial series, each row a time point.
#' @param rstmotifs List of ranked motifs, as returned by
#'   \code{\link{RankSTMotifs}} or \code{\link{CSAMiningProcess}}.
#' @param alpha Integer. The cardinality of the SAX alphabet (number of
#'   discretization levels).
#' @return A \code{\link[ggplot2]{ggplot}} object showing the heatmap with
#'   motif positions highlighted as colored squares.
#'
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices grey.colors
#'
#' @examples
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' display_motifsDataset(
#'   dataset = STMotif::example_dataset,
#'   rstmotifs[c(1:4)],
#'   5
#' )
#'
#' @export
display_motifsDataset <- function(dataset, rstmotifs, alpha) {
  
  # 1. SAX encoding -> integer values for grey scale
  color_encode <- seq_len(alpha)
  vec <- as.vector(as.matrix(dataset))
  vec_norm <- STSNormalization(vec)
  my_bin <- binning(vec_norm, alpha)
  encoded <- color_encode[my_bin$bins_factor]
  
  # Reshape to space x time matrix (transposed from original layout)
  mat_encoded <- t(matrix(encoded, nrow = nrow(dataset), ncol = ncol(dataset)))
  
  # Convert to long format
  heatmap_data <- .melt_matrix(mat_encoded)
  heatmap_data$motif <- FALSE
  heatmap_data$color <- NA_character_
  
  # 2. Mark motif positions
  if (!is.null(rstmotifs) && length(rstmotifs) > 0) {
    palette <- .safe_palette(length(rstmotifs))
    
    # Build motif positions table
    motif_positions <- do.call(rbind, lapply(seq_along(rstmotifs), function(pos) {
      data.frame(
        s     = rstmotifs[[pos]]$vecst$s,
        t     = rstmotifs[[pos]]$vecst$t,
        color = palette[pos],
        stringsAsFactors = FALSE
      )
    }))
    
    # Left join to identify which cells have motifs
    merged <- merge(
      heatmap_data, motif_positions,
      by.x = c("Var1", "Var2"), by.y = c("s", "t"),
      all.x = TRUE, suffixes = c("", ".motif")
    )
    
    has_motif <- !is.na(merged$color.motif)
    merged$motif <- has_motif
    merged$color[has_motif] <- merged$color.motif[has_motif]
    merged$color.motif <- NULL
    
    heatmap_data <- merged
  }
  
  # 3. Render
  point_colours <- ifelse(heatmap_data$motif, heatmap_data$color, NA)
  
  ggplot2::ggplot(
    data = heatmap_data,
    ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$value)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = c("white", "dimgrey"),
      values  = scales::rescale(seq_len(alpha)),
      limits  = c(1, alpha)
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("") +
    ggplot2::xlab("Space") +
    ggplot2::ylab("Time") +
    ggplot2::scale_y_reverse() +
    ggplot2::guides(fill = "none") +
    ggplot2::geom_point(
      colour      = point_colours,
      size        = 4,
      shape       = 15,
      show.legend = FALSE
    )
}
