#' Plot the selected spatial-time series with the selected motifs highlighted
#'
#' @param dataset Dataset containing numeric values
#' @param rmotifs List of ranked motifs
#' @param space Select a range of columns to plot the corresponding spatial series
#' @return Selected spatial series with the selected motifs highlighted
#' @import ggplot2
#' @import reshape2
#' @examples
#' #Launch all the workflow
#' #Plot the result
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,5)
#' stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' display_motifsSTSeries(dataset = STMotif::example_dataset,rstmotifs[c(1:4)],space = c(1:4,10:12))
#' @export
#'
display_motifsSTSeries <- function (dataset, rmotifs,space = c(1:length(dataset))){
  dataset <- as.data.frame(dataset)
  colnames(dataset) <- paste("",1:length(dataset), sep = "")

  size_motif <- nchar(rmotifs[[1]]$isaxcod)
  namesCol <- paste("ST",colnames(dataset),sep = "")
  data <- as.data.frame(dataset[,space])
  colnames(data) <- paste("ST",colnames(dataset)[space], sep = "")
  data <- data.frame(x = 1:nrow(data),data)
  data <- reshape2::melt(data,id.vars = 1)
  data <- data.frame(data, color = "black")
  palhetaCores <- brewer.pal(length(rmotifs), 'Spectral')
  levels(data$color) <- c("black", palhetaCores)

  for (position in 1:length(rmotifs)){
    for(i in 1:length(rmotifs[[position]]$vecst$s)){
      if(rmotifs[[position]]$vecst$s[i]%in%space){
        data[data$variable==namesCol[rmotifs[[position]]$vecst$s[i]],][(rmotifs[[position]]$vecst$t[i]):(rmotifs[[position]]$vecst$t[i]+(size_motif-1)),4] <- palhetaCores[position]
      }
    }
  }

  plot.series(data[1:nrow(data),])

}


#' Plot a heatmap of the dataset and highlight the selected motifs from the list
#'
#' @param dataset Numerical dataset
#' @param rankList List of ranked motifs
#' @param alpha The cardinality of the SAX alphabet
#' @return Heatmap dataset with seelected motifs
#' @import ggplot2
#' @import reshape2
#' @import scales
#' @import RColorBrewer
#' @importFrom grDevices grey.colors
#' @examples
#' #Launch all the workflow
#' #Plot the result
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,5)
#' stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' display_motifsDataset(dataset = STMotif::example_dataset, rstmotifs[c(1:4)],  5)
#' @export
#'
display_motifsDataset <- function(dataset,rankList,alpha){
  colorEncode <- 1:alpha
  datasetColor.Org <- as.matrix(dataset)
  datasetColor.Org <- as.vector(datasetColor.Org)
  datasetColor.Org <- STSNormalization(datasetColor.Org)
  mybin <- binning(datasetColor.Org, alpha)
  datasetColor.Org <- colorEncode[mybin$bins_factor]
  datasetColor.Org <- t(matrix(datasetColor.Org, nrow = nrow(dataset), ncol = ncol(dataset)))
  datasetColor.Org <- melt(datasetColor.Org)
  datasetColor.Org$motif <- FALSE

  palhetaCores <- brewer.pal(length(rankList), 'Spectral')

  motifs.plot <-data.frame("s"=NULL, "t"=NULL, "g"= NULL)
  for (pos in 1:length(rankList)){
    motifs.plot<- rbind(motifs.plot ,data.frame("s"=rankList[[pos]]$vecst$s, "t"=rankList[[pos]]$vecst$t, "g"= pos, "color"=palhetaCores[pos]))
  }

  datasetColor <- merge(datasetColor.Org, motifs.plot, by.x=c('Var1', 'Var2'), by.y=c('s', 't'), all.x = TRUE)
  datasetColor$motif[!is.na(datasetColor$g)] <- TRUE
  datasetColor$g <- NULL
  datasetColor$color <- as.character(datasetColor$color)

  ggplot(data=datasetColor, aes(x=datasetColor$Var1, y=datasetColor$Var2, fill=datasetColor$value, color=datasetColor$color))   + geom_raster() +
    scale_fill_gradientn(colours = c("white","dimgrey"), values = scales::rescale(1:alpha), limits=c(1,alpha)) +
    theme_bw() + ggtitle("") + xlab("Space") + ylab("Time") + scale_y_reverse() +
    guides(fill=FALSE, color=FALSE) +
    geom_point(colour = ifelse(datasetColor$motif,datasetColor$color,NA), size = 4, shape=15, show.legend = FALSE)
}
