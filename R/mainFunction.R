#' Adjust a Dataset
#'
#' Adjusts the dimensions of a dataset so that it can be evenly divided
#' into spatio-temporal blocks of size \code{tb} x \code{sb}.
#'
#' @param D Dataset containing numeric values.
#' @param tb Temporal block size (number of rows per block).
#' @param sb Spatial block size (number of columns per block).
#' @return Dataset with rows and columns trimmed to be divisible by
#'   \code{tb} and \code{sb}, respectively.
#' @examples
#' D <- STSADatasetAdjust(STMotif::example_dataset, 20, 12)
#' @export
STSADatasetAdjust  <- function(D, tb, sb) {
  c = ncol(D)
  r = nrow(D)
  ec = c %% sb
  er = r %% tb
  D = D[1:(r-er), 1:(c-ec)]
  return (D)
}


#' CSA Datamining Process
#'
#' Performs the complete Combined Series Approach (CSA) workflow:
#' normalization with SAX encoding, motif discovery, and ranking.
#' This is a convenience wrapper around \code{\link{NormSAX}},
#' \code{\link{SearchSTMotifs}}, and \code{\link{RankSTMotifs}}.
#'
#' @param D Dataset containing numeric values.
#' @param DS Dataset containing SAX encoded values (recomputed internally;
#'   this parameter is kept for backward compatibility).
#' @param w Word size (motif length in SAX symbols).
#' @param a Number of letters in the SAX alphabet.
#' @param sb Spatial block size (number of columns per block).
#' @param tb Temporal block size (number of rows per block).
#' @param si Minimum number of occurrences inside each block (sigma).
#' @param ka Minimum number of spatial series with occurrences inside
#'   each block (kappa).
#' @return A list of ranked motifs. Each motif contains:
#' \describe{
#'   \item{isaxcode}{Motif sequence in character format.}
#'   \item{recmatrix}{Matrix indicating which blocks contain this motif.}
#'   \item{vecst}{Data frame with columns \code{s} (spatial) and \code{t}
#'     (temporal) giving the start positions of the motif in the original
#'     dataset.}
#'   \item{rank}{List with ranking components: \code{dist}, \code{word},
#'     \code{qtd}, \code{proj}.}
#' }
#' @examples
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' rmotif <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
#' @export
CSAMiningProcess <- function(D, DS, w, a, sb, tb, si, ka) {
  DS <- NormSAX(D, a)
  stmotifs <- SearchSTMotifs(D, DS, w, a, sb, tb, si, ka)
  rstmotifs <- RankSTMotifs(stmotifs)
  return(rstmotifs)
}


#' Normalize and SAX Encode a Dataset
#'
#' Applies z-score normalization to the entire dataset and encodes
#' the values using the Symbolic Aggregate approXimation (SAX) with
#' an alphabet of size \code{a}.
#'
#' @param D Dataset containing numeric values.
#' @param a Number of letters in the SAX alphabet.
#' @return A data frame with the same dimensions as \code{D}, containing
#'   SAX letter encodings (characters from \code{a} to the \code{a}-th
#'   letter of the alphabet).
#' @examples
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' @export
NormSAX <- function(D, a) {
  vector <- as.matrix(D)
  vector <- as.vector(vector)
  vectorNorm <- (vector - mean(vector, na.rm = TRUE)) / stats::sd(vector, na.rm = TRUE)
  DS <- STSSaxEncode(D, vectorNorm, a)
  return(DS)
}


#' Search for Spatial-Time Motifs
#'
#' Discovers motifs in the spatio-temporal blocks of the dataset,
#' validates occurrence constraints, and groups motifs from neighboring
#' blocks.
#'
#' @param D Dataset containing numeric values.
#' @param DS Dataset containing SAX encoded values (as returned by
#'   \code{\link{NormSAX}}).
#' @param w Word size (motif length in SAX symbols).
#' @param a Number of letters in the SAX alphabet.
#' @param sb Spatial block size (number of columns per block).
#' @param tb Temporal block size (number of rows per block).
#' @param si Minimum number of occurrences inside each block (sigma).
#'   Default: 3.
#' @param ka Minimum number of spatial series with occurrences inside
#'   each block (kappa). Default: 3.
#' @return A list of identified motifs. Each motif contains:
#' \describe{
#'   \item{isaxcode}{Motif sequence in character format.}
#'   \item{recmatrix}{Matrix indicating which blocks contain this motif.}
#'   \item{vecst}{Data frame with columns \code{s} and \code{t} giving
#'     the start positions in the original dataset.}
#' }
#' @examples
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
#' @export
SearchSTMotifs <- function(D, DS, w, a, sb, tb, si = 3, ka = 3) {
  
  saxblocks <- STSComputeBlocks(DS, tb, sb)
  saxblocks$rectangles <- NULL
  
  blocks <- STSComputeBlocks(D, tb, sb)
  nrows = blocks$nrows
  ncols = blocks$ncols
  rectangles = blocks$rectangles
  blocks$rectangles <- NULL
  
  motifs <- list()
  size = length(blocks$datasets)
  for (i in 1:size) {
    block = blocks$datasets[[i]]
    saxblock = saxblocks$datasets[[i]]
    block = as.vector(as.matrix(block))
    saxblock = as.vector(as.matrix(saxblock))
    motifs[[i]] <- identifyMotifsInBlock(ts = block, tss = saxblock, tb = tb, w = w, a = a)
  }
  
  stmotifs <- list()
  for (i in 1:length(motifs)) {
    stmotifs <- STSIdentifySTMotif(stmotifs, motifs[[i]], nrows, ncols, rectangles[[i]], ka = ka, si = si)
  }
  
  sttightmotifs <- list()
  
  if (length(stmotifs)>0){
    for (i in 1:length(stmotifs)) {
      stmotif = stmotifs[[i]]
      s = stmotif$vecs
      t = stmotif$vect
      stmotif$vecst = data.frame(s, t)
      stmotif$vecs <- NULL
      stmotif$vect <- NULL
      stmotifs[[i]] = stmotif
    }
    
    for(stmotif in (stmotifs)) {
      sttightmotifsSplit <- STSIdentifyTightSTMotif(stmotif, rectangles)
      for (item in (sttightmotifsSplit)) {
        pos = length(sttightmotifs)+1
        sttightmotifs[[pos]] <- item
        names(sttightmotifs)[pos] = item$isaxcod
      }
    }
  }
  return (sttightmotifs)
}


#' Rank Spatial-Time Motifs
#'
#' Ranks the discovered motifs by computing a composite quality score
#' that balances spatial-temporal proximity of occurrences, entropy of
#' the SAX encoding, and quantity of occurrences.
#'
#' @param stmotifs List of identified motifs (as returned by
#'   \code{\link{SearchSTMotifs}}).
#' @return A list of motifs sorted by decreasing quality score. Each motif
#'   gains a \code{rank} component with \code{dist}, \code{word},
#'   \code{qtd}, and \code{proj} values.
#' @examples
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' @export
RankSTMotifs <- function(stmotifs) {
  rstmotifs <- list()
  if(length(stmotifs)>0){
    dataRank <- NULL
    for (i in 1:length(stmotifs)) {
      s <- stmotifs[[i]][["vecst"]][["s"]]
      t <- stmotifs[[i]][["vecst"]][["t"]]
      word <- stmotifs[[i]]$isaxcod
      occurrences <- data.frame(space = s, time = t)
      distance_rank <- comp_distance(occurrences)
      word_rank <- comp_word(stmotifs[[i]]$isaxcod)
      qtd_rank <- log(nrow(occurrences), base = 2)
      dataRank <- rbind(dataRank, data.frame(dist = distance_rank, word = word_rank, qtd = qtd_rank))
    }
    rownames(dataRank) <- c(1:length(stmotifs))
    rstmotifs <- rank_motifs(dataRank, stmotifs)
  }
  
  return(rstmotifs)
}
