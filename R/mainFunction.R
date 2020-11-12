#' Adjust a Dataset
#' Adjust the dimensions of a dataset to build the blocks
#' @param D Dataset containing numeric values
#' @param tb Temporal block size
#' @param sb Spatial block size
#' @return Dataset adjusted to build the blocks.
#' @examples
#' #Adjust a block
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


#'  CSAMiningProcess
#'
#' CSA Datamining Process
#' @param D Dataset containing numeric values
#' @param DS Dataset containing SAX encoded values
#' @param w Word Size
#' @param a Number of letters to do the encode
#' @param sb Spatial block size
#' @param tb Temporal block size
#' @param si Minimum number of occurrences inside each block
#' @param ka Minimum number of spatial-time series with occurrences inside each block
#' @return Return a list of ranked motifs. Each motif contains the information [isaxcode, recmatrix, vectst, rank], as described:
#' @return isaxcode: Motif sequences in character format
#' @return recmatrix: Matrix giving as information the blocks containing this motif
#' @return vectst: Coordinate of the start positions of the motif in the original dataset
#' @return rank: L of information used for motif ranking, as [dist, word, qtd, proj]
#' @examples
#' #CSA Datamining process
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,5)
#' rmotif <- CSAMiningProcess(D,DS,4,5,4,10,2,2)
#' @export
CSAMiningProcess <- function (D,DS,w,a,sb,tb,si,ka){
  DS <- NormSAX(D,a)
  stmotifs <- SearchSTMotifs(D,DS,w,a,sb,tb,si,ka)
  rstmotifs <- RankSTMotifs(stmotifs)
  return(rstmotifs)
}


#' Normalize the data and SAX indexing
#' @param D Dataset containing numeric values
#' @param a Number of letters use to encode
#' @return A normalized and encoded dataset for a given alphabet a
#' @examples
#' #Normalization and Sax Dataset
#' DS <- NormSAX(STMotif::example_dataset, 5)
#' @export
NormSAX <- function (D,a){
  vector <- as.matrix(D)
  vector <- as.vector(vector)
  vectorNorm <- (vector-mean(vector, na.rm = T))/stats::sd(vector, na.rm = T)
  DS <- STSSaxEncode(D, vectorNorm, a)
  return (DS)
}




#'  SearchSTMotifs
#'
#' Search for Spatial-time Motifs
#' @param D Dataset containing numeric values
#' @param DS Dataset containing SAX encoded values
#' @param w Word Size
#' @param a Number of letters to do the encode
#' @param sb "Space slice" Number of columns in each block
#' @param tb "Time slice" Number of rows in each block
#' @param si Support of Global Occurrence (GO)
#' @param ka Support for Spatial Occurrence (SO)
#' @return Return a list of identified motifs. Each motif contains the information [isaxcode, recmatrix, vectst], as described:
#' @return isaxcode: Motif sequences in character format
#' @return recmatrix: Matrix giving as information the blocks containing this motif
#' @return vectst: Coordinate of the start positions of the motif in the original dataset
#' @examples
#' #Search for Spatial-time Motifs
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,5)
#' stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
#' @export
SearchSTMotifs <- function (D,DS,w,a,sb,tb,si=3,ka=3){

  saxblocks <- STSComputeBlocks(DS, tb, sb)
  saxblocks$rectangles <- NULL

  blocks <- STSComputeBlocks(D, tb, sb)
  nrows = blocks$nrows
  ncols = blocks$ncols
  rectangles = blocks$rectangles
  blocks$rectangles <- NULL

  motifs<-list()
  size=length(blocks$datasets)
  for (i in 1:size) {
    block = blocks$datasets[[i]]
    saxblock = saxblocks$datasets[[i]]
    block = as.vector(as.matrix(block))
    saxblock = as.vector(as.matrix(saxblock))
    motifs[[i]] <- identifyMotifsInBlock(ts = block, tss = saxblock, tb = tb ,w = w, a = a)
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




#' Rank the STmotifs
#' Rank motifs by their quality
#' @param stmotifs List of identified motifs
#' @return The ranked version of the identified list of motifs
#' @examples
#' #Search for Spatial-time Motifs
#' D  <- STMotif::example_dataset
#' DS <- NormSAX(STMotif::example_dataset,5)
#' stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
#' rstmotifs <- RankSTMotifs(stmotifs)
#' @export
RankSTMotifs <- function(stmotifs) {
  rstmotifs<-list()
  if(length(stmotifs)>0){
    dataRank <- NULL
    for (i in 1:length(stmotifs)) {
      s <- stmotifs[[i]][["vecst"]][["s"]]
      t <- stmotifs[[i]][["vecst"]][["t"]]
      word <- stmotifs[[i]]$isaxcod
      occurrences<- data.frame(space = s, time = t)
      distance_rank <- comp_distance(occurrences)
      word_rank <- comp_word(stmotifs[[i]]$isaxcod)
      qtd_rank <- log(nrow(occurrences), base=2)
      dataRank <- rbind(dataRank, data.frame(dist = distance_rank, word = word_rank, qtd=qtd_rank))
    }
    rownames(dataRank) <- c(1:length(stmotifs))
    rstmotifs <- rank(dataRank,stmotifs)
  }

  return(rstmotifs)
}




