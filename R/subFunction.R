#' Binning of numeric values
#'
#' Discretizes a numeric vector into \code{a} bins using quantile-based
#' breakpoints.
#'
#' @param v Numeric vector to bin.
#' @param a Number of bins (alphabet size).
#' @return A list with components: \code{binning} (bin means), \code{bins_factor}
#'   (bin assignments), \code{q} (quantiles), \code{qf} (quantile intervals),
#'   \code{bins} (mapped values), \code{mse} (mean squared error).
#' @keywords internal
#' @noRd
binning <- function(v, a) {
  p <- seq(from = 0, to = 1, by = 1/a)
  q <- stats::quantile(v, p)
  qf <- matrix(c(q[1:(length(q)-1)],q[2:(length(q))]), ncol=2)
  vp <- cut(v, unique(q), FALSE, include.lowest=TRUE)
  m <- tapply(v, vp, mean)
  vm <- m[vp]
  mse <- mean( (v - vm)^2, na.rm = TRUE)
  return (list(binning=m, bins_factor=vp, q=q, qf=qf, bins=vm, mse=mse))
}


#' Z-score normalization
#'
#' @param vector Numeric vector.
#' @return Normalized vector with mean 0 and standard deviation 1.
#' @keywords internal
#' @noRd
STSNormalization <- function(vector) {
  return ((vector - mean(vector, na.rm = TRUE)) / stats::sd(vector, na.rm = TRUE))
}


#' SAX encoding of a dataset
#'
#' Encodes numeric values into SAX letters using quantile-based binning.
#'
#' @param dataset Original dataset (used for dimensions).
#' @param vector Numeric vector (typically flattened and normalized dataset).
#' @param a Alphabet size (number of SAX letters).
#' @return A data frame with the same dimensions as \code{dataset}, containing
#'   SAX letter encodings.
#' @keywords internal
#' @noRd
STSSaxEncode <- function(dataset, vector, a) {
  mybin <- binning(vector, a)
  myletters <- letters[1:a]
  saxvector <- myletters[mybin$bins_factor]
  saxvector = matrix(saxvector, nrow = nrow(dataset), ncol = ncol(dataset))
  saxvector = data.frame(saxvector)
  colnames(saxvector) = colnames(dataset)
  return(saxvector)
}


#' Compute spatio-temporal blocks
#'
#' Divides the dataset into rectangular blocks of size \code{tb} (temporal)
#' by \code{sb} (spatial).
#'
#' @param dataset Dataset (data frame or matrix).
#' @param tb Temporal block size (number of rows per block).
#' @param sb Spatial block size (number of columns per block).
#' @return A list with components: \code{datasets} (list of block data frames),
#'   \code{nrows} (number of block rows), \code{ncols} (number of block columns),
#'   \code{rectangles} (list of block coordinates).
#' @keywords internal
#' @noRd
STSComputeBlocks <- function(dataset, tb, sb) {
  datasets <- list()
  rectangles <- list()
  
  c = ncol(dataset)
  r = nrow(dataset)
  nc = c / sb
  nr = r / tb
  i = 1
  j = 1
  n = 1
  for (i in 1:nc) {
    sc = (i-1)*sb + 1
    ec = (i)*sb
    for (j in 1:nr) {
      sr = (j-1)*tb + 1
      er = (j)*tb
      ds = dataset[sr:er, sc:ec]
      datasets[[n]] = ds
      rect = c(sS = sc, eS = ec, sT = sr, eT = er, nr = j, nc = i)
      rectangles[[n]] = rect
      n = n + 1
    }
  }
  blocks = list(datasets = datasets, nrows = nr, ncols = nc, rectangles = rectangles)
  return (blocks)
}


#' Identify motifs within a single block
#'
#' Takes a flattened block and its SAX encoding, discovers frequent
#' subsequences of length \code{w}.
#'
#' @param ts Numeric vector (flattened block).
#' @param tss Character vector (flattened SAX block).
#' @param w Word size (motif length).
#' @param tb Temporal block size.
#' @param a Alphabet size.
#' @return A list with components: \code{Subs.SAX} (all subsequences),
#'   \code{Motif.SAX} (frequent subsequences), \code{Indices} (positions).
#' @keywords internal
#' @noRd
identifyMotifsInBlock <- function(ts, tss, w, tb, a) {
  # Generate all possible subsequences
  ts.sax <- NULL
  for (i in 1:length(tss)){
    if(floor((i-1)/tb)==floor((i-1+w-1)/tb)){ # Check boundary
      ts.sax  <- rbind(ts.sax ,c(i,tss[i:(i+w-1)]) )
    }
  }
  
  ts.sax <- stats::na.omit(ts.sax)
  ts.sax <- as.data.frame(ts.sax, stringsAsFactors = FALSE)
  
  colnames(ts.sax) <- c("StartPosition", 1:w)
  ts.sax$StartPosition <- as.numeric(ts.sax$StartPosition)
  
  # Group identical SAX words and collect their start positions
  i = j <- 1
  indices <- list()
  for (i in 1:nrow(ts.sax)){
    saxMotif <- paste(ts.sax[i,-1], collapse = "")
    indices[[saxMotif]] <- c(indices[[saxMotif]],ts.sax[i,1])
  }
  while (j <= length(indices)){ # Remove singletons
    if(length(indices[[j]])<=1){indices[[j]]<-NULL}else{j<-j+1}
  }
  
  # Build sub-matrices for each motif
  motif.sax <- NULL
  if (length(indices)>0){
    for (i in 1:length(indices)){
      motif.sax[[i]] <- ts.sax[which(ts.sax[,1] %in% indices[[i]]),]
    }
  }
  
  return(list(Subs.SAX=ts.sax, Motif.SAX=motif.sax, Indices=indices))
}


#' Handle motifs from one block
#'
#' Validates Block Occurrence (BO) and Block Spatial Occurrence (BSO)
#' constraints and accumulates motifs across blocks.
#'
#' @param stmotifs Accumulated list of spatial-time motifs.
#' @param motif Motifs found in the current block.
#' @param nrows Number of block rows.
#' @param ncols Number of block columns.
#' @param rectangle Coordinates of the current block.
#' @param ka Minimum number of spatial series with occurrences (kappa).
#' @param si Minimum number of total occurrences (sigma).
#' @return Updated list of spatial-time motifs.
#' @keywords internal
#' @noRd
STSIdentifySTMotif <- function(stmotifs, motif, nrows, ncols, rectangle, ka, si) {
  k <- length(stmotifs)
  
  # Block coordinates
  sS = rectangle["sS"]
  eS = rectangle["eS"]
  sT = rectangle["sT"]
  eT = rectangle["eT"]
  nr = rectangle["nr"]
  nc = rectangle["nc"]
  
  recMatrix = matrix(rep(0, nrows*ncols), nrow = nrows, ncol = ncols)
  tb <- eT - sT + 1
  sb <- eS - sS + 1
  
  if(length(motif$Indices)>0){
    for(a in 1:length(motif$Indices)){
      vec <- motif$Indices[[a]]
      
      # BO - Block Occurrences validation
      if(length(vec) >= si) {
        scount <- rep(0, sb)
        
        for(z in 1: length(vec)) {
          i <- as.integer(vec[z] / tb) + 1
          scount[i] <- 1
        }
        
        # BSO - Block Spatial Occurrences Validation
        if(sum(scount) >= ka) {
          isaxcod <- paste(motif$Motif.SAX[[a]][1,2:(length(motif$Subs.SAX))], collapse = "")
          
          vect <- as.integer(vec %% tb) + sT - 1
          vecs <- as.integer(vec / tb) + sS
          
          i <- match(isaxcod, names(stmotifs))
          if (is.na(i)) {
            k = k + 1
            stmotifs[[k]] <- list(isaxcod=isaxcod, vecs=vecs, vect=vect, recmatrix=recMatrix)
            stmotifs[[k]]$recmatrix[nr, nc] = 1
            names(stmotifs)[k] = isaxcod
          }
          else {
            list <- stmotifs[[i]]
            list$recmatrix[nr, nc] = max(list$recmatrix)+1
            list$vect <- c(list$vect, vect)
            list$vecs <- c(list$vecs, vecs)
            stmotifs[[i]] <- list
          }
        }
      }
    }
  }
  return (stmotifs)
}


#' Identify tight (connected) spatial-time motifs
#'
#' Removes isolated motifs by checking adjacency of blocks in the
#' recurrence matrix. Splits disconnected components into separate motifs.
#'
#' @param stmotif A single spatial-time motif with its recurrence matrix.
#' @param rectangles List of block coordinates.
#' @return A list of tight motifs (connected components).
#' @keywords internal
#' @noRd
STSIdentifyTightSTMotif <- function(stmotif, rectangles) {
  tight <- list()
  mat <- stmotif$recmatrix
  vecst <- stmotif$vecst
  
  for (i in 1:nrow(mat)) {
    for (j in 1:(ncol(mat)-1)) {
      if (mat[i,j] != 0) {
        iP <- i + 1
        jP <- j + 1
        if ((iP <= nrow(mat)) && (mat[iP,j] != 0)) {
          k <- min(mat[iP,j], mat[i,j])
          mat[mat == mat[iP,j] | mat == mat[i,j]] = k
        }
        if ((jP <= ncol(mat)) && (mat[i,jP] != 0)) {
          k <- min(mat[i,jP], mat[i,j])
          mat[mat == mat[i,jP] | mat == mat[i,j]] = k
        }
        if ((iP <= nrow(mat)) && (mat[iP,j] != 0) && (jP <= ncol(mat)) && (mat[i,jP] != 0)) {
          k <- min(mat[iP,jP], mat[i,j])
          mat[mat == mat[iP,jP] | mat == mat[i,j]] = k
        }
      }
    }
  }
  
  vec <- as.vector(mat)
  vec <- vec[vec > 0]
  vec <- unique(vec)
  k <- 1
  stmotif_org <- stmotif
  
  for (i in (vec)) {
    stmotif <- stmotif_org
    stmotif$recmatrix[mat != i] <- 0
    stmotif$recmatrix[mat == i] <- k
    vecrects <- as.vector(stmotif$recmatrix)
    rects <- rectangles[vecrects>0]
    stmotif$vecst <- vecst
    conds = rep(FALSE, nrow(stmotif$vecst))
    for (rect in (rects)) {
      sS = rect["sS"]
      eS = rect["eS"]
      sT = rect["sT"]
      eT = rect["eT"]
      conds = conds | (stmotif$vecst$s >= sS & stmotif$vecst$s <= eS & stmotif$vecst$t >= sT & stmotif$vecst$t <= eT)
    }
    stmotif$vecst <- stmotif$vecst[conds,]
    tight[[k]] <- stmotif
    k <- k + 1
  }
  return(tight)
}


#' Euclidean distance-based ranking component
#'
#' Computes a proximity score based on a greedy minimum spanning tree
#' of the motif occurrence positions.
#'
#' @param data Data frame with columns \code{space} and \code{time}.
#' @return Inverse of the mean edge weight (higher = closer occurrences).
#' @keywords internal
#' @noRd
comp_distance <- function(data) {
  nv <- nrow(data)
  na <- nrow(data)*(nrow(data)-1)/2
  
  ver <- rep(FALSE, nv)
  adj_mat <- matrix(0, nrow = na, ncol=3)
  k <- 0
  for (i in (1:(nv-1))) {
    for (j in ((i+1):nv)) {
      k <- k + 1
      adj_mat[k, 1] <- i
      adj_mat[k, 2] <- j
      adj_mat[k, 3] <- sqrt((data$space[i]-data$space[j])^2+(data$time[i]-data$time[j])^2)
    }
  }
  adj_mat <- data.frame(s = adj_mat[,1],d = adj_mat[,2], w = adj_mat[,3])
  o <- order(adj_mat$w)
  adj_mat <- adj_mat[o,]
  
  edges <- NULL
  for (k in 1:na) {
    i <- adj_mat$s[k]
    j <- adj_mat$d[k]
    if (!ver[i] | !ver[j]) {
      ver[i] <- TRUE
      ver[j] <- TRUE
      edges <- rbind(edges, adj_mat[k,])
    }
  }
  return(1/mean(edges$w))
}


#' SAX word entropy ranking component
#'
#' Computes the Shannon entropy of the character distribution in a
#' SAX-encoded motif word.
#'
#' @param str Character string (SAX word).
#' @return Shannon entropy in bits.
#' @keywords internal
#' @noRd
comp_word <- function(str) {
  x <- strsplit(str, "^")
  x <- x[[1]]
  n <- length(x)
  x <- table(x)
  x <- x / n
  y <- 0
  for (i in 1:length(x)) {
    y <- y - x[i]*log(x[i],2)
  }
  return(y)
}


#' Min-max normalization
#'
#' Scales numeric columns to the [0, 1] range.
#'
#' @param data Data frame with numeric columns.
#' @param norm.set Optional pre-computed normalization parameters.
#' @return A list with \code{data} (normalized) and \code{norm.set} (parameters).
#' @keywords internal
#' @noRd
normalize_minmax <- function(data, norm.set = NULL)
{
  data = data.frame(data)
  nums = unlist(lapply(data, is.numeric))
  data = data[ , nums]
  
  if(is.null(norm.set))
  {
    minmax = data.frame(t(sapply(data, max, na.rm = TRUE)))
    minmax = rbind(minmax, t(sapply(data, min, na.rm = TRUE)))
    colnames(minmax) = colnames(data)
    rownames(minmax) = c("max", "min")
  }
  else {
    minmax = norm.set
  }
  for (i in 1:ncol(data))
    data[,i] = (data[,i] - minmax["min", i]) / (minmax["max", i] - minmax["min", i])
  return (list(data=data, norm.set=minmax))
}


#' Rank motifs by multidimensional projection
#'
#' Projects the ranking metrics (distance, word entropy, quantity) onto
#' a single score and sorts the motifs in decreasing order.
#'
#' @param dataRank Data frame with columns \code{dist}, \code{word}, \code{qtd}.
#' @param stmotifs List of motifs to rank.
#' @return Ordered list of motifs with \code{rank} component added.
#' @keywords internal
#' @noRd
rank_motifs <- function(dataRank, stmotifs)
{
  dataRankOrg <- dataRank
  
  dataRank <- normalize_minmax(dataRank)$data
  
  dataRank = as.matrix(dataRank)
  
  transf <- rep(sqrt(0.5), ncol(dataRank))
  
  dataRankOrg$proj = dataRank %*% transf
  
  # Order by decreasing projection score
  o <- order(dataRankOrg$proj, decreasing = TRUE)
  stmotifsRank <- list()
  for (i in 1:length(stmotifs)) {
    indice <- o[i]
    stmotifs[[indice]][["rank"]] <- c(dataRankOrg[indice,]['dist'], dataRankOrg[indice,]['word'], dataRankOrg[indice,]['qtd'], dataRankOrg[indice,]['proj'])
    stmotifsRank[[i]] <- stmotifs[[indice]]
  }
  return (stmotifsRank)
}
