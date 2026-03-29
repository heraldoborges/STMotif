# ============================================================================
# Core functionality tests
# ============================================================================

# --- Test data setup --------------------------------------------------------

D  <- STMotif::example_dataset
DS <- NormSAX(D, 5)


# --- NormSAX ----------------------------------------------------------------

test_that("NormSAX returns correct dimensions", {
  result <- NormSAX(D, 5)
  expect_equal(dim(result), dim(D))
})

test_that("NormSAX returns character (SAX) values", {
  result <- NormSAX(D, 5)
  expect_true(all(sapply(result, is.character) | sapply(result, is.factor)))
})

test_that("NormSAX uses correct alphabet size", {
  for (a in c(3, 5, 7)) {
    result <- NormSAX(D, a)
    unique_letters <- unique(unlist(result))
    expect_true(all(unique_letters %in% letters[1:a]),
                info = paste("alphabet size:", a))
  }
})

test_that("NormSAX preserves column names", {
  result <- NormSAX(D, 5)
  expect_equal(colnames(result), colnames(D))
})


# --- STSADatasetAdjust -----------------------------------------------------

test_that("STSADatasetAdjust returns correct dimensions", {
  result <- STSADatasetAdjust(D, 10, 4)
  expect_equal(nrow(result) %% 10, 0)
  expect_equal(ncol(result) %% 4, 0)
})

test_that("STSADatasetAdjust is identity when dimensions match", {
  result <- STSADatasetAdjust(D, 20, 12)
  expect_equal(dim(result), dim(D))
})


# --- SearchSTMotifs ---------------------------------------------------------

test_that("SearchSTMotifs returns a list", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  expect_type(stmotifs, "list")
})

test_that("SearchSTMotifs finds at least one motif", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  expect_true(length(stmotifs) > 0)
})

test_that("Each motif has required components", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  for (i in seq_along(stmotifs)) {
    expect_true("isaxcod" %in% names(stmotifs[[i]]),
                info = paste("motif", i, "missing isaxcod"))
    expect_true("recmatrix" %in% names(stmotifs[[i]]),
                info = paste("motif", i, "missing recmatrix"))
    expect_true("vecst" %in% names(stmotifs[[i]]),
                info = paste("motif", i, "missing vecst"))
  }
})

test_that("Motif vecst has s and t columns", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  for (i in seq_along(stmotifs)) {
    expect_true(all(c("s", "t") %in% colnames(stmotifs[[i]]$vecst)),
                info = paste("motif", i))
  }
})

test_that("Motif SAX code has correct length (w)", {
  w <- 4
  stmotifs <- SearchSTMotifs(D, DS, w, 5, 4, 10, 2, 2)
  for (i in seq_along(stmotifs)) {
    expect_equal(nchar(stmotifs[[i]]$isaxcod), w,
                 info = paste("motif", i))
  }
})


# --- RankSTMotifs -----------------------------------------------------------

test_that("RankSTMotifs returns a list with rank component", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  rstmotifs <- RankSTMotifs(stmotifs)
  expect_type(rstmotifs, "list")
  expect_true(length(rstmotifs) > 0)
  expect_true("rank" %in% names(rstmotifs[[1]]))
})

test_that("RankSTMotifs rank has required components", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  rstmotifs <- RankSTMotifs(stmotifs)
  rank_info <- rstmotifs[[1]]$rank
  expect_true(all(c("dist", "word", "qtd", "proj") %in% names(rank_info)))
})

test_that("RankSTMotifs preserves motif count", {
  stmotifs <- SearchSTMotifs(D, DS, 4, 5, 4, 10, 2, 2)
  rstmotifs <- RankSTMotifs(stmotifs)
  expect_equal(length(rstmotifs), length(stmotifs))
})

test_that("RankSTMotifs handles empty input", {
  rstmotifs <- RankSTMotifs(list())
  expect_equal(length(rstmotifs), 0)
})


# --- CSAMiningProcess -------------------------------------------------------

test_that("CSAMiningProcess returns ranked motifs", {
  rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  expect_type(rstmotifs, "list")
  expect_true(length(rstmotifs) > 0)
  expect_true("rank" %in% names(rstmotifs[[1]]))
})

test_that("CSAMiningProcess result matches step-by-step execution", {
  # Step-by-step
  ds_step <- NormSAX(D, 5)
  stm_step <- SearchSTMotifs(D, ds_step, 4, 5, 4, 10, 2, 2)
  rst_step <- RankSTMotifs(stm_step)
  
  # All-in-one
  rst_all <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  
  expect_equal(length(rst_all), length(rst_step))
  expect_equal(rst_all[[1]]$isaxcod, rst_step[[1]]$isaxcod)
})


# --- Visualization ----------------------------------------------------------

test_that("display_motifsSTSeries returns a ggplot", {
  rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  n <- min(4, length(rstmotifs))
  p <- display_motifsSTSeries(D, rstmotifs[seq_len(n)], space = c(1:4, 10:12))
  expect_s3_class(p, "ggplot")
})

test_that("display_motifsDataset returns a ggplot", {
  rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  n <- min(4, length(rstmotifs))
  p <- display_motifsDataset(D, rstmotifs[seq_len(n)], 5)
  expect_s3_class(p, "ggplot")
})

test_that("display_motifsSTSeries handles NULL rstmotifs", {
  p <- display_motifsSTSeries(D, NULL, space = 1:4)
  expect_s3_class(p, "ggplot")
})

test_that("display_motifsSTSeries handles empty rstmotifs", {
  p <- display_motifsSTSeries(D, list(), space = 1:4)
  expect_s3_class(p, "ggplot")
})

test_that("display_motifsSTSeries works with single series", {
  rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  p <- display_motifsSTSeries(D, rstmotifs[1], space = 1)
  expect_s3_class(p, "ggplot")
})

test_that("display_motifsDataset works with single motif", {
  rstmotifs <- CSAMiningProcess(D, DS, 4, 5, 4, 10, 2, 2)
  p <- display_motifsDataset(D, rstmotifs[1], 5)
  expect_s3_class(p, "ggplot")
})


# --- Internal helpers (accessed via :::) ------------------------------------

test_that(".wide_to_long produces correct output", {
  df <- data.frame(x = 1:3, A = c(10, 20, 30), B = c(40, 50, 60))
  result <- STMotif:::.wide_to_long(df)
  expect_equal(nrow(result), 6)
  expect_equal(colnames(result), c("x", "variable", "value"))
  expect_equal(result$value, c(10, 20, 30, 40, 50, 60))
  expect_equal(as.character(result$variable), c("A", "A", "A", "B", "B", "B"))
})

test_that(".melt_matrix produces correct output", {
  mat <- matrix(1:6, nrow = 2, ncol = 3)
  result <- STMotif:::.melt_matrix(mat)
  expect_equal(nrow(result), 6)
  expect_equal(colnames(result), c("Var1", "Var2", "value"))
  expect_equal(result$Var1, c(1, 2, 1, 2, 1, 2))
  expect_equal(result$Var2, c(1, 1, 2, 2, 3, 3))
  expect_equal(result$value, 1:6)
})

test_that(".safe_palette handles n < 3", {
  p1 <- STMotif:::.safe_palette(1)
  expect_length(p1, 1)
  expect_true(grepl("^#", p1))
  
  p2 <- STMotif:::.safe_palette(2)
  expect_length(p2, 2)
  
  p0 <- STMotif:::.safe_palette(0)
  expect_length(p0, 0)
})

test_that(".safe_palette returns n colors for n >= 3", {
  for (n in 3:8) {
    p <- STMotif:::.safe_palette(n)
    expect_equal(length(p), n, info = paste("n =", n))
  }
})
