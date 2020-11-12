## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ---- echo=FALSE---------------------------------------------------------
source(file = "../R/mainFunction.R")
source(file = "../R/subFunction.R")
source(file = "../R/visualization.R")
library(stats)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("STMotif")

## ---- eval = FALSE-------------------------------------------------------
#  library(STMotif)

## ---- echo=TRUE----------------------------------------------------------

# The process is launched on the provided example dataset
dim(D <- STMotif::example_dataset)

# Normalizartion and SAX indexing
DS <- NormSAX(D = STMotif::example_dataset,a =5)

# Information of the normalized and SAX indexing dataset 
# The candidates built 
head(NormSAX(D = STMotif::example_dataset, a = 5)[,1:10])


## ---- echo=TRUE----------------------------------------------------------
# The list of motifs 
# stmotifs <- SearchSTMotifs(D,DS,w,a,sb,tb,si,ka)
stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
stmotifs[[1]]

## ---- echo=TRUE----------------------------------------------------------
# The rank list of stmotifs 
rstmotifs <- RankSTMotifs(stmotifs)
rstmotifs[[1]]

## ---- echo=TRUE----------------------------------------------------------
# CSAMiningProcess
rstmotifs <- RankSTMotifs(stmotifs)
rstmotifs[[1]]

## ----fig, fig.height = 4, fig.width = 5, fig.align = "center"------------
display_motifsDataset(dataset = STMotif::example_dataset, rstmotifs[c(1:4)],  5)

## ----fig1, fig.height = 4, fig.width = 5, fig.align = "center"-----------
display_motifsSTSeries(dataset = STMotif::example_dataset,rstmotifs[c(1:4)],space = c(1:4,10:12))

