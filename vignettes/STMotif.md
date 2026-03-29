---
title: "STMotif R Package"
author: "Heraldo Borges, Amin Bazaz, Eduardo Ogasawara"
date: "2020-11-12"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Spatial-Time Motif Discovery with STMotif}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  
---





The goal of the `STSMotif` R package is to allows the discovery and ranking of a motif in spatial-time series quickly and efficiently.





### Introduction

A pattern that significantly occurs in a time series is called a motif. In spatial time series data, these patterns may not be substantially present in a single time series but dispersed over several times series, limited in both space and time. The `STMotif R package` was developed to simplify the Spatio-temporal data mining on the search for these motifs. We present the functions available in `STMotif package` through the sample dataset, also available in this package.

First, install the package by typing:


```r
install.packages("STMotif")
```

Then, load the package by typing:


```r
library(STMotif)
```

It provides two categories of functions: for discovering and ranking motifs (CSAMiningProcess) and functions for viewing the identified motifs.


### 1. CSAMiningProcess

 1. The function `NormSAX` allows the normalization and SAX indexing of the dataset.
 

```r

# The process is launched on the provided example dataset
dim(D <- STMotif::example_dataset)
#> [1] 20 12

# Normalizartion and SAX indexing
DS <- NormSAX(D = STMotif::example_dataset,a =5)

# Information of the normalized and SAX indexing dataset 
# The candidates built 
head(NormSAX(D = STMotif::example_dataset, a = 5)[,1:10])
#>                      
#> 1 a c c c c c c c e c
#> 2 a a e c e e e c c e
#> 3 c e e e c e d e e e
#> 4 e e b e e d e e d b
#> 5 e c c b b c b c a e
#> 6 b d c a a a b e a d
```



2. The function `SearchSTMotifs` allows to check and filter the stmotifs, grouping the motifs from the neighboring block. 



```r
# The list of motifs 
# stmotifs <- SearchSTMotifs(D,DS,w,a,sb,tb,si,ka)
stmotifs <- SearchSTMotifs(D,DS,4,5,4,10,2,2)
stmotifs[[1]]
#> $isaxcod
#> [1] "ceeb"
#> 
#> $recmatrix
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    0    0
#> 
#> $vecst
#>   s t
#> 1 1 3
#> 2 3 1
#> 3 4 2
```



3. The function `RankSTMotifs` allows to rank the stmotifs list, making a balance between distance among the occurrences of a motif with the encoded information on the motif itself and his quantity. 


```r
# The rank list of stmotifs 
rstmotifs <- RankSTMotifs(stmotifs)
rstmotifs[[1]]
#> $isaxcod
#> [1] "bded"
#> 
#> $recmatrix
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    1    1    1
#> 
#> $vecst
#>    s  t
#> 1  1 11
#> 2  2 11
#> 3  4 17
#> 4  5 17
#> 5  8 15
#> 6 10 15
#> 7 12 12
#> 
#> $rank
#> $rank$dist
#> [1] 0.5259316
#> 
#> $rank$word
#> [1] 1.5
#> 
#> $rank$qtd
#> [1] 2.807355
#> 
#> $rank$proj
#>       [,1]
#> 3 1.522208
```


4.All this process can be summarized in the function `CSAMiningProcess` which performs all the steps listed above.


```r
# CSAMiningProcess
stmotifs <- CSAMiningProcess(D,DS,4,5,4,10,2,2)
rstmotifs[[1]]
#> $isaxcod
#> [1] "bded"
#> 
#> $recmatrix
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    1    1    1
#> 
#> $vecst
#>    s  t
#> 1  1 11
#> 2  2 11
#> 3  4 17
#> 4  5 17
#> 5  8 15
#> 6 10 15
#> 7 12 12
#> 
#> $rank
#> $rank$dist
#> [1] 0.5259316
#> 
#> $rank$word
#> [1] 1.5
#> 
#> $rank$qtd
#> [1] 2.807355
#> 
#> $rank$proj
#>       [,1]
#> 3 1.522208
```


### 2. Visualization

- Plot a heatmap of the dataset and highlight the selected motifs from the list


```r
display_motifsDataset(dataset = STMotif::example_dataset, rstmotifs[c(1:4)],  5)
```

<img src="figure/fig-1.png" title="plot of chunk fig" alt="plot of chunk fig" style="display: block; margin: auto;" />



- Plot the selected spatial-time series with the selected motifs highlighted


```r
display_motifsSTSeries(dataset = STMotif::example_dataset,rstmotifs[c(1:4)],space = c(1:4,10:12))
```

<img src="figure/fig1-1.png" title="plot of chunk fig1" alt="plot of chunk fig1" style="display: block; margin: auto;" />


