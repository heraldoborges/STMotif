# STMotif 2.0.3
* Removed dependency on `reshape2` (retired from CRAN). 
* Added `rlang` to Imports; set minimum ggplot2 >= 3.4.0.
* Fixed `brewer.pal()` crash with fewer than 3 motifs.
* Fixed default `space` parameter in `display_motifsSTSeries()`.
* Replaced deprecated ggplot2 idioms (`guides(fill = FALSE)`, `aes()` with `$`).
* Replaced `T`/`F` with `TRUE`/`FALSE` throughout. 
* Renamed internal helpers to avoid masking `base::rank()` and false S3 method detection (`normalize.minmax`, `plot.series`).
* Replaced broad `@import` with selective `@importFrom`.

# STMotif 2.0.2

* Fix some functions.

# STMotif 2.0.1

* Paper about the technique implemented in the package published in Journal https://content.iospress.com/articles/intelligent-data-analysis/ida194759

# STMotif 2.0.0

* New example dataset.
* Fix some functions.
* New visualization functions.


# STMotif 1.0.1

* Fix some functions.
* Improvements of the documentation.

# STMotif 1.0.0

* Redefining the functions names.
* Improvements of the documentation.

# STMotif 0.1.1

* Improvements of the documentation.
* Fix the plot function

# STMotif 0.1.0

* Fist version of the STMotif package.
