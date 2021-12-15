# LOMAR

An R package to deal with point sets from single molecule localization microscopy.
This package provides three sets of functionalities:
  - data input: read SMLM data as point sets either from csv files or from TIFF images
  - registration: point sets registration using different algorithms
  - topological data analysis: compute similarity between point sets using persistent homology

# Installation

To install this package, run (from within R):

``` R
library(devtools)
install_git('https://git.embl.de/heriche/lomar')
```

This package depends on these other packages:
  * data.table
  * TDA
  * foreach
  * parallel
  * doParallel
  * proxy
  * reshape2
  * pracma
  * transport
  * RANN
  * ff
  * aws
  * dbscan
  * EBImage (from Bioconductor)
  
