# insect

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/insect)](https://cran.r-project.org/package=insect)
[![CRAN_Downloads_Badge](http://cranlogs.r-pkg.org/badges/grand-total/insect)](https://cran.r-project.org/package=insect)
[![DOI](https://zenodo.org/badge/87808693.svg)](https://zenodo.org/badge/latestdoi/87808693)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--7332--7931-brightgreen.svg)](http://orcid.org/0000-0002-7332-7931) 
[![Build_Status](https://travis-ci.org/shaunpwilkinson/insect.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/insect)
[![Project_Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)


### Informatic sequence classification trees

`insect` is an R package for taxonomic identification of amplicon 
sequence variants generated during DNA meta-barcoding analysis. 
The learning and classification algorithms implemented in the 
package are based on full probabilistic models (profile hidden Markov models) 
and offer highly accurate taxon IDs, albeit at a relatively high computational cost.

The package also contains functions for searching and downloading reference 
sequences and taxonomic information from NCBI, 
a "virtual PCR" tool for sequence trimming, 
a function for "purging" erroneously labeled reference sequences, 
and several other handy tools.  

`insect` is designed to be used in conjunction with the 
[dada2](https://benjjneb.github.io/dada2/index.html) pipeline or any other
de-noising tool that produces a list of amplicon sequence variants (ASVs). 
While unfiltered sequences can also be processed with high accuracy, 
the **insect** classification algorithm is relatively slow, 
since it uses a computationally intensive dynamic
programming algorithm to find the likelihood values
of each sequence given the models at each node of the classification tree. 
Hence an appropriately filtered input dataset will generally be 
much faster to process.



### Installation

To download **insect** from CRAN and load the package, run

```R
install.packages("insect")
library(insect)
```

To download the latest development version from GitHub, run:

```R
devtools::install_github("shaunpwilkinson/insect", build_vignettes = TRUE) 
library("insect")
```


### Classifying sequences

An overview of the package and its functions can be found by running

```R
vignette("insect-vignette")
```

If you experience a problem using this package please feel free to
raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/insect/issues).


### Acknowledgements

This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.


