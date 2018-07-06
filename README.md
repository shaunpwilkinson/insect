# insect

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/insect)](https://cran.r-project.org/package=insect)
[![CRAN_Downloads_Badge](http://cranlogs.r-pkg.org/badges/grand-total/insect)](https://cran.r-project.org/package=insect)
[![DOI](https://zenodo.org/badge/87808693.svg)](https://zenodo.org/badge/latestdoi/87808693)
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0002--7332--7931-brightgreen.svg)](http://orcid.org/0000-0002-7332-7931) 
[![Build_Status](https://travis-ci.org/shaunpwilkinson/insect.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/insect)
[![Project_Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

### Informatic sequence classification trees

`insect` is an R package for DNA meta-barcoding analysis. It provides a 
bioinformatics pipeline that automates the process from HTS sequence 
de-multiplexing and quality filtering to probabilistic taxonomic 
assignment using informatic sequence classification trees. 
It also contains functions for searching and downloading sequences 
from GenBank, a "virtual PCR" tool, a new algorithm for classification tree learning, 
and a method for classifying DNA barcode sequences using pre-computed trees.


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


### Help

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


