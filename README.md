# insect

[![Build Status](https://travis-ci.org/shaunpwilkinson/insect.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/insect)

### Informatic sequence classification trees

`insect` is an R package for DNA meta-barcoding analysis. It provides a 
bioinformatics pipeline that automates the process from HTS sequence 
de-multiplexing and quality filtering to probabilistic taxonomic 
assignment using informatic sequence classification trees. 
It also contains functions for searching and downloading sequences 
from GenBank, a "virtual PCR" tool, a new algorithm for tree learning, 
and a method for classifying DNA barcode sequences using pre-computed trees.

### Installation
`insect` is currently only available as a development version, with a stable
release available on CRAN shortly. To download the package from 
GitHub you will first need to ensure you have a C/C++ compliler and the 
[devtools](https://github.com/hadley/devtools) R package installed. 
Linux users will generally have a compiler such as `gcc` installed by default; 
however Windows users will need to download 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac 
OSX users will need [Xcode](https://developer.apple.com/xcode) 
(note that Rtools and Xcode are not R packages). To download and install 
devtools, run 
```R
install.packages("devtools")
``` 
and then install and load `insect` by running 
```R
devtools::install_github("shaunpwilkinson/insect") 
library("insect")
```

### Help
An overview of the package and it's functions can be found by running
```R
?insect
```
If you experience a problem using this package please feel free to
raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/insect/issues).

### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.


