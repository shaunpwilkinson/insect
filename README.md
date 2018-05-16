# insect


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


