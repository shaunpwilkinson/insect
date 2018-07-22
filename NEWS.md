# insect 1.1.1

Patch release addressing a bug in `classify` causing rearrangement
of columns in the output table. Thanks to Rachel Wade for the report.

Also, the `classify` function no longer returns factors in the output 
table.



# insect 1.1.0

Minor release with following improvements: 

* Taxonomy information now added to top level of classification trees.
* Prevented conversion to cladogram following tree-learning step to save time.
* searchGB now appends taxonomy IDs to sequence names.
* Sequence naming change in line with qiime pipeline format.
* classify now outputs a data frame.
* Obsoleted tabulize function.

# insect 1.0.0

Submitted to CRAN 2018-05-13.
