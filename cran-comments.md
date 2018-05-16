# insect version 1.0.0

This is the second CRAN submission addressing the following comments from Swetlana:

Thanks, please add a reference for the method in the 'Description' field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.

Please replace \dontrun{} by \donttest{} in your Rd-files.

Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace. That is not allow by CRAN policies. Please only write/save files if the user has specified a directory. In your examples/vignettes/tests you can write to tempdir().

Please fix and resubmit.

Best,
Swetlana Herbrandt


Thank you very much for the feedback, these points have all now been remedied.

## Test environments

 * local ubuntu 16.04.4 x86_64-pc-linux-gnu; R version 3.4.4
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.5.0
 * winbuilder devel R version using R version 3.5.0 beta (2018-04-10 r74581)

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Shaun Wilkinson <shaunpwilkinson@gmail.com>'

Possibly mis-spelled words in DESCRIPTION:
  Informatic (3:8)
  barcoding (9:62)
  bioinformatics (9:25)
  de (10:47)
  informatic (11:46)

These have been checked and are all ok. 

## Downstream dependencies

None.
