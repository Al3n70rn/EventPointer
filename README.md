# EventPointer: An effective identification of alternative splicing events using junction microarrays
Software to identify alternative splicing events using junction microarrays

EventPointer R package can be applied to complex experimental designs by giving the required contrast and design matrices.
The software provides a list of the detected splicing events indicating: gene name, type of event (Cassette, Alternative 3',...,etc),
genomic position, statistical significance and affected protein domains.

The algorithm requires low amounts of RAM memory and performs the test in less than 1 minute (Depends on the PC)

## Installation
R package EventPointer can be installed in R as:
```r
library(devtools)
install_github("jpromeror/EventPointer")
```
For problems with the installation of "dcGOR" R Package
refer to:
http://dcgor.r-forge.r-project.org/install.html


## Vignette
The EventPointer Vignette can be visualized in the following URL
https://rawgit.com/jpromeror/EventPointer/master/EP_Vignette.html
