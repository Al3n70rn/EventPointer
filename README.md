# EventPointer
R package to identify alternative splicing events using junction microarrays

EventPointer R package can be applied to complex experimental designs by giving the required contrast and design matrices.
The software provides a list of the detected splicing events indicating: gene name, type of event (Cassette, Alternative 3',...,etc),
genomic position, statistical significance and affected protein domains.

The algorithm requires low amounts of RAM memory and performs the test in less than 1 minute (Depends on the PC)

The R package is available on [CRAN] and can be installed in R as:
```r
install.packages('EventPointer')
```
[CRAN]: https://cran.r-project.org/
