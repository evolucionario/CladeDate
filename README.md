### CladeDate: Empirical calibration density generator for divergence time estimation

Santiago Claramunt

Department of Natural History, Royal Ontario Museum, and
Department of Ecology and Evolutionary Biology, University of Toronto, Ontario, Canada.

E-mail: claramunt.bio@gmail.com

CladeDate is an R package for the generation of empirical calibration information from the fossil record. CladeDate uses simple mathematical models to estimate the age of a clade and its uncertainty based on fossil ages. Using a Monte Carlo approach, CladeDate generates empirical densities representing the uncertainty associated with the age of the clade and fits standard probability density functions that can be used in time-tree inference software such as BEAST2, MrBayes, and MCMCtree.

## ATTENTION users of MCMCtree:
Note that the parameterization of the skew-student distribution provided by clade.date is different from the one used in MCMCtree, thus you can't use those parameters directly when setting a prior. I will fix this in the next version. 

### Instalation

From the R console:

````
library(devtools)
install_github("evolucionario/CladeDate")
library(CladeDate)
````

### References

Claramunt, S. 2022. CladeDate: calibration information generator for divergence time estimation. Methods in Ecology and Evolution  [Early View](https://doi.org/10.1111/2041-210X.13977).
