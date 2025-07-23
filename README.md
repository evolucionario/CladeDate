## CladeDate: Empirical calibration density generator for divergence time estimation

### Version: 1.2

Santiago Claramunt

Department of Biological Sciences, University of New Orleans, New Orleans, Louisiana, USA.

E-mail: sclaramu@uno.edu

CladeDate is an R package for the generation of empirical calibration information from the fossil record. CladeDate uses fossil ages together with estimators of upper bounds of truncated distributions to estimate the age of clades and their uncertainty. Using a Monte Carlo approach, CladeDate generates empirical densities representing the uncertainty associated with the age of the clade and fits standard probability density functions that can be used in time-tree estimation software such as BEAST2, MrBayes, and MCMCtree.

## Instalation

From the R console, first install required packaged if they are not installed already:
```
install.packages("devtools")
install.packages("phangorn")
install.packages("sn")
```
Then istall CladeDate from GitHub:
```
library(devtools)
install_github("evolucionario/CladeDate")
library(CladeDate)
```

## References

Claramunt, S. (2022) CladeDate: calibration information generator for divergence time estimation. Methods in Ecology and Evolution 13, 2331-2338. https://doi.org/10.1111/2041-210X.13977
