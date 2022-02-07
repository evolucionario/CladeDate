# CladeDate
CladeDate: Empirical calibration density generator for divergence time estimation

CladeDate is an R package for the generation of empirical calibration information from the fossil record. CladeDate uses simple mathematical models to estimate the age of a clade and its uncertainty based on fossil ages. Using a Monte Carlo approach, CladeDate generates empirical densities representing the uncertainty associated with the age of the clade and fits standard probability density functions that can be used in time-tree inference software such as BEAST2, MrBayes, and MCMCtree.

# Instalation

From the R console:

library(devtools)
install_github("evolucionario/CladeDate")
library(CladeDate)
