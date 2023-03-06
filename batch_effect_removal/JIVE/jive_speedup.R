library(tictoc)
library(r.jive)
library(RSpectra)
library(Rcpp)

sourceCpp("matrix_multiplication.cpp", showOutput = FALSE)

####################################################################################################

source("jive_v2.R")
source("jive_iter.R")
source("jive_perm.R")
source("jive_bic.R")
