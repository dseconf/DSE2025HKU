# create initial ml for starting points

rm(list=ls())

# do we do multicore or not
multicore = TRUE

library(randtoolbox)
library(bayesm)
library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
if(multicore) {

  library(RcppParallel)

  setThreadOptions(stackSize = 1024*1024*1024)
}

fit.index = 9  #index of simulation

# screen width

options("width"=200)

# paths to data and code files
datapath = ""
fpath = ""

set.seed(821511)

# load in R functions

source(paste(fpath,"ijc_functions_new.r",sep=""))

source(paste(fpath,"fit_setup",fit.index,".r",sep=""))

source(paste(fpath,"run_fit_data",fit.index,".r",sep=""))
