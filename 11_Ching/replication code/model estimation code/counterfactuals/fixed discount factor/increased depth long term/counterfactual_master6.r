# master file to run counterfactuals

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

cf.index = 6 # index of counterfactual

spec.index = 3  # index of estimation specification

# screen width

options("width"=200)

# paths to data, code files, and output
datapath = ""
fpath = ""

set.seed(821511)

# load in R functions

source(paste(fpath,"ijc_functions_new.r",sep=""))

source(paste(fpath,"estimation_setup",spec.index,".r",sep=""))

source(paste(fpath,"run_cf_data",cf.index,".r",sep=""))

