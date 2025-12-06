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

desktop = TRUE # windows desktop or laptop

cf.index = 1 # index of counterfactual

spec.index = 2  # index of estimation specification (see word document)

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

