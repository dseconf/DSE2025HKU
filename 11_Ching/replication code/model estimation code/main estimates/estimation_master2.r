# master file to run estimation

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

spec.index = 2  # index of estimation specification

# screen width

options("width"=200)

# paths to data, code files, and output
datapath = ""
fpath = ""

set.seed(821511)

# load in R functions

source(paste(fpath,"ijc_functions_new.r",sep=""))

source(paste(fpath,"estimation_setup",spec.index,".r",sep=""))

source(paste(fpath,"run_initial_ml_data",spec.index,".r",sep=""))

source(paste(fpath,"run_ijc_data",spec.index,".r",sep=""))

