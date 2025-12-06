# code to fit purchase hazard/purchase probabilities

rm(list=ls())

# do we do multicore or not
multicore = TRUE

library(randtoolbox)
library(bayesm)
library(Rcpp)
library(RcppArmadillo)
#library(numDeriv)
library(MASS)
library(Matrix)
if(multicore) {

  library(RcppParallel)

  setThreadOptions(stackSize = 1024*1024*1024)
}

fit.index = 3  #index of simulation
clust.index = 1  #index of cluster we run
nclust = 3  # number of clusters

# screen width

options("width"=200)

# paths to data, code files, and output

datapath = ""
fpath = ""

set.seed(821511)

# load in R functions

source(paste(fpath,"ijc_functions_new.r",sep=""))

source(paste(fpath,"fit_setup_clust",fit.index,".r",sep=""))

source(paste(fpath,"run_fit_data_clust",fit.index,".r",sep=""))
