# master file to create our data

rm(list=ls())

library(Rcpp)

# read in store and panel data from the original data.
# if either are FALSE, then just load them
readstore = TRUE
readpanel = TRUE

# run the first step of making the estimation data - this is
# a little slow
makeestdata1 = TRUE

# product category
product = "laundet"

# paths to data files
rawdatapath = #Put the path to the raw IRI files here  #"C:/Users/MattOsborne/Documents/iri data/Academic Dataset External/"
outputpath = #Put the path where the estimation files will be stored here #"C:/Users/MattOsborne/Dropbox/discount factor work/data/"

# List describing when, and what gets dropped

droplist = vector("list",length=1)
dropindx = 1   # index of droplist

# function to update droplist
update.droplist <- function(data,reason,olddata,newdata,idvar) {
  nid1 = length(unique(olddata[[idvar]]))
  nid2 = length(unique(newdata[[idvar]]))
  d1 = droplist
  di1 = dropindx
  droplist[[di1]] = list(data=data,reason=reason,
          nobs1=nrow(olddata),nobs2=nrow(newdata),
          nid1=nid1,nid2=nid2)
  di1=di1+1
  assign("droplist",droplist,envir=globalenv())
  assign("dropindx",di1,envir=globalenv())
}

# variables that are relevant when we create the estimation data

# which demographic variables do we keep in the merged panel
# original dataset names
demonames.orig = c("Combined.Pre.Tax.Income.of.HH","Family.Size",
"HH_RACE","HH_AGE","HH_EDU","HH_OCC",
"Education.Level.Reached.by.Male.HH","Education.Level.Reached.by.Female.HH",
"Children.Group.Code")
# new names
demovars = c("HH_INC","HH_SIZE","HH_RACE","HH_AGE","HH_EDU","HH_OCC",
"maleedu","fedu","child")

# if this flag is true we will remove powders
# a tabulation of the entire dataset shows very little overlap between
# powder and liquid purchasers
#       0    1
#  0   10  507
#  1 4851   24
remove.powder = TRUE

# if we want to fill in missing prices in hhmerge set to true
fillprices = TRUE

# maximum allowable inventory
maxinv = 300

# percentiles of total quantity purchased outside of which we drop observations
quprobs = c(0,1)

# minimum number of units need to purchase to be in sample
minunits = 5

# maximum number of units in a purchase occasion we will allow
bigJ = 5

# maximum acceptable gap in length between purchases

maxpgap = 40

# years we will use for estimation
estyears = 5:7

# number of package sizes we include (based on the HH panel).
# starting with the most popular, then next most popular, etc

nsize = 5  #3

# number of brands to include, starting with most popular, etc

nbrand = 20

# bounds on household consumption rates - sometimes we see really high
# or low rates, probably due to missing data or something else strange
# going on

cratebounds =  c(0.25, 2)

# use the c code for aggregation

use.ccode = TRUE

# method for computing prices when multiple store visits occur
# 1 - use min price
# 2 - use max price
# 3 - use random price (equal weights)
# 4 - use random price (higher weight to more visited store)
# 5 - nonrandom prioritizing stores in order of visits
# 6 - prioritize stores where you actually buy

pricemethod = 6

# flag to load in imputed price arrays
# imputing the prices takes some time, even if it's done in C++

loadprices = FALSE

# rule about price imputation
imprule = 1 # compute the average of prices in periods t+1 and t-1


keepobjs = c(ls(),"keepobjs")

# load store data and store characteristics data
if(readstore) {
  for(bigyear in 1:7) {
    yearnum = paste("Year",bigyear,sep="")
    source("read_storedata.r")
    allobjs = ls()
    rm(list=setdiff(allobjs,c(keepobjs,"yearnum")))
  }
}

# load household panel data, and remove individuals who don't make the static
if(readpanel) {
  for(bigyear in 1:7) {
    yearnum = paste("Year",bigyear,sep="")
    source("read_paneldata.r")
    cat("Number of observations dropped in year ",bigyear," for not making static\n",sep="")
    cat("Observations: ",nobs1-nobs2," of ",nobs1,"\n",sep="")
    cat("IDs: ",nid1-nid2," of ",nid1,"\n\n",sep="")
    allobjs = ls()
    rm(list=setdiff(allobjs,c(keepobjs,"yearnum")))
  }
}

# create estimation data
if(makeestdata1) {
  source("create_estimation_data.r")
} else {
  load(paste(outputpath,"mergedpanel_",product,"_allyears.RData",sep=""))
}

source("create_estimation_data_2.r")


