# Setup parameters for estimation on IRI data

################################################################################

# load ML estimates with ML covariance matrix

loadml = FALSE

# which code to compute hessian.  Slow code is more accurate

slowhess = FALSE

# set to TRUE if we are calculating DIC

computeDIC = FALSE

# number of iterations for DIC (set to 1 for myopic model and more for fwd-looking to guarantee convergence)

niter.DIC = 10

# spec where we only estimation the brand specific coefficients.

brandonly = FALSE

# simulate brand choices (when brandonly = TRUE)

brandsim = FALSE

# estimate the model with heterogeneous parameters? (deprecated??)

h.params = TRUE

# set to TRUE to toggle using consumption rate from data, rather than calibrating it

data.crate = TRUE

# set to TRUE if we want to enforce a lower bound on consumption rate that is from the data

crhhlb = FALSE

# set to true if we add in shifters for bigger sizes (put at end of par)

sizeshifter = TRUE

# set to true if we use consumption rate in initial estimation

init.crate = FALSE

# redefine package sizes in terms of estimated laundry load
# size, and round up so that we fit the state space exactly

redefsize = TRUE

# true for simulated data

simulation = FALSE

# enforce myopic model

myopic = FALSE

# true if we do reduced form tests - this takes memory so don't do it for estimation runs
rform = FALSE

# load stage 1 or not

loadfirststage = FALSE

# enforce zero storage costs up to capacity constraint

zeroscost = TRUE

# number of demographic vars we keep

ndvars = 1

if(rform) {
  library(plm)
}

# if we want to fill in missing prices set to true
fillprices = TRUE

# true if we are doing em algo

do.em = FALSE

# true if we are doing ijc

do.ijc = TRUE

# true if we are doing icj and stage 1 is estimated

do.ijc.stage1 = TRUE

# number of cut points in inclusive value spline
# for 0 we use a linear regression
ncutiv = as.integer(1)

# type of spline we are using
# 1: linear spline (deprecated)
# 2: b spline
# 3: b spline, but only with 1 size regressed on its own basis
# note if ncutiv = 0, the spline type doesn't matter
splinetype = as.integer(2)

# this is only true for debugging/artificial data - it will set the
# consumption shock draws to be predetermined
# if true cdraws will come from info$cdraws
usecdraws = FALSE

# integration type  (I think this is deprecated)
# 1 - quadrature
# 2 - monte carlo

inttype = as.integer(2)

# method for drawing random grid (stored in Mcmc)
# 1: grid over inclusive values, is randomly drawn each iteration
# 2: grid over prices, is drawn in advance and fixed

gridmethod = as.integer(2)

# deprecated - was used when we did the polynomial approximation
necoef = as.integer(0)

# homogeneous consumption rate

crfix = TRUE

# force homogeneity across store visits, inclusive value transitions
# use this only when all parameters are fixed across the population
# this is much faster to estimate

h.model = FALSE

# number of package sizes we include (based on the HH panel).
# starting with the most popular, then next most popular, etc

nsize = 5

# number of brands to include, starting with most popular, etc
# this comes out of the size of brsize
# nbrand = 20

# max number of units we will allow
# this is now read in from the data file since the dropping of data occurs there
# bigJ = 5

# number of grid points for IV interpolation
niv = 10

# inventory model
# 1 - our initial model with bound on inventory
# 2 - the FIFO model
invmodel = as.integer(1)
maxbottles = as.integer(4)  #number of bottles you can hold in the FIFO model

# simulation initial choices (only relevant in choicesimijc so we can set to false here)

siminitchoices = FALSE

# show plots of inventory states and value functions

doplots = FALSE

# round prices to integer - can help with state space

intprices = TRUE

# number of grid points for inventory interpolation
if(invmodel == 1) {
  ninv = 100
} else {
  if(redefsize) {

    invd1 = 0.2
    invd2 = 1

    break1 = 1
    break2 = 60*0.2

    ninv = 1 + break1/invd1 + (break2-break1)/invd2

  } else {
    invd1 = 0.125 #0.25
    invd2 = 0.25 #0.5

    break1 = 3.25
    break2 = 12.5

    ninv = 1 + break1/invd1 + (break2-break1)/invd2  #can we automate this?

  }

}

# this is where we put a break in inventory interpolation
if(invmodel==1) {
  first20 = as.integer(0.2*ninv)
} else {
  first20 = as.integer(1+break1/invd1) #as.integer(27)
}
# number of simulation draws in second step MLE (deprecated)
if(crfix) {
  nsim = 1
} else {
  nsim = 500 #20
}

# how we draw initial inventories. 1 for using data and 2 for drawing choices
idrawtype = as.integer(1)

# number of periods to draw initial inventories
ninitt = as.integer(52)

# storage cost model type
if(invmodel == 1) {
  contmodel = as.integer(2)
} else {
  contmodel = as.integer(1)
}
#0: discrete quadratic (default case)
#1: continuous quadratic
#2: linear starting at omega0

# TRUE if we want to estimate inventory bounds (and assume 0 inventory)

if(contmodel==2) {
  h.invbound = TRUE
} else {
  h.invbound = FALSE
}

# error distribution type... not yet finished
errortype = 1
#1: standard logit for all choices
#2: consumption shock on buy/not buy decision and logit everywhere else

# maximum number of switches we record if errortype > 1
maxswitch = as.integer(10)

# flag for heterogeneous parameters
hparams = TRUE

# number of random grid draws - used for heterogeneous parameters model
nrgrid = as.integer(100)

# generate random draws within likelihood evaluation.  TRUE will be faster
genrnd = TRUE
if(do.ijc) {
  genrnd = FALSE
}

# return individual inventories in estimation (should this be deprecated??)
retinv = FALSE

# debug flag - not really used now but might be useful later
debug = FALSE

# product

product = "laundet"

# integration type
# 1 - quadrature
# 2 - monte carlo

inttype = as.integer(2)

# do iv dependence regressions on fixed grid, rather than data
# I think this is deprecated since we don't use the dependence grid anymore
doivdepgrid = TRUE
pctdg = 0.1  # percent of data we will keep

# maximum allowable inventory
maxinv = 300

# setup up number of cost draws, quad points, etc.
# I think this may not be necessary anymore except for checking

if(crfix) {
  ncost = 1
} else {
  ncost = 10

  quad <- .Fortran("gaussq", as.integer(1), as.integer(ncost), as.double(0), as.double(0),
  as.integer(0), as.double(c(-1, 1)), double(ncost), t = double(ncost), w = double(ncost), PACKAGE = "gss")
}

if(inttype == 1) {
  nqpts = 15

  quadiv <- .Fortran("gaussq", as.integer(4), as.integer(nqpts), as.double(0), as.double(0),
  as.integer(0), as.double(c(-1, 1)), double(nqpts), t = double(nqpts), w = double(nqpts), PACKAGE = "gss")

  qpts = cbind(quadiv$t*sqrt(2),quadiv$w/sqrt(pi))

} else {

  nqpts = 100

  hdraws = halton(nqpts+100,dim=nsize)

  qpts = qnorm(hdraws[100+1:nqpts,])

}


# data stuff - may not be needed so check it

yearlist = paste("Year",1:7,sep="")

weekstub = c("1114_1165","1166_1217","1218_1269","1270_1321","1322_1373","1374_1426","1427_1478")

ystub = c("1","2","3","4","5","6","7")

estyears = c(5:7)    # years we will use for estimation
simyears = c(1:4)    # years we will use for simulation
allyears = c(simyears,estyears)

load(paste(datapath,"estimationdata_",product,"_allyears.RData",sep=""))

# cut down data to remove households that disappear for a while

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]
lpanid = c(0,hhmerge$PANID[1:(nrow(hhmerge)-1)])
first = hhmerge$PANID != lpanid

hhcrate = hhmerge$crate[first]

hhinitmerge = hhinitbig[!is.na(hhinitbig$ppunit),]

stvisitinit = !is.na(hhinitbig$ppunit)

vol1 = hhmerge$hhvol[first]
initvolhh = aggregate(hhinitmerge$UNITS,by=list(hhinitmerge$PANID),FUN=sum)

firstweek = aggregate(hhinitmerge$WEEK,by=list(hhinitmerge$PANID),FUN=min)
lastweek = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=max)

hh1 = data.frame(cbind(hhmerge$PANID[first],vol1,initvolhh$x,firstweek$x,lastweek$x,hhcrate))
colnames(hh1) = c("PANID","vol","initvol","firstweek","lastweek","crate1")

hh1$crate2 = (hh1$vol+hh1$initvol)/(hh1$lastweek-hh1$firstweek+1)

hhcrate = hh1$crate2
#hhcrate = hhmerge$crate[first]

inv1 = getinv(hhinitbig,hhbig,hh1$crate1,brsize,vols)

inv2 = getinv(hhinitbig,hhbig,hh1$crate2,brsize,vols)

#hhcrate2 = getnewcrate(hhcrate,inv1$initx,inv1$xbig) #recalibrate cons rates
#hhcrate = hhcrate2

maxnzero = apply(inv2$nzero,2,max,na.rm=TRUE)  # this seems to give slightly more sensible inventory dist

allids = unique(hhinitbig$PANID)

allids = allids[order(allids)]

maxinv2 = apply(inv2$invbig,2,max,na.rm=TRUE)

avginv2 = apply(inv2$invbig,2,mean,na.rm=TRUE)
avginv2pre = apply(inv2$initinv,2,mean,na.rm=TRUE)

temp = abs(avginv2/6.25-avginv2pre/6.25) > 1

crate.new1 = hhcrate
crate.new1[temp] = hh1$crate1[temp]

inv.new1 = getinv(hhinitbig,hhbig,crate.new1,brsize,vols)
maxinv.new1 = apply(inv.new1$invbig,2,max,na.rm=TRUE)

keep1 = maxnzero <= 10 & maxinv.new1/6.25 < 15
nclust = 3
crclust = hclust(dist(crate.new1[keep1]),method="ward.D2")
crassign = cutree(crclust,nclust)

temp = cumsum(keep1)
#temp = temp[keep1]

crassign = crassign[-c(temp[136],temp[399])]

keep1[c(136,399)] = FALSE #this guy caused problems

idkeep = allids[keep1]  #missing brindex 4, 14, I think is ok

hhbig = hhbig[hhbig$PANID %in% idkeep,]
hhmerge = hhmerge[hhmerge$PANID %in% idkeep,]
hhmergebr = hhmergebr[hhmergebr$PANID %in% idkeep,]

hhinitbig = hhinitbig[hhinitbig$PANID %in% idkeep,]
hhinitmerge = hhinitmerge[hhinitmerge$PANID %in% idkeep,]
#stop()
nbrand = nrow(brsize)

hhbig = hhbig[order(hhbig$PANID,hhbig$WEEK),]
hhmergebr = hhmergebr[order(hhmergebr$PANID,hhmergebr$WEEK),]
hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]
lpanid = c(0,hhmerge$PANID[1:(nrow(hhmerge)-1)])

first = hhmerge$PANID != lpanid

hhcrate = crate.new1[keep1]

cfactor = 16/3.125

cat("Consumption rate distribution (New Sample):\n")

print(summary(hhcrate*cfactor))

cat("Number of households in sample: ",length(hhcrate),"\n")

# starting parameter values - some of these will be based on data

# consumption rate (lower and upper bounds)
if(crfix) {
  crate = rep(mean(hhmerge$crate[first]),2)
} else {
  crate = c(min(hhmerge$crate),max(hhmerge$crate))
}

# price coefficient
alfa = -0.5

# consumption utility - normalize to zero
gama = 0

# stockout cost
nu = 1

# discount factor
if(myopic) {
  beta = 0
} else {
  beta = 0.95
}

# logit variance - normalize
logitf = 1

# fixed cost of purchase
fcost = -2

# probability of store visit - assume fixed across households for now
stvisitprob = 1-sum(hhbig$totunits < 0)/nrow(hhbig)

# storage costs
if(contmodel==2) {
  omega = c(299,0)
} else {
  omega = c(0.1,0.01)
}
# starting dynamic parameters
xdynamic = c(crate,alfa,gama,nu,beta,logitf,fcost,stvisitprob,omega)

# package size
packsize = vols[1:nsize]

if(redefsize & invmodel == 2) {

  convfactor = 0.2
  packsize1 = packsize/convfactor
  packsize1 = (invd2/convfactor)*as.integer(packsize1/(invd2/convfactor))
  packsize = packsize1*convfactor

}

# I think ncut and ncutiv are redundant here - they were different when
# we did the IV dependence grid

nhhs = as.integer(length(unique(hhmerge$PANID)))
nobs = nrow(hhbig)

# add in a variable to track bottle choice to init

hhinitbig$sizechoice = as.integer(0)
hhinitbig$nbottles = as.integer(0)
hhinitbig$multi = 0
for(i in 1:nsize) {
  x = hhinitbig[[paste("Volunits",i,sep="")]]>0
  hhinitbig$sizechoice[x] = as.integer(i)
  hhinitbig$nbottles[x] = as.integer(hhinitbig[[paste("Volunits",i,sep="")]][x])
  hhinitbig$multi[x] = hhinitbig$multi[x] + 1
}

# check how many sizes are purchased in the data

sizechoice = brsize[hhmergebr$brindex,2]
nsizehh1 = aggregate(sizechoice,by=list(hhmergebr$PANID,sizechoice),FUN=length)
colnames(nsizehh1) = c("PANID","size","nobs")
nsizehh = aggregate(nsizehh1$nobs,by=list(nsizehh1$PANID),FUN=length)
tobshh = aggregate(nsizehh1$nobs,by=list(nsizehh1$PANID),FUN=sum)
colnames(tobshh) = c("PANID","tobshh")
nsizehh1 = merge(nsizehh1,tobshh)
nsizehh1$frac = nsizehh1$nobs/nsizehh1$tobs
nsizehh2 = aggregate(nsizehh1$frac*nsizehh1$frac,by=list(nsizehh1$PANID),FUN=mean)
print(table(nsizehh$x))/nhhs
if(doplots) {hist(nsizehh2$x)}
nsizehh3 = aggregate(nsizehh1$frac,by=list(nsizehh1$PANID),FUN=max)
if(doplots) {hist(nsizehh3$x)}

# construct elements of info list

# put constants into info

info = list(nbrand = as.integer(nrow(brsize)), prindex = as.integer(5), print = FALSE,
bigJ = as.integer(bigJ), ncut = as.integer(ncutiv),
ngrid = as.integer(c(rep(niv,nsize),ninv)),
ncutiv = as.integer(ncutiv), retindll = FALSE, nsize = as.integer(nsize),
nhhs = as.integer(nhhs), crfix=crfix, data.crate=data.crate, usecdraws = usecdraws,
inttype = as.integer(inttype), nrgrid = as.integer(nrgrid), necoef=as.integer(necoef),
nsim = as.integer(nsim), retinv = retinv, myopic = myopic,
contmodel = as.integer(contmodel), debug = debug, idrawtype = as.integer(idrawtype),
first20 = as.integer(first20), h.invbound = h.invbound, genrnd = genrnd,
hhprint = as.integer(0), repprint = as.integer(0), ncost = as.integer(ncost),
qpts = qpts, maxinv = maxinv, h.model = h.model, packsize = packsize, tol = 1e-5,
brsize = brsize, splinetype = splinetype, crhhlb = crhhlb, brandonly = brandonly,
siminitchoices = siminitchoices, sizeshifter = sizeshifter)

info$brsize = matrix(as.integer(brsize),nrow=nrow(brsize),ncol=ncol(brsize))

# make matrix of demographics

info$dvarnames = c("Constant","highincome","agedummy","college","fsize2")

info$dvars = as.matrix(cbind(rep(1,nhhs),hhmerge[first,info$dvarnames[-1]]),
                       nrow=nhhs,ncol=length(info$dvarnames))

colnames(info$dvars) = info$dvarnames

for(i in 2:length(info$dvarnames)) {
  cat("table of demographic variable ",info$dvarnames[i],"\n")
  print(table(info$dvars[,i]))
}

info$invmodel = invmodel
info$maxbottles = maxbottles

# bstates - later we may use this
info$bstates = as.integer(1)
info$revarray = as.integer(1)
info$nb = as.integer(1)

# to do :

# state space grid
# the first 3 rows don't matter, they've been replaced with rgrid
info$pgrid = matrix(0,nrow=nsize+1,ncol=max(nrgrid,ninv))

posinv = hhmerge$invnew[hhmerge$invnew > 0]

pr = 0.2
qinv = quantile(posinv,probs=pr)
while(qinv == 0) {
  pr=pr+0.1
  qinv = quantile(posinv,probs=pr)
}

# compute and index inventory states
if(invmodel == 1) {
  # old version
  info$pgrid[nsize+1,1:first20] = qinv*seq(0,first20-1)/(first20-1)
  info$pgrid[nsize+1,(first20+1):ninv] = qinv + (maxinv-qinv)*seq(1,ninv-first20)/(ninv-first20)
} else {
  #info$pgrid[nsize+1,1:ninv] = 0:(ninv-1)*0.25  #think about whether or not we should start at 0.
  info$pgrid[nsize+1,1:ninv] = c(seq(0,break1,invd1),seq(break1+invd2,break2,invd2))
  # define all the b states
  temp = getbstates(nsize,maxbottles)
  info$nb = temp[[1]]
  info$bstates = temp[[2]]
  info$revarray = temp[[3]]
  #info$first20 = as.integer(ninv)  #this should force the interpolation to use an evenly spaced grid

}

# data containing initial inventories

info$initx = hhinitbig$UNITS
info$initb = hhinitbig$nbottles
info$inits = hhinitbig$sizechoice
if(idrawtype==1) {
  ninitt = length(unique(hhinitbig$WEEK))
  info$ninitt = ninitt
}

# consumption draws and initial consumption draws if debugging

if(genrnd) {
  info$cdraws = as.double(0)
} else {
  info$cdraws = runif(nobs*nsim)
}

if(genrnd) {
  info$initcdraws = as.double(0)
} else {
  info$initcdraws = runif(nsim*nhhs*ninitt)
}

# these don't get used as far as I can tell - remove them?
info$initivdraws = as.double(0)
info$logitdraws = as.double(0)


# gmean, gvari, dg - these only matter if gridmethod is 1
# put in filler values and fill in later, using a preliminary estimate

info$gmean = rep(0,nsize)
info$gvari = diag(nsize)
info$dg = 1.0

# this is redefinition/resimulation of the data

if(intprices) {

  hhmerge[,colnames(hhmerge)[1:info$nbrand + info$prindex - 1]] =
    round(hhmerge[,colnames(hhmerge)[1:info$nbrand + info$prindex - 1]],2)
  hhmergebr[,colnames(hhmerge)[1:info$nbrand + info$prindex - 1]] =
    round(hhmergebr[,colnames(hhmerge)[1:info$nbrand + info$prindex - 1]],2)

}

#hhcrate = hhmerge$crate[first]

if(data.crate) {
  info$crate = hhcrate
}

if(multicore) {
  sourceCpp("ijc_Rcpp_par.cpp")
} else {
  sourceCpp("ijc_Rcpp.cpp")
}

# indexes of observations for each household in the big and merged data

beginobs = aggregate(1:nrow(hhbig),by=list(hhbig$PANID),FUN=min)

info$obshhinds = beginobs$x

beginobs1 = aggregate(1:nrow(hhmerge),by=list(hhmerge$PANID),FUN=min)

info$obshhindsmerge = beginobs1$x

# expansion vectors

#expansion vector for hhmerge to hhbig

id1 = 1:nrow(hhmerge)

d1 = data.frame(hhmerge$WEEK,hhmerge$PANID,id1)

colnames(d1) = c("WEEK","PANID","id1")

d2 = merge(hhbig,d1,all.x=TRUE)

d2$id1[is.na(d2$id1)] = 1

info$expandbig = d2$id1

# expansion vector for hh vector to hhbig - for a/r

hhids = unique(hhbig$PANID)

id1 = 1:length(hhids)

d1 = data.frame(hhids,id1)

colnames(d1) = c("PANID","id1")

d2 = merge(hhbig,d1,all.x=TRUE)

info$expandhh = d2$id1

# expansion vectory for hh vector to hhmerge

d2 = merge(hhmerge,d1,all.x=TRUE)

info$expandmerge = d2$id1

rm(id1,d1,d2)

# default - all states considered

if(invmodel == 1) {

  badmissable = TRUE

  info$badmissable = badmissable
  info$nbsmall = as.integer(sum(badmissable))
  info$bindex = as.integer(cumsum(badmissable))

} else {

  badmissable = rep(TRUE,ncol(info$bstates))

  info$badmissable = badmissable
  info$nbsmall = as.integer(sum(badmissable))
  info$bindex = as.integer(cumsum(badmissable))
}

if(invmodel == 2) {
  # check which b states get visited


  if(nsize == 2) {

    #c(1,2,3,4,7,8,15)
    badmissable = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,
                    FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)  # all states with only 1 size

    info$badmissable = badmissable
    info$nbsmall = as.integer(sum(badmissable))
    info$bindex = as.integer(cumsum(badmissable))
    info$bindex[!badmissable] = -1

  } else if(nsize == 3) {

    badmissable = rep(FALSE,ncol(info$bstates))
    badmissable[c(1:5,9,13,14,27,40)] = TRUE

    info$badmissable = badmissable
    info$nbsmall = as.integer(sum(badmissable))
    info$bindex = as.integer(cumsum(badmissable))
    info$bindex[!badmissable] = -1

  } else {

    hhinitmerge = hhinitbig[!is.na(hhinitbig$ppunit),]

    stvisitinit = !is.na(hhinitbig$ppunit)

    #iiout = checkinvstate(hh1$crate2,hhbig,info,as.integer(c(1,nhhs)),
    #                      hhinitmerge,stvisitinit)

    iiout = checkinvstate(hhcrate,hhbig,info,as.integer(c(1,nhhs)),
                          hhinitmerge,stvisitinit)

    # figure out total volume held - for check against data

    tvol = iiout$inv
    for(i in 2:maxbottles) {
      inds = info$bstates[i,iiout$b]>0
      tvol[inds] = tvol[inds] + info$packsize[info$bstates[i,iiout$b[inds]]]
    }
    inv1 = getinv(hhinitbig,hhbig,hhcrate,brsize,vols)

    hhcrate2 = getnewcrate(hhcrate,inv1$initx,inv1$xbig)

    inv2 = getinv(hhinitbig,hhbig,hhcrate2,brsize,vols)

    # note: hh1$hhcrate2 and hhcrate2 aren't hugely different, second reduces state space a lot tho
    #hhcrate = hhcrate2
    #info$crate = hhcrate2
    #iiout = checkinvstate(hhcrate2,hhbig,info,as.integer(c(1,nhhs)),
    #                      hhinitmerge,stvisitinit)

    save(iiout,file=paste("iiout",maxbottles,".RData",sep=""))

    #iiout = checkinvstate(rep(mean(hhcrate),nhhs),hhbig,info,as.integer(c(1,nhhs)),
    #                      hhinitmerge,stvisitinit)

    temp = table(iiout$b)
    temp = temp[order(temp,decreasing=TRUE)]/nrow(hhbig)
    print(cbind(temp))

    stfrac = cumsum(temp)

    print(cbind(1:sum(stfrac<0.99),t(info$bstates)[as.integer(names(temp[stfrac<0.99])),],stfrac[stfrac<0.99]))

    # get a sense of how often we see full inventory states

    tt = t(info$bstates)
    tt1 = tt[as.integer(names(temp)),]
    print(sum(temp[tt1[,maxbottles]>0]))

    print(sum(hhbig$totunits > 0 & tt[iiout$b,maxbottles] > 0)/nrow(hhbig))

    # comparisons
    if(1==0) {
      iiout3 = iiout
      for(i in 4:5) {
        load(paste("iiout",i,".RData",sep=""))
        print(summary(abs(iiout$inv-iiout3$inv)))
        print(sum(iiout$inv==iiout3$inv)/nrow(hhbig))
        inds = iiout$inv!=iiout3$inv
        print(summary(abs(iiout$inv[inds]-iiout3$inv[inds])*16))
      }
      load("iiout4.RData")
      iiout4 = iiout
      load("iiout5.RData")
      print(summary(abs(iiout$inv-iiout4$inv)))
      print(sum(iiout$inv==iiout4$inv)/nrow(hhbig))
      inds = iiout$inv!=iiout4$inv
      print(summary(abs(iiout$inv[inds]-iiout4$inv[inds])*16))
      stop()
    }

    badmissable = rep(FALSE,info$nb)
    #badmissable[1:6] = TRUE
    badmissable[as.integer(names(temp[stfrac<0.95]))] = TRUE

    info$badmissable = badmissable
    info$nbsmall = as.integer(sum(badmissable))
    info$bindex = as.integer(cumsum(badmissable))
    info$bindex[!badmissable] = -1

  }

}


# construct importance grid for prices.  This is only
# used for gridmethod == 2

bylist = NULL

i=0
for(v in colnames(hhmerge)[1:info$nbrand + info$prindex - 1]) {
  i=i+1
  bylist[[i]] = hhmerge[[v]]
}

pricetab = aggregate(hhmerge$PANID,by=bylist,FUN=length)

colnames(pricetab) = c(colnames(hhmerge)[1:info$nbrand + info$prindex - 1],"Freq")

pricetab$prob = pricetab$Freq/sum(pricetab$Freq)

psampleinds = sample(1:nrow(pricetab),size=nrgrid,replace=TRUE,prob=pricetab$prob)

# note that the price grid is different here.
# the ROWS should have prices, and the COLUMNS should index draws

info$rgrid = t(as.matrix(pricetab[psampleinds,1:info$nbrand]))

info$impdist = pricetab$prob[psampleinds]

# make indexes of all the prices
# indexes should be at household level

bylist1 = list(hhmerge$PANID)
i=0
for(v in colnames(hhmerge)[1:info$nbrand + info$prindex - 1]) {
  i=i+1
  bylist1[[i+1]] = hhmerge[[v]]
}

pricetab2 = aggregate(hhmerge$PANID,by=bylist1,FUN=length)
colnames(pricetab2) = c("PANID",colnames(hhmerge)[1:info$nbrand + info$prindex - 1],
        "nobspr")
idtab = aggregate(pricetab2$PANID,by=list(pricetab2$PANID),FUN=length)
colnames(idtab) = c("PANID","nobsid")

pricetab2 = merge(pricetab2,idtab)
pricetab2 = pricetab2[order(pricetab2$PANID),]
pricetab2$prnum = 1:nrow(pricetab2)

idlag = c(0,pricetab2$PANID[2:length(pricetab2$PANID)-1])
idfirst = pricetab2$PANID != idlag

idpriceind = pricetab2$prnum-cumsum(idfirst*pricetab2$nobsid)+pricetab2$nobsid

pricetab2$idpriceind = idpriceind

temp = merge(hhmerge,pricetab2,all.x=TRUE)

temp = temp[order(temp$PANID,temp$WEEK),]

temp = temp[,c("PANID","WEEK","idpriceind")]

temp1 = merge(hhbig,temp,all.x=TRUE)

temp1 = temp1[order(temp1$PANID,temp1$WEEK),]

temp1$idpriceind[is.na(temp1$idpriceind)] = 0

info$idpriceindmerge = temp$idpriceind
info$idpriceindbig = temp1$idpriceind

info$maxpind = as.integer(max(info$idpriceindmerge))

pricetab$prindex = 1:nrow(pricetab)

temp = merge(hhmerge,pricetab,all.x=TRUE)
temp = temp[order(temp$PANID,temp$WEEK),]
temp = temp[,c("PANID","WEEK","prindex")]

temp1 = merge(hhbig,temp,all.x=TRUE)
temp1 = temp1[order(temp1$PANID,temp1$WEEK),]
temp1$prindex[is.na(temp1$prindex)] = 0

info$hpriceindbig = as.integer(temp1$prindex)
info$hmaxpind = as.integer(nrow(pricetab))

# set up full parameter vector - some things may be fixed here
if(sizeshifter) {
  info$xfull = c(rep(0,info$nbrand),xdynamic,rep(0,nsize-1))
} else {
  info$xfull = c(rep(0,info$nbrand),xdynamic)
}
npbig = length(info$xfull)
nbrand = info$nbrand

# set up fixed and nonfixed parameters
info$fixed = rep(TRUE,npbig)

# decide which parameters should be set equal to other parms
# here we will set size intercepts on brands to be equal

if(1==1) {
  brsize.l = c(0,brsize[1:(nrow(brsize)-1),1])

  brsize.e = (brsize.l != brsize[,1])*(1:nrow(brsize))

  for(i in 2:nrow(brsize)) {
    if(brsize.e[i] == 0) {
      brsize.e[i] = brsize.e[i-1]
    }
  }

  info$paramequal = as.integer( c(brsize.e*(brsize.l == brsize[,1]),rep(0,length(xdynamic))) )
  if(sizeshifter) {
    info$paramequal = c(info$paramequal,rep(as.integer(0),nsize-1))
  }
} else {
  info$paramequal = as.integer( rep(0,npbig) )
}
# cons rate
info$fixed[nbrand+1] = FALSE
info$fixed[nbrand+2] = crfix

# fix cons rate if it comes from the data
if(data.crate) {
  info$fixed[nbrand + 1:2] = TRUE
}

# price coefficient
info$fixed[nbrand+3] = FALSE

# stockout cost
info$fixed[nbrand+5] = FALSE

# discount factor
info$fixed[nbrand+6] = myopic

# fixed cost of purchase
info$fixed[nbrand+8] = FALSE

# inventory bound
if(contmodel==2) {
  if(h.invbound) {
    #info$fixed[nbrand+10] = FALSE
  }
} else {
  info$fixed[nbrand+10:11] = FALSE
}

if(sizeshifter) {
  info$fixed[(npbig-(nsize-2)):npbig] = FALSE
}
# upper and lower bounds on parameters

#info$lbounds = c(rep(-Inf,nbrand),0,0,-Inf,0,0,0,0,-Inf,0,0,0)
#info$ubounds = c(rep(Inf,nbrand),Inf,Inf,0,Inf,Inf,1,Inf,0,1,info$maxinv,Inf)

crub = Inf #max(c(max(hhmerge$crate),info$packsize))

info$lbounds = c(rep(-Inf,nbrand),0,0,-Inf,0,0,0,0,-Inf,0,0,0)
info$ubounds = c(rep(Inf,nbrand),crub,Inf,0,Inf,Inf,1,Inf,0,1,info$maxinv,Inf)

# parameter transformations
# 0: no transformation
# 1: exponential
# 2: lbound + (ubound-lbound)*exp(x)/(1+exp(x))
if(data.crate) {
  info$tform = c(rep(0,nbrand),c(0,0,-1,0,1,2,0,-1,0,2,1))
} else {
  if(!is.finite(crub)) {
    info$tform = c(rep(0,nbrand),c(1,1,-1,0,1,2,0,0,0,2,1))
  } else {
    info$tform = c(rep(0,nbrand),c(2,1,-1,0,1,2,0,0,0,2,1))
  }
}
if(contmodel==1) {
  lx = length(info$xfull)
  info$tform[c(lx-1,lx)] = c(1,1)
}

if(sizeshifter) {
  info$lbounds = c(info$lbounds,rep(-Inf,nsize-1))
  info$ubounds = c(info$ubounds,rep(Inf,nsize-1))
  info$tform = c(info$tform,rep(0,nsize-1))
}
# simulate artificial data

# parameter setup

param = matrix(rep(info$xfull,nhhs),nrow=length(info$xfull),ncol=nhhs)

#nsave = 10

#ntilde = 3
#indexes = 1:10
#nextind = 1
#vfsave = rep(0,nsave*info$nrgrid*info$nb*info$ngrid[length(info$ngrid)]*nhhs)
#oldrgrid = rep(0,nsize*nrgrid*nsave)
#oldparam = rep(0,(length(info$xfull)+info$nbrand)*nhhs*nsave)

rep = 1

# the parameters we estimate should be un-fixed

# only need to focus on brand pars, the others should be ok

info$fixed[1:nrow(info$brsize)] = TRUE
info$fixed[1:nrow(info$brsize)][info$brsize[,1]>1 & info$paramequal[1:nrow(info$brsize)]==0] = FALSE

info$dvarnames = info$dvarnames[1:ndvars]

if(ndvars == 1) {
  info$dvars = matrix(info$dvars[,1:ndvars],nrow = info$nhhs,ncol=1)
} else {
  info$dvars = info$dvars[,1:ndvars]
}

# free up some of the bigger brands - helps with fits to shares

if(sizeshifter) {
  info$fixed[2] = FALSE
  info$paramequal[2] = as.integer(0)
  info$fixed[6] = FALSE
  info$paramequal[6] = as.integer(0)
  info$fixed[8] = FALSE
  info$paramequal[8] = as.integer(0)
  info$sizebrand = rep(TRUE,info$nbrand)
  info$sizebrand[c(2,6,8)] = FALSE
}

nbcoef = sum(info$brsize[,1]>1 & info$paramequal[1:nrow(info$brsize)]==0)

if(data.crate) {
  info$fixed[info$nbrand+1] = TRUE
}

npsmall = sum(!info$fixed)

# Mcmc contains information relevant to the MCMC sampler
# for the simulation we need to know what params are fixed across
# individuals and which ones aren't

Mcmc = list(R=20000,keep=1,nprint=100,nsave = 1,nsave.mcmc=1000,dostage1=do.ijc.stage1,
rho=rep(0.1^2,2),rho1=rep(0.05^2,2),mhtype=2,npropsave=25,
s1drawtype = 2,mhtype2=2,ntilde=1,splitfixed=FALSE,
gridmethod=as.integer(gridmethod),
nitervf=1,niterstop=1,vfinterpmethod = as.integer(2),nwrite=as.integer(1000))

# set up population-fixed coefficients
# we have to redefine this below I think
Mcmc$fixedcoef = rep(TRUE,npbig)
# unfix brand 1 and the price coefficient
#Mcmc$fixedcoef[c(3,info$nbrand+3)] = FALSE
Mcmc$fixedcoef = Mcmc$fixedcoef[!info$fixed]

# unfix everything except for brand pars 17 and 18 which were problematic
#Mcmc$fixedcoef[c(1:15,19,21)] = FALSE
#Mcmc$fixedcoef[17:18] = TRUE
#Mcmc$fixedcoef[c(19,21,22)] = FALSE

# degrees of freedom parameter on the prior variance
nu0 = sum(!Mcmc$fixedcoef) + 3
Mcmc$S0 = diag(max(1,sum(!Mcmc$fixedcoef)))*nu0*0.1
Mcmc$nu = nu0

# parameter names

pnkeep = c(paste("Product",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("Feat",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("Disp",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("PR",info$brsize[,1],"_",info$brsize[,2],sep=""))

pnames = c(pnkeep,"clb","cub","price","gamma","stockout","discount factor","??","FixedCost","StVisit","Omega1","Omega2")

pnames1 = pnames[c(1:info$nbrand,length(pnames)-10:0)]

if(sizeshifter) {
  pnames1 = c(pnames1,paste(info$packsize[2:nsize]*100/6.25,"OZ"))
}

Mcmc$pnames = pnames1

ind = which(pnames1[!info$fixed]=="price")

rr = which(pnames1[!info$fixed]=="Product10_5")

#Mcmc$fixedcoef[c(2,4,6:15,(ind-2):(ind+1),(ind+2):sum(!info$fixed))] = FALSE
Mcmc$fixedcoef[c(2,4,6:(rr-1),(rr+1):15,ind:(ind+2))] = FALSE

