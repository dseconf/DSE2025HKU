# run counterfactual analysis

# create new counterfactual price data

# create cfdata and cfinit

sourceCpp(paste(fpath,"ijc_Rcpp_par.cpp",sep=""))

prcols = which(substr(colnames(hhmerge),1,7) == "Product")

cfdata = as.matrix(hhmerge[,prcols],nrow=nrow(hhmerge),ncol=length(prcols))

prcolsi = which(substr(colnames(hhinitmerge),1,7) == "Product")

cfinit = as.matrix(hhinitmerge[,prcolsi],nrow=nrow(hhinitmerge),ncol=length(prcolsi))

info.old = info

# re-index prices

bylist = NULL

i=0
for(v in colnames(hhmerge)[1:info$nbrand + info$prindex - 1]) {
  i=i+1
  bylist[[i]] = cfdata[,i]
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
  bylist1[[i+1]] = cfdata[,i]
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

# redefine dropmat

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
xstart[22-as.numeric(data.crate)] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

load("bdrawstart.RData")

ind = which(pnames1[!info$fixed]=="price")

betamean = betamean[-(ind+2)]
bdrawm = bdrawm[-(ind+2),]

xx = info$xfull
xx[!info$fixed] = betamean

hhmerge1 = hhmerge
hhmerge1[,prcols] = cfdata

test = fitiv.hh.r(matrix(xx,nrow=npbig,ncol=nhhs),hhmerge1,info,1,1,nhhs)

rm(hhmerge1)

info$dropmat = test[[2]]

info$ivdrop = test[[3]] > 1e12

cat("Number of households with static IV process:\n")
print(table(apply(info$dropmat,1,sum)))

cat("Number of households with simplified IV process:\n")
print(table(info$ivdrop))

cat("Summary of max condition numbers:\n")
print(summary(test[[3]]))

# run code to estimate counterfactuals

load(paste("est",spec.index,"_mcmcout.RData",sep=""))

nkeep = Mcmc$R/Mcmc$keep

betadraw = array(mcmc.out$betadraw,c(nkeep,npsmall,nhhs))

Mcmc$nsave = 10
Mcmc$ntilde = 10

# I'm not sure if we need to bother with this stuff, but I'll copy it anyways

load(paste("estml",spec.index,".RData",sep=""))

# define proposals

if(init.crate & !info$fixed[info$nbrand+1]) {
  tempfixed = info$fixed
  tempfixed[info$nbrand+1] = TRUE
  tempfc = tempfc[!tempfixed]
  Mcmc$proposal = ml.vcov[tempfc,tempfc]
} else {
  Mcmc$proposal = ml.vcov[Mcmc$fixedcoef,Mcmc$fixedcoef]
}

nu0 = sum(!Mcmc$fixedcoef) + 1
Mcmc$S0 = diag(max(1,sum(!Mcmc$fixedcoef)))
Mcmc$nu = nu0

info$dvarnames = info$dvarnames[1:ndvars]

if(ndvars == 1) {
  info$dvars = matrix(info$dvars[,1:ndvars],nrow = info$nhhs,ncol=1)
} else {
  info$dvars = info$dvars[,1:ndvars]
}

Mcmc$betabar = matrix(0,nrow=sum(!Mcmc$fixedcoef),ncol=ndvars)
Mcmc$betaA = 1000*diag(sum(!Mcmc$fixedcoef)*ndvars)
if(ndvars > 1) {
  Mcmc$bA = 1000*diag(ndvars)
}

Mcmc$bbarAfixed = 1000*diag(sum(Mcmc$fixedcoef))

# set prior for stockout and df

Mcmc$betabar[sum(!Mcmc$fixedcoef)-1,1] = log(0.4)
Mcmc$betaA[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1] = 1

Mcmc$S0[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1] = Mcmc$S0[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1]/10

pstart1 = info$xfull
if(init.crate) {
  pstart1[!tempfixed] = out.ml$par
} else {
  pstart1[!info$fixed] = out.ml$par
}

inputs = list(pstart = matrix(pstart1,nrow=length(pstart1),ncol=nhhs))

if(init.crate) {
  if(!is.finite(crub)) {
    inputs$pstart[info$nbrand+1,] = log(hhcrate)
    inputs$pstart[info$nbrand+2,] = log(hhcrate)
  } else {
    inputs$pstart[info$nbrand+1,] = log(hhcrate/crub) - log(1.0-hhcrate/crub)
    inputs$pstart[info$nbrand+2,] = log(hhcrate/crub) - log(1.0-hhcrate/crub)
  }
}

cf.out = MCMCcounterfactual(betadraw,Mcmc,info,hhbig,hhmerge,hhmergebr,
                            cfdata,cfinit,as.integer(1:nhhs),stvisitinit)

save(cf.out,file=paste("cfout",cf.index,".RData",sep=""))


