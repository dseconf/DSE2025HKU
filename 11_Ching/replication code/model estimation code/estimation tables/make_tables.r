# create tables of estimates for the paper

rm(list=ls())

# do we do multicore or not
multicore = TRUE

library(randtoolbox)
library(bayesm)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)
library(Matrix)
library(MASS)
if(multicore) {

  library(RcppParallel)

  setThreadOptions(stackSize = 1024*1024*1024)
}

# screen width

options("width"=200)

desktop = TRUE

indests = FALSE  # show estimates based on averages of individual level draws

# paths to data, code files, and output

datapath = ""
fpath = ""
tablepath = ""
figpath = ""


source(paste(fpath,"ijc_functions_new.r",sep=""))

# specification index

spec.index = 2

source(paste("estimation_setup",spec.index,".r",sep=""))

load(paste(fpath,"mcmc_starting_params",spec.index,".RData",sep=""))

# load in the data with the brand information in it

yind = 0
for(yearnum in yearlist) {
  yind = yind + 1

  load(paste(datapath,"brands_",yearnum,"_",product,".RData",sep=""))

  brandinfo = brandinfo[,-7]

  load(paste(datapath,"hhpanel_",product,"_",yearnum,".RData",sep=""))

  hhpanel = hhpanel[,-15]

  if(yind == 1) {
    binfo.big = brandinfo
    hhpanel.big = hhpanel
  } else {
    if(yind == 7) {binfo.big = rbind(binfo.big,brandinfo)}
    hhpanel.big = rbind(hhpanel.big,hhpanel)
  }

}

# merge together brand info with purchases to back out brand names

hhm1 = hhmergebr[,c("PANID","WEEK","brindex")]
hhp1 = hhpanel.big[,c("PANID","WEEK","L5","L9","VOL_EQ")]
hhm1$VOL_EQ = info$packsize[brsize[hhm1$brindex,2]]

hhp1 = merge(hhp1,hhm1,all.x=TRUE)
hhp1 = hhp1[!is.na(hhp1$brindex),]

brnames = rep("",info$nbrand)
brvols = rep("",info$nbrand)

for(i in 1:info$nbrand) {

  temp = table(hhp1$L5[hhp1$brindex==i])

  if(length(temp) > 0) {

    brnames[i] = names(temp)[which.max(temp)]

    brvols[i] = paste(info$packsize[brsize[i,2]]*16," OZ",sep="")
  }

}

fixname <- function(x) {paste(toupper(substr(x,1,1)),tolower(substr(x,2,nchar(x))),sep="")}

fixsentence <- function(x) {

  y = strsplit(x,split=" ")[[1]]
  res = ""
  for(i in 1:length(y)) {
    if(i < length(y)) {
      res = paste(res,fixname(y[i])," ",sep="")
    } else {
      res = paste(res,fixname(y[i]),sep="")
    }
  }
  res

}

for(i in 1:length(brnames)) {
  brnames[i] = fixsentence(brnames[i])
}

save(brnames,file="brnames.RData")

# add in sizes for freed up brands

brnames[!info$sizebrand] = paste(brnames[!info$sizebrand],", ",info$packsize[info$brsize[!info$sizebrand,2]]*100/6.25," OZ",sep="")

pnkeep = c(brnames,
paste("Feat",brnames,"_",brvols,sep=""),
paste("Disp",brnames,"_",brvols,sep=""),
paste("PR",brnames,"_",brvols,sep=""))

pnkeep = pnkeep[1:info$nbrand]

if(invmodel==1) {
  pnames = c(pnkeep,"Consumption Rate","Consumption Rate (Upper Bound)","Price Coefficient",
             "gamma","Stockout Cost","Discount Factor","??","Fixed Cost of Purchase","StVisit","Inventory Bound","Omega2")
} else {
  pnames = c(pnkeep,"Consumption Rate","Consumption Rate (Upper Bound)","Price Coefficient",
             "gamma","Stockout Cost","Discount Factor","??","Fixed Cost of Purchase","StVisit","SC Linear","SC Quadratic","Inadmissable State Penalty")
}

# load in data and make tables of results

load(paste(fpath,"est",spec.index,"_mcmcout.RData",sep=""))

if(info$sizeshifter) {
  pnames = c(pnames,paste(info$packsize[2:nsize]*100/6.25,"OZ"))
}

nhhs = info$nhhs

burnin = as.integer(0.5*Mcmc$R)

nkeep = Mcmc$R/Mcmc$keep

npsmall = sum(!info$fixed)

betabig = array(mcmc.out$betadraw,c(nkeep,npsmall,nhhs))

# need to divide price coeff by info$bigJ
pindx = which(Mcmc$pnames[!info$fixed]=="price")

betabig[,pindx,] = betabig[,pindx,]/info$bigJ

# rescale storage costs to be in 100 oz bottle equivalents

if(info$invmodel==2) {
  scindx = which(Mcmc$pnames[!info$fixed] == "Omega1")
  betabig[,scindx,] = betabig[,scindx,]*6.25
  betabig[,scindx+1,] = betabig[,scindx+1,]*6.25*6.25
}


qvalues = c(0.33,0.5,0.66)  # values of quantiles for population distribution of params
nq = length(qvalues)
estq1 = array(0,c(Mcmc$R-burnin+1,npsmall,nq))

if(sum(!info$fixed)==1) {
  estmean1 = apply(array(betabig[burnin:Mcmc$R,,],c(Mcmc$R-burnin+1,1,nhhs)),c(1,2),mean)
  estsd1 = apply(array(betabig[burnin:Mcmc$R,,],c(Mcmc$R-burnin+1,1,nhhs)),c(1,2),sd)
  for(i in 1:nq) {
    estq1[,,i] = apply(array(betabig[burnin:Mcmc$R,,],c(Mcmc$R-burnin+1,1,nhhs)),c(1,2),quantile,probs=qvalues[i])
  }
} else {
  estmean1 = apply(betabig[burnin:Mcmc$R,,],c(1,2),mean)
  estsd1 = apply(betabig[burnin:Mcmc$R,,],c(1,2),sd)
  for(i in 1:nq) {
    estq1[,,i] = apply(betabig[burnin:Mcmc$R,,],c(1,2),quantile,probs=qvalues[i])
  }
  estmean2 = apply(betabig[burnin:Mcmc$R,,],c(2,3),mean)

  if(spec.index != 62 & spec.index != 67) {
    estmode2 = matrix(0,nrow=sum(!info$fixed),ncol=nhhs)
    betabigu = betabig
    for(i in 1:sum(!info$fixed)) {
      if(abs(info$tform[!info$fixed][i]) == 1) {
        betabigu[,i,] = log(info$tform[!info$fixed][i]*betabigu[,i,])
      } else if (abs(info$tform[!info$fixed][i]) == 2) {
        betabigu[,i,] = log(betabigu[,i,])-
          log(1-betabigu[,i,])
      }
      for(n in 1:nhhs) {
        d = density(betabigu[burnin:Mcmc$R,i,n])
        estmode2[i,n] = d$x[which.max(d$y)]
      }
    }
    estmode2 = tform(estmode2,info)[!info$fixed,]
  }
}

estmean = apply(estmean1,2,mean)
estsd = apply(estsd1,2,mean)
estq = array(0,c(npsmall,nq))
for(i in 1:nq) {
  estq[,i] = apply(estq1[,,i],2,mean)
}

estmeanci = apply(estmean1,2,quantile,probs=c(0.025,0.975))
estsdci = apply(estsd1,2,quantile,probs=c(0.025,0.975))
estqci = array(0,c(npsmall,2,nq))

for(i in 1:nq) {
  estqci[,,i] = t(apply(estq1[,,i],2,quantile,probs=c(0.025,0.975)))
}

# confidence intervals for each individual

estci.ind.l = apply(betabig[burnin:Mcmc$R,,],c(2,3),quantile,probs=c(0.025))
estci.ind.u = apply(betabig[burnin:Mcmc$R,,],c(2,3),quantile,probs=c(0.975))

# save estimates for DIC calculation

save(estmean2,file=paste("param_means",spec.index,".RData",sep=""))

sink(paste("bayesian_estimates",spec.index,".txt",sep=""))
cat("Estimates (Population Mean and SD of Individual Draws)\n")

outest = cbind(estmean,t(estmeanci),estsd,t(estsdci))
colnames(outest) = c("Average (Est)", "0.025%", "0.975%",
        "S.D. (Est)", "0.025%", "0.975%")
rownames(outest) = pnames[!info$fixed]
print(outest)
cat("\n")

outq = NULL
for(i in 1:nq) {
  outq = cbind(outq,estq[,i],estqci[,,i])
}

colnames(outq) = paste("X",1:(3*nq),sep="")

for(i in 1:nq) {
  colnames(outq)[1:3 + 3*(i-1)] = c(paste(qvalues[i],sep=""),"0.025%", "0.975%")
}
rownames(outq) = pnames[!info$fixed]

cat("Estimates (Population Percentiles of Individual Draws)\n")
print(outq)
cat("\n")

bbig = array(mcmc.out$bdraw,c(Mcmc$R,sum(!info$fixed),ncol(info$dvars)))

if(ndvars == 1) {
  estbmean1 = matrix(apply(bbig[burnin:Mcmc$R,,],2,mean),nrow=sum(!info$fixed),ncol=1)
  estbsd1 = matrix(apply(bbig[burnin:Mcmc$R,,],2,sd),nrow=sum(!info$fixed),ncol=1)
} else {
  estbmean1 = apply(bbig[burnin:Mcmc$R,,],c(2,3),mean)
  estbsd1 = apply(bbig[burnin:Mcmc$R,,],c(2,3),sd)
}
estbout = matrix(0,nrow=2*nrow(estbmean1),ncol=ncol(estbmean1))
colnames(estbout) = c("Constant",colnames(info$dvars)[-1])
rownames(estbout) = rep("",nrow(estbout))

for(i in 1:nrow(estbmean1)) {
  estbout[2*i-1,] = estbmean1[i,]
  estbout[2*i,] = estbsd1[i,]
  rownames(estbout)[2*i-1] = pnames[!info$fixed][i]
  rownames(estbout)[2*i] = paste("SD: ",pnames[!info$fixed][i],sep="")
}

cat("Estimates (Underlying Demographic Coefficients)\n")
print(estbout)
sink()

# latex tables

if(invmodel==1) {
  sigfig = 2
} else {
  sigfig = 4
}
if(spec.index==62 | spec.index == 69) {
  sigfig=4
}

n1 = gsub("&","\\\\&",pnames)[!info$fixed]
allesttab = rep("",2*sum(!info$fixed))
for(i in 1:(length(allesttab)/2)) {
  allesttab[2*i-1] = n1[i]
}

allesttab.sd = rep("",2*sum(!info$fixed))
for(i in 1:(length(allesttab.sd)/2)) {
  allesttab.sd[2*i-1] = n1[i]
}


nqmid = 2

for(i in 1:(length(allesttab)/2)) {

                                        # do the sd version
  allesttab.sd[2*i-1] = paste(allesttab.sd[2*i-1]," & ",round(estmean[i],sigfig),sep="")
  allesttab.sd[2*i] = paste(allesttab.sd[2*i]," & [",round(estmeanci[1,i],sigfig),", ",round(estmeanci[2,i],sigfig),"]",sep="")
  if(Mcmc$fixedcoef[i]) {
    allesttab.sd[2*i-1] = paste(allesttab.sd[2*i-1]," & -",sep="")
    allesttab.sd[2*i] = paste(allesttab.sd[2*i]," & ",sep="")
  } else {
    allesttab.sd[2*i-1] = paste(allesttab.sd[2*i-1]," & ",round(estsd[i],sigfig),sep="")
    allesttab.sd[2*i] = paste(allesttab.sd[2*i]," & [",round(estsdci[1,i],sigfig),", ",round(estsdci[2,i],sigfig),"]",sep="")
  }

  if(Mcmc$fixedcoef[i]) {
    for(k in 1:nqmid) {
      allesttab[2*i-1] = paste(allesttab[2*i-1]," & -",sep="")
      allesttab[2*i] = paste(allesttab[2*i]," & ",sep="")
    }
  } else {
    for(k in 1:nqmid) {
      allesttab[2*i-1] = paste(allesttab[2*i-1]," & ",round(estq[i,k],sigfig),sep="")
      allesttab[2*i] = paste(allesttab[2*i]," & [",round(estqci[i,1,k],sigfig),", ",round(estqci[i,2,k],sigfig),"]",sep="")
    }
  }
  allesttab[2*i-1] = paste(allesttab[2*i-1]," & ",round(estmean[i],sigfig),sep="")
  allesttab[2*i] = paste(allesttab[2*i]," & [",round(estmeanci[1,i],sigfig),", ",round(estmeanci[2,i],sigfig),"]",sep="")
  if(Mcmc$fixedcoef[i]) {
    for(k in (nqmid+1):nq) {
      allesttab[2*i-1] = paste(allesttab[2*i-1]," & -",sep="")
      allesttab[2*i] = paste(allesttab[2*i]," & ",sep="")
    }
  } else {
    for(k in (nqmid+1):nq) {
      allesttab[2*i-1] = paste(allesttab[2*i-1]," & ",round(estq[i,k],sigfig),sep="")
      allesttab[2*i] = paste(allesttab[2*i]," & [",round(estqci[i,1,k],sigfig),", ",round(estqci[i,2,k],sigfig),"]",sep="")
    }
  }
  allesttab[2*i-1] = paste(allesttab[2*i-1]," \\\\",sep="")
  allesttab[2*i] = paste(allesttab[2*i]," \\\\",sep="")
}


lltotbar = mcmc.out$llbrsave + mcmc.out$llqsave

margll = mean(lltotbar[burnin:Mcmc$R])

allesttab = c(allesttab,"\\hline",paste("Log-likelihood & ",round(margll,sigfig)," \\\\",sep=""))

if(info$sizeshifter) {
  nbrcoef = sum(!info$fixed[1:info$nbrand])
  n = length(allesttab)

  brandtab = allesttab[c(1:(2*nbrcoef),(n-2*(info$nsize-1)-1):(n-2))]

  dyntab = allesttab[c((2*nbrcoef+1):(n-2*(info$nsize-1)-2),(n-1):n)]

  brandtab.sd = allesttab.sd[c(1:(2*nbrcoef),(n-2*(info$nsize-1)-1):(n-2))]

} else {
  nbrcoef = sum(!info$fixed[1:info$nbrand])

  brandtab = allesttab[1:(2*nbrcoef)]

  dyntab = allesttab[(2*nbrcoef+1):length(allesttab)]

  brandtab.sd = allesttab.sd[1:(2*nbrcoef)]

}

brandtab.sd2 = brandtab.sd[1:26]
for(i in 1:24) {
  brandtab.sd2[i] = paste(brandtab.sd2[i]," & ",brandtab.sd[26+i]," \\\\",sep="")
}
for(i in 25:26) {
  brandtab.sd2[i] = paste(brandtab.sd2[i]," & & & \\\\",sep="")
}

write.table(brandtab.sd2,file=paste(tablepath,"brandests",spec.index,"_sd.tex",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(brandtab,file=paste(tablepath,"brandests",spec.index,".tex",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

write.table(dyntab,file=paste(tablepath,"dynests",spec.index,".tex",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

if(!info$fixed[info$nbrand+6]) {

  # discount factor density plot

  dfind = sum(!info$fixed[1:(info$nbrand+6)])

  df.est = (estmean2[dfind,])

  pdf(paste(figpath,"dfindests",spec.index,".pdf",sep=""))
  plot(density(df.est,from=0,to=1),main="",xlab = "Discount Factor Estimate")
  dev.off()

}

