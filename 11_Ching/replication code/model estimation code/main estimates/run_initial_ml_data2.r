# do initial estimation: maximum likelihood to get proposals for fixed parameters and starting points

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
xstart[which(pnames1[!info$fixed]=="discount factor")] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

# this was in the earlier code and loaded in the original output from run 45.  The file is huge
# so I replaced it with bdrawstart.RData
#load("est1_mcmcout.RData")

#Rsave = 10000
#Rstart = 5000
#ndtrue = 5

#bdrawsave = array(mcmc.out$bdraw,c(Rsave,npsmall,ndtrue))

#betabig = array(mcmc.out$betadraw,c(Rsave,npsmall,nhhs))

#betamean = apply(betabig[Rstart:Rsave,,],2,FUN=mean)

#bdrawm = apply(bdrawsave[Rstart:Rsave,,],c(2,3),FUN=mean)

#rm(mcmc.out,bdrawsave,betabig)

load("bdrawstart.RData")

xx = info$xfull
xx[!info$fixed] = betamean

test = fitiv.hh.r(matrix(xx,nrow=npbig,ncol=nhhs),hhmerge,info,1,1,nhhs)

info$dropmat = test[[2]]

info$ivdrop = test[[3]] > 1e12

cat("Number of households with static IV process:\n")
print(table(apply(info$dropmat,1,sum)))

cat("Number of households with simplified IV process:\n")
print(table(info$ivdrop))

cat("Summary of max condition numbers:\n")
print(summary(test[[3]]))

if(loadml) {

  load(paste("estml",spec.index,".RData",sep=""))

} else {

  if(info$h.model) {

    out.ml = optim(xstart,ll,control=list(maxit=5000,reltol=1e-4),Mcmc=Mcmc,info=info,
                   hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)

    if(slowhess) {
      temphh = ll.hess(out.ml$par,Mcmc,info,hhbig,hhmerge,hhmergebr)
      save(temphh,file=paste("numderivhessian",spec.index,".RData",sep=""))
      ml.vcov = chol2inv(chol(temphh))
    } else {
      par.vcov = llcov(out.ml$par,Mcmc,info,hhbig,hhmerge,hhmergebr,gradtol=1e-3)
      ml.vcov = par.vcov[[1]]
    }

  } else {

    Mcmc.ll = Mcmc
    Mcmc.ll$fixedcoef = rep(TRUE,length(Mcmc$fixedcoef))

    info$h.model = TRUE

    if(data.crate) {
      info$crate[1:nhhs] = mean(hhcrate)
    }

    out.ml = list(par = bdrawm[,1])

    rm(mcmc.out)

    if(slowhess) {
      temphh = ll.hess(out.ml$par,Mcmc.ll,info,hhbig,hhmerge,hhmergebr)
      save(temphh,file=paste("numderivhessian",spec.index,".RData",sep=""))
      ml.vcov = chol2inv(chol(temphh))
    } else {
      par.vcov = llcov(out.ml$par,Mcmc.ll,info,hhbig,hhmerge,hhmergebr,gradtol=1e-3)
      ml.vcov = par.vcov[[1]]
    }

    if(spec.index == 12) {
      ml.vcov[23,23] = 1  #I'm doing this because it gives me zero in this specific instance
    }

    if(data.crate) {
      info$crate[1:nhhs] = hhcrate
    }

    info$h.model = FALSE

    rm(Mcmc.ll)

  }

  save(out.ml,ml.vcov,file=paste("estml",spec.index,".RData",sep=""))

}



