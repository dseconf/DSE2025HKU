# do initial estimation: maximum likelihood to get proposals for fixed parameters and starting points

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
#xstart[22-as.numeric(data.crate)] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

load("bdrawstart.RData")

ind = which(pnames1[!info$fixed]=="price")

betamean = betamean[-(ind+2)]
bdrawm = bdrawm[-(ind+2),]

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



