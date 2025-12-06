# do initial estimation: maximum likelihood to get proposals for fixed parameters and starting points

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
xstart[which(pnames1[!info$fixed]=="discount factor")] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

load("ml_fit.RData")

xout2 = rep(0,sum(!info$fixed))
ind = which(pnames1[!info$fixed]=="discount factor")
xout2[1:(ind-1)] = out.ml$par[1:(ind-1)]
xout2[ind] = 0
xout2[(ind+1):sum(!info$fixed)] = out.ml$par[ind:length(out.ml$par)]

xx = info$xfull
#if(sizeshifter) {
#  xx[!info$fixed] = c(xout2,rep(0,info$nsize-1))
#} else {
xx[!info$fixed] = xout2
#}
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

    #out.ml = optim(xstart,ll,control=list(maxit=5000,reltol=1e-4),Mcmc=Mcmc,info=info,
    #               hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)

    #load("out.br.RData")
    #xout2[19] = 0
    #xout2[20] = log(0.4)
    #xout2[21] = log(0.8)-log(0.2)
    #xout2[22] = log(3.5)

    if(sizeshifter) {
      out.ml = list(par=c(xout2,rep(0,info$nsize-1)))
    } else {
      out.ml = list(par=xout2)
    }

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

    #out.ml = optim(xstart,ll,control=list(maxit=5000,reltol=1e-4),Mcmc=Mcmc.ll,info=info,
    #               hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)

    #load("ml_bfgsA2_fix1.RData")

    #out.ml = out.ml2
    #xin = out.ml$par

    #load("out.br.RData")
    #xout2[9] = -4
    #xout2[19] = 0
    #xout2[20] = log(0.4)
    #xout2[21] = log(0.8)-log(0.2)
    #xout2[22] = log(3.5)

    #if(sizeshifter) {
    #  out.ml = list(par=c(xout2,rep(0,info$nsize-1)))
    #} else {
    out.ml = list(par=xout2)
    #}

    #xx = matrix(info$xfull, nrow=npbig, ncol=nhhs)
    #xx[!info$fixed,] = xin
    #out.ml = list(par = untform(xx,info)[,1],value = 28986.62)

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

#stop()


