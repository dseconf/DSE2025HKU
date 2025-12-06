# do initial estimation: maximum likelihood to get proposals for fixed parameters and starting points

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
xstart[which(pnames1[!info$fixed]=="discount factor")] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

#load("est1_mcmcout.RData")

Rsave = 10000
Rstart = 5000
ndtrue = 5
nhhbig = 312

#bdrawsave = array(mcmc.out$bdraw,c(Rsave,npsmall,ndtrue))

#betabig = array(mcmc.out$betadraw,c(Rsave,npsmall,nhhbig))

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

if(h.model) {
  info$crate[1:nhhs] = mean(hhcrate)
}

# check log-likelihood
#test1 = ll(xstart,Mcmc,info,hhbig,hhmerge,hhmergebr,indflag=TRUE)


# compare fits

ibds = seq(0,50,0.5)

purch1 = emp.probs(hhbig,hhmerge,hhinitbig,hhcrate,brsize,vols,first)

pprobs1 = purch1$pprobs

ind=which(pnames1[!info$fixed]=="discount factor")

# try this - max over small set of parms instead of all

xxdf = bdrawm[ind,1]
xf = rep(TRUE,length(bdrawm[,1]))
xf[21] = FALSE

info1 = info
info1$xfull = tform(matrix(bdrawm[,1],nrow=npsmall,ncol=nhhs),info)[,1]
info1$fixed = rep(TRUE,npbig)
info1$fixed[info$nbrand+6] = FALSE

Mcmc1 = Mcmc
Mcmc1$fixedcoef = rep(TRUE,1)

out.mldf = optimize(ll,c(log(0.01)-log(0.99),log(0.995)-log(1-0.995)),Mcmc=Mcmc1,info=info1,
               hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)
save(out.mldf,file=paste("out_mldf_A_fix1_",fit.index,"_",clust.index,".RData",sep=""))

info1 = info
info1$xfull = tform(matrix(bdrawm[,1],nrow=npsmall,ncol=nhhs),info)[,1]
info1$fixed = rep(TRUE,npbig)
info1$fixed[c(info$nbrand+5,info$nbrand+6,info$nbrand+8)] = FALSE

Mcmc1 = Mcmc
Mcmc1$fixedcoef = rep(TRUE,3)

xstart = bdrawm[(ind-1):(ind+1),1]
xstart[2] = out.mldf$minimum

out.ml = optim(xstart,ll,method="BFGS",
               control=list(maxit=10000),Mcmc=Mcmc1,info=info1,
               hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)
save(out.ml,file=paste("out_ml_A_fix1_",fit.index,"_",clust.index,".RData",sep=""))
