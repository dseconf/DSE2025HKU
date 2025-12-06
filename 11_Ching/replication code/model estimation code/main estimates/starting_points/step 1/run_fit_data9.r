# This code will construct purchase probabilities/value functions for the initial plots in the paper

xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
xstart[which(pnames1[!info$fixed]=="discount factor")] = 0  #set the df starting point to 0.5
info$paramstart = NULL

source(paste(fpath,"test_functions.r",sep=""))

test = fitiv.hh.r(matrix(info$xfull,nrow=npbig,ncol=nhhs),hhmerge,info,1,1,nhhs)

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

#ll1 = ll(xstart,Mcmc,info,hhbig,hhmerge,hhmergebr,indflag=TRUE)

info1=info
info1$xfull[info$nbrand+6] = 0.001

info1$fixed[info$nbrand+6] = TRUE

x1=xstart
x1=x1[-which(pnames1[!info$fixed]=="discount factor")]

# try fitting everything, without dynamics or fwd-looking stuff

out.ml = optim(x1,ll,method="BFGS",
               control=list(maxit=10000,reltol=1e-5),Mcmc=Mcmc,info=info1,
               hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr)

save(out.ml,file="ml_fit.RData")
