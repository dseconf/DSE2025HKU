# estimate parameters after setup

#xstart = untform(matrix(info$xfull,nrow=npbig,ncol=nhhs),info)[,1]
#xstart[22] = -10
#info$paramstart = NULL

#out.ml = list(par=xstart)
#ml.vcov = diag(sum(!info$fixed))

#Mcmc$nprint = 1

cat("initial ML estimates\n")

tfmlpar = tform(matrix(out.ml$par,nrow=length(out.ml$par),ncol=nhhs),info)

print(cbind(tfmlpar[!info$fixed,1],sqrt(diag(ml.vcov))))

#info$fixed[length(info$fixed)-1] = TRUE

#stop()
# MCMC part

Mcmc$nsave = 10
Mcmc$ntilde = 10

# population-varying vs fixed coefficients

# this should be left to the setup file
#Mcmc$fixedcoef = rep(TRUE,npbig)
#if(init.crate & !info$fixed[info$nbrand+1]) {
#  tempfc = Mcmc$fixedcoef
#}
#Mcmc$fixedcoef = Mcmc$fixedcoef[!info$fixed]

# define proposals

if(init.crate & !info$fixed[info$nbrand+1]) {
  tempfixed = info$fixed
  tempfixed[info$nbrand+1] = TRUE
  tempfc = tempfc[!tempfixed]
  Mcmc$proposal = ml.vcov[tempfc,tempfc]
} else {
  Mcmc$proposal = ml.vcov[Mcmc$fixedcoef,Mcmc$fixedcoef]
}

#Mcmc$proposal = Mcmc$proposal[1:10,1:10]
#Mcmc$fixedcoef = Mcmc$fixedcoef[1:21]

nu0 = sum(!Mcmc$fixedcoef) + 1
Mcmc$S0 = diag(max(1,sum(!Mcmc$fixedcoef)))
Mcmc$nu = nu0

cat("Prior variance diagonal means:\n")
print(diag(Mcmc$S0/(Mcmc$nu-sum(!Mcmc$fixedcoef)-1)))

cat("Prior variance diagonal standard dev:\n")
print(sqrt(diag(Mcmc$S0)*diag(Mcmc$S0)*(2*(Mcmc$nu-sum(!Mcmc$fixedcoef)))/((Mcmc$nu-sum(!Mcmc$fixedcoef))*((Mcmc$nu-sum(!Mcmc$fixedcoef)-1)^2)*(Mcmc$nu-sum(!Mcmc$fixedcoef)-3))))

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

Mcmc$betabar[sum(!Mcmc$fixedcoef)-2,1] = log(0.4)
Mcmc$betaA[sum(!Mcmc$fixedcoef)-2,sum(!Mcmc$fixedcoef)-2] = 1

Mcmc$betabar[sum(!Mcmc$fixedcoef)-1,1] = log(0.85) - log(1-0.85)
Mcmc$betaA[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1] = 0.6^2

Mcmc$S0[sum(!Mcmc$fixedcoef)-2,sum(!Mcmc$fixedcoef)-2] = Mcmc$S0[sum(!Mcmc$fixedcoef)-2,sum(!Mcmc$fixedcoef)-2]/10
Mcmc$S0[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1] = Mcmc$S0[sum(!Mcmc$fixedcoef)-1,sum(!Mcmc$fixedcoef)-1]/5


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

#hhest = which(hhcrate > 0.25 & hhcrate < 2)

hhest = 1:nhhs

save(info,Mcmc,hhest,inputs,file=paste("mcmc_starting_params",spec.index,".RData",sep=""))

mcmc.out = MCMCLoops(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr,hhest)

save(mcmc.out,file=paste("est",spec.index,"_mcmcout.RData",sep=""))

#load("sim1_mcmcout.RData")

#stop()

# print output

burnin = as.integer(0.5*Mcmc$R)

nkeep = Mcmc$R/Mcmc$keep

betabig = array(mcmc.out$betadraw,c(nkeep,sum(!info$fixed),nhhs))
if(sum(!info$fixed)==1) {
  estmean1 = apply(array(betabig[burnin:Mcmc$R,,],c(Mcmc$R-burnin+1,1,nhhs)),c(1,2),mean)
  estsd1 = apply(array(betabig[burnin:Mcmc$R,,],c(Mcmc$R-burnin+1,1,nhhs)),c(1,2),sd)
} else {
  estmean1 = apply(betabig[burnin:Mcmc$R,,],c(1,2),mean)
  estsd1 = apply(betabig[burnin:Mcmc$R,,],c(1,2),sd)
}

estmean = apply(estmean1,2,mean)
estsd = apply(estsd1,2,mean)

estmeanci = apply(estmean1,2,quantile,probs=c(0.025,0.975))
estsdci = apply(estsd1,2,quantile,probs=c(0.025,0.975))

pnkeep = c(paste("Product",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("Feat",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("Disp",info$brsize[,1],"_",info$brsize[,2],sep=""),
paste("PR",info$brsize[,1],"_",info$brsize[,2],sep=""))

pnkeep = pnkeep[1:info$nbrand]

pnames = c(pnkeep,"clb","cub","price","gamma","stockout","discount factor","??","FixedCost","StVisit","Omega1","Omega2")

cat("Estimates (Individual Draws)\n")

outest = cbind(estmean,t(estmeanci),estsd,t(estsdci))
colnames(outest) = c("Average (Est)", "0.025%", "0.975%",
        "S.D. (Est)", "0.025%", "0.975%")
rownames(outest) = pnames[!info$fixed]

print(outest)

if(1==0) {

  # print moments of true draws alongsize moments of estimates

  tdtf = tform(truedraws,info)
  tdtf = tdtf[!info$fixed,]

  td.mean = apply(tdtf,1,mean)
  td.sd = apply(tdtf,1,sd)

  outmean = cbind(outest[,1:3],td.mean)
  colnames(outmean) = c(colnames(outest)[1:3],"True Mean")
  rownames(outmean) = pnames[!info$fixed]

  outsd = cbind(outest[,4:6],td.sd)
  colnames(outsd) = c(colnames(outest)[4:6],"True SD")
  rownames(outsd) = pnames[!info$fixed]

  print(outmean)
  print(outsd)

}

# heirarchical parameter moments

if(length(info$dvarnames)==1) {
    bbig = array(mcmc.out$bdraw,c(Mcmc$R,sum(!info$fixed)))

    estbmean1 = array(apply(bbig[burnin:Mcmc$R,],2,mean),c(sum(!info$fixed),1))
    estbsd1 = array(apply(bbig[burnin:Mcmc$R,],2,sd),c(sum(!info$fixed),1))

} else {
    bbig = array(mcmc.out$bdraw,c(Mcmc$R,sum(!info$fixed),ncol(info$dvars)))

    estbmean1 = apply(bbig[burnin:Mcmc$R,,],c(2,3),mean)
    estbsd1 = apply(bbig[burnin:Mcmc$R,,],c(2,3),sd)

}

estbout = matrix(0,nrow=2*nrow(estbmean1),ncol=ncol(estbmean1))
colnames(estbout) = c("Constant",colnames(info$dvars)[-1])
rownames(estbout) = rep("",nrow(estbout))

if(1==0) {
  truep2 = matrix(0,nrow=2*nrow(estbmean1),ncol=ncol(estbmean1))
  colnames(truep2) = paste("TR_",c("Constant",colnames(info$dvars)[-1]),sep="")
  for(i in 1:nrow(estbmean1)) {
    truep2[2*i-1,] = truep[i,]
  }
}

for(i in 1:nrow(estbmean1)) {
  estbout[2*i-1,] = estbmean1[i,]
  estbout[2*i,] = estbsd1[i,]
  rownames(estbout)[2*i-1] = pnames[!info$fixed][i]
  rownames(estbout)[2*i] = paste("SD: ",pnames[!info$fixed][i],sep="")
}

cat("Heirarchical mean coefficients\n")
print(estbout)

#cat("Estimates (Mean Coefficients) and True Parameters\n")
#print(cbind(estbout,truep2))

#cat("standard deviations\n")

#if(sum(!Mcmc$fixedcoef)==1) {
#  estsdout1 = mean(estsd1[,!Mcmc$fixedcoef])
#} else {
#  estsdout1 = apply(estsd1[,!Mcmc$fixedcoef],1,mean)
#}
#print(cbind(estsdout1,truesd[!Mcmc$fixedcoef]))

#colnames(estbmean1) = c("Constant",colnames(info$dvars)[-1])
#colnames(estbsd1) = c("Constant",colnames(info$dvars)[-1])

if(sum(!Mcmc$fixedcoef) > 0) {

  Wbig = array(mcmc.out$Wdraw,c(Mcmc$R,sum(!info$fixed),sum(!info$fixed)))

  Wsmall = Wbig[,!Mcmc$fixedcoef,!Mcmc$fixedcoef]

  Wdiag = matrix(0,nrow=Mcmc$R,sum(!Mcmc$fixedcoef))

  for(i in 1:sum(!Mcmc$fixedcoef)) {

    Wdiag[,i] = Wsmall[,i,i]

  }

  estW = apply(sqrt(Wdiag[burnin:Mcmc$R,]),2,mean)
  ciW = apply(sqrt(Wdiag[burnin:Mcmc$R,]),2,quantile,probs=c(0.025,0.975))

  cat("est standard deviations\n")
  print(cbind(estW,t(ciW)))

}

# marginal log-likelihoods

cat("Marginal log-likelihood (overall): ",mean(mcmc.out$llqsave[burnin:Mcmc$R]+mcmc.out$llbrsave[burnin:Mcmc$R]),"\n")

cat("Marginal log-likelihood (quantity): ",mean(mcmc.out$llqsave[burnin:Mcmc$R]),"\n")

cat("Marginal log-likelihood (brand): ",mean(mcmc.out$llbrsave[burnin:Mcmc$R]),"\n")


Rsave = 10000
Rstart = 5000
ndtrue = 5
nhhbig = 312

bdrawsave = array(mcmc.out$bdraw,c(Rsave,npsmall,ndtrue))

betabig = array(mcmc.out$betadraw,c(Rsave,npsmall,nhhbig))

betamean = apply(betabig[Rstart:Rsave,,],2,FUN=mean)

bdrawm = apply(bdrawsave[Rstart:Rsave,,],c(2,3),FUN=mean)

save(betamean, bdrawm, file = "bdrawstart.RData")
