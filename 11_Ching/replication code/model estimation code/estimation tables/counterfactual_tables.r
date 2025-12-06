# make counterfactual tables

rm(list=ls())


options("width"=200)
datapath = ""
fpath = ""
cfpath = ""       #path to counterfactual output

# do we do multicore or not
multicore = TRUE

library(randtoolbox)
library(bayesm)
library(Rcpp)
library(RcppArmadillo)
library(Matrix)
if(multicore) {

  library(RcppParallel)

  setThreadOptions(stackSize = 1024*1024*1024)
}

spec.index=2

source(paste(fpath,"ijc_functions_new.r",sep=""))

source("estimation_setup2.r")

burnin = 10000
Nsave = 20000
nbrand = info$nbrand

includedprod = 1:nbrand

avgq = vector("list",length=8)
avgr = vector("list",length=8)

ciq = vector("list",length=8)
cir = vector("list",length=8)

cf = vector("list",length=8)

# these counterfactuals are
# 1 - data prices, est discount factor
# 2 - higher freq, est discount factor
# 3 - higher depth, est discount factor
# 4 - data prices, fixed discount factor
# 5 - higher freq, fixed discount factor
# 6 - higher depth, fixed discount factor
# 7: 2, with p expectations from data
# 8: 3, with p expectations from data
# 9: 5, with p expectations from data
# 10: 6, with p expectations from data

cflist = 1:10

for(i in 1:length(cflist)) {

  load(paste(cfpath,"cfout",cflist[i],".RData",sep=""))
  avgq[[i]] = apply(cf.out$quantities[burnin:Nsave,includedprod],2,mean)
  avgr[[i]] = apply(cf.out$revenues[burnin:Nsave,includedprod],2,mean)
  ciq[[i]] = apply(cf.out$quantities[burnin:Nsave,includedprod],2,quantile,probs=c(0.025,0.975))
  cir[[i]] = apply(cf.out$revenues[burnin:Nsave,includedprod],2,quantile,probs=c(0.025,0.975))

  cf[[i]] = cf.out

}

qchange = vector("list",length=6)

revchange = vector("list",length=6)

qchange[[1]] = (cf[[2]]$quantities - cf[[1]]$quantities)
qchange[[2]] = (cf[[3]]$quantities - cf[[1]]$quantities)
qchange[[3]] = (cf[[5]]$quantities - cf[[4]]$quantities)
qchange[[4]] = (cf[[6]]$quantities - cf[[4]]$quantities)
qchange[[5]] = (cf[[7]]$quantities - cf[[1]]$quantities)
qchange[[6]] = (cf[[8]]$quantities - cf[[1]]$quantities)
qchange[[7]] = (cf[[9]]$quantities - cf[[4]]$quantities)
qchange[[8]] = (cf[[10]]$quantities - cf[[4]]$quantities)

revchange[[1]] = (cf[[2]]$revenues - cf[[1]]$revenues)
revchange[[2]] = (cf[[3]]$revenues - cf[[1]]$revenues)
revchange[[3]] = (cf[[5]]$revenues - cf[[4]]$revenues)
revchange[[4]] = (cf[[6]]$revenues - cf[[4]]$revenues)
revchange[[5]] = (cf[[7]]$revenues - cf[[1]]$revenues)
revchange[[6]] = (cf[[8]]$revenues - cf[[1]]$revenues)
revchange[[7]] = (cf[[9]]$revenues - cf[[4]]$revenues)
revchange[[8]] = (cf[[10]]$revenues - cf[[4]]$revenues)

avgcq = vector("list",length=8)
avgcr = vector("list",length=8)

cicq = vector("list",length=8)
cicr = vector("list",length=8)

rchangetot = vector("list",length=8)
rchangeci = vector("list",length=8)
rchangesd = vector("list",length=8)

qchangetot = vector("list",length=8)
qchangeci = vector("list",length=8)

for(i in 1:8) {

  avgcq[[i]] = apply(qchange[[i]][burnin:Nsave,includedprod],2,mean)
  avgcr[[i]] = apply(revchange[[i]][burnin:Nsave,includedprod],2,mean)
  cicq[[i]] = apply(qchange[[i]][burnin:Nsave,includedprod],2,quantile,probs=c(0.025,0.975))
  cicr[[i]] = apply(revchange[[i]][burnin:Nsave,includedprod],2,quantile,probs=c(0.025,0.975))

  qchangetot[[i]] = mean(apply(qchange[[i]][burnin:Nsave,includedprod],1,sum))
  qchangeci[[i]] = quantile(apply(qchange[[i]][burnin:Nsave,includedprod],1,sum),probs=c(0.025,0.975))

  rchangetot[[i]] = mean(apply(revchange[[i]][burnin:Nsave,includedprod],1,sum))
  rchangeci[[i]] = quantile(apply(revchange[[i]][burnin:Nsave,includedprod],1,sum),probs=c(0.025,0.975))
  rchangesd[[i]] = sd(apply(revchange[[i]][burnin:Nsave,includedprod],1,sum))

}

# show category changes, for the main text

cat("Promotional depth increase (long-term):\n")
cis.tot.q = matrix(0,nrow=4,ncol=2)
cis.tot.r = matrix(0,nrow=4,ncol=2)
j=0
for(i in c(2,4,6,8)) {
  j=j+1
  temp = apply(qchange[[i]][burnin:Nsave,includedprod],1,sum)
  cis.tot.q[j,] = quantile(temp,probs=c(0.025,0.975))
  temp = apply(revchange[[i]][burnin:Nsave,includedprod],1,sum)
  cis.tot.r[j,] = quantile(temp,probs=c(0.025,0.975))
}

cat("Baseline:\n")
print(c(sum(avgcq[[2]]),cis.tot.q[1,]))
print(c(sum(avgcr[[2]]),cis.tot.r[1,]))
cat("Rational Expectations:\n")
print(c(sum(avgcq[[4]]),cis.tot.q[2,]))
print(c(sum(avgcr[[4]]),cis.tot.r[2,]))

cat("Promotional depth increase (short-term):\n")

cat("Baseline:\n")
print(c(sum(avgcq[[6]]),cis.tot.q[3,]))
print(c(sum(avgcr[[6]]),cis.tot.r[3,]))
cat("Rational Expectations:\n")
print(c(sum(avgcq[[8]]),cis.tot.q[4,]))
print(c(sum(avgcr[[8]]),cis.tot.r[4,]))


emeans1 = rbind(c(avgcq[[2]][1],avgcr[[2]][1],avgcq[[4]][1],avgcr[[4]][1]),
               c(avgcq[[6]][1],avgcr[[6]][1],avgcq[[8]][1],avgcr[[8]][1]))
ecislb1 = rbind(c(cicq[[2]][1,1],cicr[[2]][1,1],cicq[[4]][1,1],cicr[[4]][1,1]),
               c(cicq[[6]][1,1],cicr[[6]][1,1],cicq[[8]][1,1],cicr[[8]][1,1]))
ecisub1 = rbind(c(cicq[[2]][2,1],cicr[[2]][2,1],cicq[[4]][2,1],cicr[[4]][2,1]),
               c(cicq[[6]][2,1],cicr[[6]][2,1],cicq[[8]][2,1],cicr[[8]][2,1]))

emeans2 = rbind(c(avgcq[[1]][1],avgcr[[1]][1],avgcq[[3]][1],avgcr[[3]][1]),
               c(avgcq[[5]][1],avgcr[[5]][1],avgcq[[7]][1],avgcr[[7]][1]))
ecislb2 = rbind(c(cicq[[1]][1,1],cicr[[1]][1,1],cicq[[3]][1,1],cicr[[3]][1,1]),
               c(cicq[[5]][1,1],cicr[[5]][1,1],cicq[[7]][1,1],cicr[[7]][1,1]))
ecisub2 = rbind(c(cicq[[1]][2,1],cicr[[1]][2,1],cicq[[3]][2,1],cicr[[3]][2,1]),
               c(cicq[[5]][2,1],cicr[[5]][2,1],cicq[[7]][2,1],cicr[[7]][2,1]))

sigfig = 2

cftab1 = c("Long-Term","","Short-Term","")

for(i in 1:2) {
  for(j in 1:4) {
    cftab1[2*i-1] = paste(cftab1[2*i-1]," & ",round(emeans1[i,j],sigfig),sep="")
  }
  cftab1[2*i-1] = paste(cftab1[2*i-1]," \\\\",sep="")
  for(j in 1:4) {
    cftab1[2*i] = paste(cftab1[2*i]," & [",round(ecislb1[i,j],sigfig),", ",round(ecisub1[i,j],sigfig),"]",sep="")
  }
  cftab1[2*i] = paste(cftab1[2*i]," \\\\",sep="")
}

write.table(cftab1,file="cf_table_depth.tex",quote=FALSE,row.names=FALSE,col.names=FALSE)

cftab2 = c("Long-Term","","Short-Term","")

for(i in 1:2) {
  for(j in 1:4) {
    cftab2[2*i-1] = paste(cftab2[2*i-1]," & ",round(emeans2[i,j],sigfig),sep="")
  }
  cftab2[2*i-1] = paste(cftab2[2*i-1]," \\\\",sep="")
  for(j in 1:4) {
    cftab2[2*i] = paste(cftab2[2*i]," & [",round(ecislb2[i,j],sigfig),", ",round(ecisub2[i,j],sigfig),"]",sep="")
  }
  cftab2[2*i] = paste(cftab2[2*i]," \\\\",sep="")
}

write.table(cftab2,file="cf_table_freq.tex",quote=FALSE,row.names=FALSE,col.names=FALSE)
