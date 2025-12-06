# run counterfactual analysis

# create new counterfactual price data

# create cfdata and cfinit

# flag off price reductions

# we don't want to overwrite stuff like hhmerge so I use attach here

info$pexpect = TRUE  # use data price expectations instead of counterfactual

attach(paste(datapath,"mergedpanel_",product,"_allyears.RData",sep=""))

# check rates of promotions

# use prdf and re-aggregate prices

prdf$WEEK = as.integer(prdf$WEEK)

nweek = length(min(prdf$WEEK):max(prdf$WEEK))
temp = rep(1:4,1+nweek/4)

temp1 = data.frame(as.integer(min(prdf$WEEK):max(prdf$WEEK)),as.integer(temp[1:nweek]))
colnames(temp1) = c("WEEK","month")

prdf = merge(prdf,temp1)
prdf = prdf[order(prdf$IRI_KEY,prdf$WEEK),]

# some brands in prdf get removed in later cuts, so need to fix this

prcols = which(substr(colnames(hhmerge),1,7) == "Product")

predcols = which(substr(colnames(hhmerge),1,2) == "PR")

prdf = prdf[,c("IRI_KEY","WEEK",colnames(hhmerge)[c(prcols,predcols)],"month")]

sourceCpp("find_deals.cpp")

out = findmodes(prdf,as.integer(info$nbrand),FALSE)$modalprice

pct = 0.00

prcols = which(substr(colnames(prdf),1,7) == "Product")

dealprice = prdf[,prcols] < (1-pct)*out & prdf[,prcols] > 0

# use promotion flag at the store-year level

stids = unique(prdf$IRI_KEY)

minweekstub = c(1114,1166,1218,1270,1322,1374,1427)
maxweekstub = c(1165,1217,1269,1321,1373,1426,1478)

prnames = names(hhmerge)[which(substr(colnames(hhmerge),1,7) == "Product")]
prednames = names(hhmerge)[predcols]

j=0
for(prod in prnames) {
  j=j+1
  prdf[[paste("AVN",prednames[[j]],sep="")]] = NA
  prdf[[paste("AVD",prednames[[j]],sep="")]] = NA
  prdf[[paste("PR",prednames[[j]],sep="")]] = NA
  for(st in stids) {
    for(y in 1:7) {
      inds = prdf$IRI_KEY == st & prdf$WEEK >= minweekstub[y] & prdf$WEEK <= maxweekstub[y] & prdf[[prednames[j]]] != -1
      if(sum(inds) > 0) {
        #browser()
        prsale = sum(prdf[[prednames[j]]][inds] == 1)/sum(inds)
        prdf[[paste("PR",prednames[[j]],sep="")]][inds] = prsale
        prdf[[paste("AVN",prednames[[j]],sep="")]][inds] = mean(prdf[[prod]][inds & prdf[[prednames[j]]]==0])
        prdf[[paste("AVD",prednames[[j]],sep="")]][inds] = mean(prdf[[prod]][inds & prdf[[prednames[j]]]==1])
      }
      if(1==1) {
        inds = prdf$IRI_KEY == st & prdf$WEEK >= minweekstub[y] & prdf$WEEK <= maxweekstub[y] & prdf[[prednames[j]]] == -1 & prdf[[prod]] > 0
        if(sum(inds) > 0) {
          prsale = sum(dealprice[inds,j])/sum(inds)
          prdf[[paste("PR",prednames[[j]],sep="")]][inds] = prsale
          prdf[[paste("AVN",prednames[[j]],sep="")]][inds] = mean(prdf[[prod]][inds & !dealprice[,j]])
          prdf[[paste("AVD",prednames[[j]],sep="")]][inds] = mean(prdf[[prod]][inds & dealprice[,j]])
        }
      }
    }
  }
}

# increase depth of current promotions

prdf.old = prdf

chobs = matrix(FALSE,nrow=nrow(prdf),ncol=info$nbrand)

includedprod = prnames[1]

j=0
for(prod in includedprod) {
  j=which(prod == prnames)
  for(st in stids) {
    for(y in 1:7) {
      inds = prdf$IRI_KEY == st & prdf$WEEK >= minweekstub[y] & prdf$WEEK <= maxweekstub[y] & prdf[[prod]] > -1
      if(sum(inds & prdf[[prednames[j]]] != -1) > 0) {
        saleobs = inds & prdf[[prednames[j]]] == 1
        prdf[[prod]][saleobs] = prdf[[prod]][saleobs]*0.5
        chobs[saleobs,j] = TRUE
      } else if (sum(inds & prdf[[prednames[j]]] == -1) > 0) {
        temp = max(prdf[[paste("PR",prednames[[j]],sep="")]][inds & prdf[[prednames[j]]] == -1])
        if(is.finite(temp) & temp > 0) {
          saleobs = inds & dealprice[,j]
          prdf[[prod]][saleobs] = prdf[[prod]][saleobs]*0.5
          chobs[saleobs,j] = TRUE
        }
      }
    }
  }
}


updatepr = apply(chobs,1,max)

prdf$updatepr = updatepr

ii = which(colnames(prdf) == "updatepr")

cfdata = hhmerge[,c("PANID","WEEK")]
cfdata = merge(cfdata,irikeymerge,all.x=TRUE)
cfdata = merge(cfdata,prdf[,c(1:2,prcols,ii)],all.x=TRUE)
cfdata = cfdata[order(cfdata$PANID,cfdata$WEEK),]

cfinit = hhinitmerge[,c("PANID","WEEK")]
cfinit = merge(cfinit,irikeyinitmerge,all.x=TRUE)
cfinit = merge(cfinit,prdf[,c(1:2,prcols,ii)],all.x=TRUE)
cfinit = cfinit[order(cfinit$PANID,cfinit$WEEK),]

# update observations where a change is observed

prcols = which(substr(colnames(hhmerge),1,7) == "Product")

cfdata1 = hhmerge[,c(1:2,prcols)]

for(prod in prnames) {
  cfdata1[[prod]][cfdata1[[prod]]==1000] = -1
}

cfdata2 = cfdata1

prcolsi = which(substr(colnames(hhinitmerge),1,7) == "Product")

cfinit1 = hhinitmerge[,c(1:2,prcolsi)]
for(prod in prnames) {
  cfinit1[[prod]][cfinit1[[prod]]==1000] = -1
}
cfinit2 = cfinit1

cfdata1$updatepr = cfdata$updatepr
cfinit1$updatepr = cfinit$updatepr

# this may not matter but I'm doing it to be safe

cfdata1 = cfdata1[!as.logical(cfdata1$updatepr),]
cfinit1 = cfinit1[!as.logical(cfinit1$updatepr),]

cfdata = cfdata[as.logical(cfdata$updatepr),]
cfinit = cfinit[as.logical(cfinit$updatepr),]

cfdata$IRI_KEY = NULL
cfinit$IRI_KEY = NULL

cfdata = rbind(cfdata,cfdata1)
cfinit = rbind(cfinit,cfinit1)

cfdata = cfdata[order(cfdata$PANID,cfdata$WEEK),]
cfinit = cfinit[order(cfinit$PANID,cfinit$WEEK),]

# sometimes for other products we see changes.  Not sure why but let's stick to products where there
# are changes

# replace missing when missing in original data

for(prod in includedprod) {
  cfdata[[prod]] = pmin(cfdata[[prod]],cfdata2[[prod]])
  cfinit[[prod]] = pmin(cfinit[[prod]],cfinit2[[prod]])
}

for(prod in setdiff(prnames,includedprod)) {
  cfdata[[prod]] = cfdata2[[prod]]
  cfinit[[prod]] = cfinit2[[prod]]
}

# replace -1 with 1000

for(prod in prnames) {
  cfdata[[prod]][cfdata[[prod]]==-1] = 1000
  cfinit[[prod]][cfinit[[prod]]==-1] = 1000
}

cfdata$updatepr = NULL
cfinit$updatepr = NULL

cfdata = as.matrix(cfdata[,3:ncol(cfdata)],nrow=nrow(cfdata),ncol=info$nbrand)
cfinit = as.matrix(cfinit[,3:ncol(cfinit)],nrow=nrow(cfinit),ncol=info$nbrand)

hhmerge[,prcols] = cfdata  # I think we have to do this or the indexing below messes up
save(cfdata,file=paste("cfdata",cf.index,".RData",sep=""))
stop()
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
if(!info$pexpect) {
  info$rgrid = t(as.matrix(pricetab[psampleinds,1:info$nbrand]))

  info$impdist = pricetab$prob[psampleinds]
}
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
fit.index = 3
for(i in 1:3) {
  load(paste("out_ml_A_fix1_",fit.index,"_",i,".RData",sep=""))
  inputs$pstart[c(info$nbrand+5,info$nbrand+6),crassign==i] = out.ml$par[1:2]
}

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


