# Create estimation datasets

set.seed(813345)

yearlist = paste("Year",1:7,sep="")

weekstub = c("1114_1165","1166_1217","1218_1269","1270_1321","1322_1373","1374_1426","1427_1478")
minweekstub = c(1114,1166,1218,1270,1322,1374,1427)
maxweekstub = c(1165,1217,1269,1321,1373,1426,1478)

ystub = c("1","2","3","4","5","6","7")

yind = 0
for(yearnum in yearlist) {
  yind = yind + 1

  load(paste(outputpath,"brands_",yearnum,"_",product,".RData",sep=""))

  brandinfo = brandinfo[,-7]

  load(paste(outputpath,"hhpanel_",product,"_",yearnum,".RData",sep=""))

  hhpanel = hhpanel[,-15]

  load(paste(outputpath,"stores_",yearnum,"_",product,".RData",sep=""))

  if(yind == 1) {
    binfo.big = brandinfo
    hhpanel.big = hhpanel
    stdata.big = stdata
  } else {
    if(yind == 7) {binfo.big = rbind(binfo.big,brandinfo)}
    hhpanel.big = rbind(hhpanel.big,hhpanel)
    stdata.big = rbind(stdata.big,stdata)
  }

}

brandinfo = binfo.big
hhpanel = hhpanel.big
stdata = stdata.big
rm(list=c("binfo.big","hhpanel.big","stdata.big"))
# take store data and create UPCs and prices
# generate upc variable

sy = as.character(stdata$SY)
ge = as.character(stdata$GE)
vend = as.character(stdata$VEND)
item = as.character(stdata$ITEM)

lvend = nchar(vend)   #this is always 5 so we're good here
litem = nchar(item)
for(i in 1:4) {
  item[litem==i] = paste(paste(rep("0",5-i),collapse=""),item[litem==i],sep="")
}
lsy = nchar(sy)
sy[lsy==1] = paste("0",sy[lsy==1],sep="")
upc1 = paste(sy,"0",ge,vend,item,sep="")

# compute price per unit and merge in brand data

ppunit = stdata$DOLLARS/stdata$UNITS

stdata$ppunit = ppunit
stdata$COLUPC = upc1

hhpanel$ppunit = hhpanel$DOLLARS/hhpanel$UNITS

stdata$SY <- NULL
stdata$GE <- NULL
stdata$VEND <- NULL
stdata$ITEM <- NULL

# remove duplicate UPCS from brandinfo so the merge doesn't mess up

brandinfo = brandinfo[order(brandinfo$COLUPC),]
lupc = c("0",brandinfo$COLUPC[1:(nrow(brandinfo)-1)])
brandinfo = brandinfo[brandinfo$COLUPC != lupc,]

#stold = stdata

stdata = merge(stdata,brandinfo)  # in laundry detergent we lose 3 UPCS here. initial 2 digit code is 08 - what are these?

# remove powders
if(remove.powder) {
  update.droplist("hhsmall","bought powders",hhpanel,hhpanel[hhpanel$FORM == "LIQUID",],"PANID")
  hhpanel = hhpanel[hhpanel$FORM == "LIQUID",]
}

# array of price data

sizetab = table(hhpanel$VOL_EQ)
sizetab = sizetab[order(sizetab,decreasing=TRUE)]
vols = as.numeric(names(sizetab))

brandtab = table(hhpanel$L5)
brandtab = brandtab[order(brandtab,decreasing=TRUE)]

nweek = length(unique(hhpanel$WEEK))

allweek = unique(stdata$WEEK)
allweek = allweek[order(allweek)]

allstore = union(unique(stdata$IRI_KEY),unique(hhpanel$IRI_KEY))
allstore = allstore[order(allstore)]

nstore = length(allstore)

# first, check through all the brand and size combos to make sure that all are actually available
# something is not then drop it

brsize = matrix(0,nrow=1,ncol=2)

index = 0

for(b in 1:nbrand) {
  for(s in 1:nsize) {
    if(sum(stdata$L5 == names(brandtab)[b] & stdata$VOL_EQ == vols[s]) > 0) {
      index = index + 1
      if(index == 1) {
        brsize[1,1] = b
        brsize[1,2] = s
      } else {
        index = index + 1
        brsize = rbind(brsize,c(b,s))
      }
    }
  }
}

nbrsize = nrow(brsize)

stdata$brindex = 0
for(i in 1:nbrsize) {
  stdata$brindex[stdata$L5 == names(brandtab)[brsize[i,1]] & stdata$VOL_EQ == vols[brsize[i,2]]] = i
}

stdata$wkindex = 0
for(t in 1:nweek) {
  stdata$wkindex[stdata$WEEK == allweek[t]] = t
}

stdata$stindex = 0
for(s in 1:nstore) {
  stdata$stindex[stdata$IRI_KEY == allstore[s]] = s
}

meanprice.st = aggregate(stdata$ppunit,by=list(stdata$brindex,stdata$wkindex,stdata$stindex),FUN=mean)

colnames(meanprice.st) = c("brindex","wkindex","stindex","ppunit")

maxfeat.st = aggregate(stdata$F != "NONE",by=list(stdata$brindex,stdata$wkindex,stdata$stindex),FUN=max)

colnames(maxfeat.st) = c("brindex","wkindex","stindex","Feat")

maxdisp.st = aggregate(stdata$D,by=list(stdata$brindex,stdata$wkindex,stdata$stindex),FUN=max)

colnames(maxdisp.st) = c("brindex","wkindex","stindex","Disp")

maxpr.st = aggregate(stdata$PR,by=list(stdata$brindex,stdata$wkindex,stdata$stindex),FUN=max)

colnames(maxpr.st) = c("brindex","wkindex","stindex","PR")

stdata$stindex <- NULL
stdata$wkindex <- NULL

hhpanel$brindex = 0
for(i in 1:nbrsize) {
  hhpanel$brindex[hhpanel$L5 == names(brandtab)[brsize[i,1]] & hhpanel$VOL_EQ == vols[brsize[i,2]]] = i
}

hhpanel$wkindex = 0
for(t in 1:nweek) {
  hhpanel$wkindex[hhpanel$WEEK == allweek[t]] = t
}

hhpanel$stindex = 0
for(s in 1:nstore) {
  hhpanel$stindex[hhpanel$IRI_KEY == allstore[s]] = s
}

meanprice.hh = aggregate(hhpanel$ppunit,by=list(hhpanel$brindex,hhpanel$wkindex,hhpanel$stindex),FUN=mean)

colnames(meanprice.hh) = c("brindex","wkindex","stindex","ppunit")

hhpanel$brindex <- NULL
hhpanel$wkindex <- NULL
hhpanel$stindex <- NULL

# fill in price from the store data

if(loadprices) {
  load(paste(outputpath,"pricesave.RData",sep=""))
} else {
  if(use.ccode) {

    sourceCpp("aggregate_prices.cpp")

    out = agprices(meanprice.st, meanprice.hh, maxfeat.st, maxdisp.st, maxpr.st, as.integer(nbrsize), as.integer(nweek), as.integer(nstore))

    prarray = out[[1]]
    prchange = out[[2]]
    farray = out[[3]]
    darray = out[[4]]
    tprarray = out[[5]]
    rm(out)

  } else {

    # this is extremely slow - use the C++ version instead

    prarray = array(-1,c(nweek,nstore,nbrsize))
    prchange = array(-1,c(nweek,nstore,nbrsize))

    farray = array(-1,c(nweek,nstore,nbrsize))
    darray = array(-1,c(nweek,nstore,nbrsize))
    tprarray = array(-1,c(nweek,nstore,nbrsize))

    for(i in 1:nbrsize) {
      for(t in 1:nweek) {
        for(s in 1:nstore) {
          ind1 = which(meanprice.st$brindex == i & meanprice.st$wkindex == t & meanprice.st$stindex == s)
          if(length(ind1) > 0) {
            prarray[t,s,i] = meanprice.st$ppunit[ind1]
            farray[t,s,i] = maxfeat.st$Feat[ind1]
            darray[t,s,i] = maxdisp.st$Disp[ind1]
            tprarray[t,s,i] = maxpr.st$PR[ind1]
          }

          ind1 = which(meanprice.hh$brindex == i & meanprice.hh$wkindex == t & meanprice.hh$stindex == s)

          if(length(ind1) > 0) {
            pnew = meanprice.hh$ppunit[ind1]
            prchange[t,s,i] = prarray[t,s,i] - pnew
            prarray[t,s,i] = pnew
          }

        }
      }
    }

  }

  # fill in missing prices from household panel

  # impute missing prices
  # there might be a better way to do this... probably want to check into this

  for(i in 1:nbrsize) {
    for(s in 1:nstore) {
      # if missing in period 1, find first nonmissing price and impute backwards
      # if the product is always missing, we will assume it is not available
      if(prarray[1,s,i] == -1) {
        if(sum(prarray[,s,i]==-1) != nweek) {
          t=2
          while(prarray[t,s,i] == -1) {
            t=t+1
          }
          prarray[1:(t-1),s,i] = prarray[t,s,i]
        }
      }
      t=2
      while(t <= nweek-1) {
        increment = 1
        if(prarray[t,s,i] == -1) {
          pold = prarray[t-1,s,i]
          while(t+increment < nweek & prarray[t+increment,s,i] == -1) {
            increment=increment+1
            #cat(increment,"\n")
          }
          if(t+increment == nweek & prarray[t+increment,s,i] == -1) {
            if(imprule==1) {
              prarray[t:nweek,s,i] = pold
            }
          } else {
            if(imprule==1) {
              pnext = prarray[t+increment,s,i]
              prarray[t:(t+increment-1),s,i] = 0.5*(pold+pnext)
            }
          }
          increment=increment+1
        }
        t=t+increment
      }
    }
  }

  save(prarray,prchange,farray,darray,tprarray,file=paste(outputpath,"pricesave.RData",sep=""))

}



# now that we have prices, create choice data we can use for estimation

# first, remove households who make purchases that are outside the acceptable choice set
# we may want to include other brands as an aggregate and only exclude bad sizes
# on a quick look this doesn't make a big difference - roughly 130 households

hhinit = hhpanel[hhpanel$WEEK < minweekstub[min(estyears)],]
hhpanel = hhpanel[hhpanel$WEEK >= minweekstub[min(estyears)],]

goodbrand = hhpanel$L5 %in% names(brandtab[brsize[,1]])
goodsize = hhpanel$VOL_EQ %in% vols[brsize[,2]]

goodpurch = goodbrand & goodsize

hhgoodpurch = aggregate(goodpurch,by=list(hhpanel$PANID),FUN=min)

hhsmall = hhpanel

hhsmall = merge(hhsmall,hhgoodpurch,by.x="PANID",by.y="Group.1")

hhsmall = hhsmall[hhsmall$x==1,]
hhsmall$x <- NULL

update.droplist("hhsmall","bought outside choice set",hhpanel,hhsmall,"PANID")

hhsmall = hhsmall[,c("IRI_KEY","WEEK","PANID","UNITS","DOLLARS","Trip_Count","L5","VOL_EQ",demovars)]

hhsmall$BRAND = hhsmall$L5
hhsmall$L5 <- NULL

# next, check into multiple purchases within a trip

hhsmall$brindex = 0

for(i in 1:nrow(brsize)) {
  hhsmall$brindex[hhsmall$BRAND == names(brandtab)[brsize[i,1]] & hhsmall$VOL_EQ == vols[brsize[i,2]]] = i
}

# remove households that purchase a brand not available in the store sample - very rare

badid = unique(hhsmall$PANID[hhsmall$brindex == 0])

update.droplist("hhsmall","bought brand not available in store sample",hhsmall,hhsmall[!(hhsmall$PANID %in% badid),],"PANID")

hhsmall = hhsmall[!(hhsmall$PANID %in% badid),]

brmin = aggregate(hhsmall$brindex,by=list(hhsmall$PANID,hhsmall$WEEK),FUN=min)
brmax = aggregate(hhsmall$brindex,by=list(hhsmall$PANID,hhsmall$WEEK),FUN=max)

badid = unique(hhsmall$PANID[brmin$x != brmax$x]) #guys who buy multiple brands in one week

update.droplist("hhsmall","bought multiple brands in a single week",hhsmall,hhsmall[!(hhsmall$PANID %in% badid),],"PANID")

hhsmall = hhsmall[!(hhsmall$PANID %in% badid),]

stmin = aggregate(hhsmall$IRI_KEY,by=list(hhsmall$PANID,hhsmall$WEEK),FUN=min)
stmax = aggregate(hhsmall$IRI_KEY,by=list(hhsmall$PANID,hhsmall$WEEK),FUN=max)

weeklen1 = aggregate(hhsmall$PANID,by=list(hhsmall$PANID,hhsmall$WEEK),FUN=length)

badid = union(badid,unique(hhsmall$PANID[stmin$x != stmax$x])) #guys who make purchases from multiple stores in one week

update.droplist("hhsmall","bought from different stores in a single week",hhsmall,hhsmall[!(hhsmall$PANID %in% badid),],"PANID")

hhsmall = hhsmall[!(hhsmall$PANID %in% badid),]

# aggregate up multiple purchases of same product into single purchase

purchsum = aggregate(cbind(hhsmall$UNITS,hhsmall$DOLLARS),by=list(hhsmall$PANID,hhsmall$WEEK),FUN=sum)

colnames(purchsum) = c("PANID","WEEK","totunits","totdollars")

othermax = aggregate(cbind(hhsmall$IRI_KEY,hhsmall[,demovars],
                           hhsmall$brindex),by=list(hhsmall$PANID,hhsmall$WEEK),FUN=max)

colnames(othermax) = c("PANID","WEEK","IRI_KEY",demovars,"brindex")

purchsum = purchsum[order(purchsum$PANID,purchsum$WEEK),]
othermax = othermax[order(othermax$PANID,othermax$WEEK),]

hhsmall = cbind(purchsum,othermax[,3:ncol(othermax)])

hhsmall$ppunit = hhsmall$totdollars/hhsmall$totunits

hhsmall$flag = 1

# next, merge in data on household trips

inittrip = NULL
yind = 0
for(yr in ystub) {
  if(as.integer(yr) %in% estyears) {
    yind = yind + 1
    hhtrip = read.csv(paste(rawdatapath,"demos trips external/trips",yr," jul08.csv",sep=""),stringsAsFactors=FALSE)
    if(is.null(hhtrip$IRI_KEY)) {
      if(!is.null(hhtrip$IRI_Key)) {
        names(hhtrip)[names(hhtrip)=="IRI_Key"] <- "IRI_KEY"
      }
    }
    if(!is.null(hhtrip$KRYSCENTS)) {
      hhtrip$KRYSCENTS <- NULL
    }
    if(yind == 1) {
      hhtrip.big = hhtrip
    } else {
      hhtrip.big = rbind(hhtrip.big,hhtrip)
    }
  } else {
    itrip = read.csv(paste(rawdatapath,"demos trips external/trips",yr," jul08.csv",sep=""),stringsAsFactors=FALSE)
    if(is.null(itrip$IRI_KEY)) {
      if(!is.null(itrip$IRI_Key)) {
        names(itrip)[names(itrip)=="IRI_Key"] <- "IRI_KEY"
      }
    }
    if(!is.null(itrip$KRYSCENTS)) {
      itrip$KRYSCENTS <- NULL
    }
    inittrip = rbind(inittrip,itrip)
  }
}

hhtrip = hhtrip.big
rm(hhtrip.big)

allid = unique(hhpanel$PANID)
hhtrip = hhtrip[hhtrip$PANID %in% allid,]
inittrip = inittrip[inittrip$PANID %in% allid,]

if(is.null(hhtrip$IRI_KEY)) {
  if(!is.null(hhtrip$IRI_Key)) {
    hhtrip$IRI_KEY = hhtrip$IRI_Key
    hhtrip$IRI_Key = NULL
  }
}

if(is.null(inittrip$IRI_KEY)) {
  if(!is.null(inittrip$IRI_Key)) {
    inittrip$IRI_KEY = inittrip$IRI_Key
    inittrip$IRI_Key = NULL
  }
}

tripminiri2 = aggregate(hhtrip$IRI_KEY,by=list(hhtrip$PANID,hhtrip$WEEK),FUN=min)
tripmaxiri2 = aggregate(hhtrip$IRI_KEY,by=list(hhtrip$PANID,hhtrip$WEEK),FUN=max)

weeklen2 = aggregate(hhtrip$PANID,by=list(hhtrip$PANID,hhtrip$WEEK),FUN=length)

minminiri2 = aggregate(hhtrip$MINUTE,by=list(hhtrip$PANID,hhtrip$WEEK),FUN=min)
minmaxiri2 = aggregate(hhtrip$MINUTE,by=list(hhtrip$PANID,hhtrip$WEEK),FUN=max)

# collapse down trips that occur in the same store and same week
# I don't think we can match trips to individual purchases

hhtrip.ag = aggregate(rep(0,nrow(hhtrip)),by=list(hhtrip$PANID,hhtrip$WEEK,hhtrip$IRI_KEY),FUN=sum)

colnames(hhtrip.ag) = c("PANID","WEEK","IRI_KEY","x")
hhtrip.ag$x <- NULL

hhtrip.ag = hhtrip.ag[order(hhtrip.ag$PANID,hhtrip.ag$WEEK,hhtrip.ag$IRI_KEY),]

tripminiri = aggregate(hhtrip.ag$IRI_KEY,by=list(hhtrip.ag$PANID,hhtrip.ag$WEEK),FUN=min)
tripmaxiri = aggregate(hhtrip.ag$IRI_KEY,by=list(hhtrip.ag$PANID,hhtrip.ag$WEEK),FUN=max)

hhmerge = merge(hhtrip.ag,hhsmall,all.x=TRUE)



inittrip.ag = aggregate(rep(0,nrow(inittrip)),by=list(inittrip$PANID,inittrip$WEEK,inittrip$IRI_KEY),FUN=sum)

colnames(inittrip.ag) = c("PANID","WEEK","IRI_KEY","x")
inittrip.ag$x <- NULL

inittrip.ag = inittrip.ag[order(inittrip.ag$PANID,inittrip.ag$WEEK,inittrip.ag$IRI_KEY),]

inittripminiri = aggregate(inittrip.ag$IRI_KEY,by=list(inittrip.ag$PANID,inittrip.ag$WEEK),FUN=min)
inittripmaxiri = aggregate(inittrip.ag$IRI_KEY,by=list(inittrip.ag$PANID,inittrip.ag$WEEK),FUN=max)

hhinit$flag = 1

hhinitmerge = merge(inittrip.ag,hhinit[,c("IRI_KEY","WEEK","PANID","UNITS","DOLLARS","Trip_Count","L5","VOL_EQ",demovars,"flag")],all.x=TRUE)

hhinit$flag = NULL

# identify the first store visited.  since initial prices are just
# for simulation just pick a random one

#i.ag1 = c(0,inittrip.ag$PANID[1:(nrow(inittrip.ag)-1)])
#i.ag2 = c(0,inittrip.ag$WEEK[1:(nrow(inittrip.ag)-1)])

#inittrip.ag1 = inittrip.ag[inittrip.ag$PANID != i.ag1 | inittrip.ag$WEEK != i.ag2,]

#hhinitmerge = merge(hhinit,inittrip.ag1,all.x=TRUE)

# a lot of panel ids are in hhmerge, but not hhsmall.  Remove them here
# leaving them in makes putting the data together a lot slower.

hhmerge = hhmerge[hhmerge$PANID %in% unique(hhsmall$PANID),]

hhinitmerge = hhinitmerge[hhinitmerge$PANID %in% unique(hhsmall$PANID),]

# remove observations where a household bought from one store but visited another

hhmerge$totunits[is.na(hhmerge$totunits)] = 0

# these should be the same

stboughtmin = aggregate(hhtrip.ag$IRI_KEY[hhmerge$totunits>0],by=list(hhtrip.ag$PANID[hhmerge$totunits>0],hhtrip.ag$WEEK[hhmerge$totunits>0]),FUN=min)
stboughtmax = aggregate(hhtrip.ag$IRI_KEY[hhmerge$totunits>0],by=list(hhtrip.ag$PANID[hhmerge$totunits>0],hhtrip.ag$WEEK[hhmerge$totunits>0]),FUN=max)

colnames(stboughtmin) = c("PANID","WEEK","IRI_BOUGHT")

hhmerge = merge(hhmerge,stboughtmin,all.x=TRUE)

hhmerge = hhmerge[is.na(hhmerge$IRI_BOUGHT) | hhmerge$IRI_BOUGHT == hhmerge$IRI_KEY, ]

hhmerge$IRI_BOUGHT <- NULL

prdf <- data.frame(rep(allstore,nweek),kronecker(allweek,rep(1,nstore)),matrix(0,nrow=nweek*nstore,ncol=nrow(brsize)*4))
colnames(prdf) = c("IRI_KEY","WEEK",paste("Product",brsize[,1],"_",brsize[,2],sep=""),
                   paste("Feat",brsize[,1],"_",brsize[,2],sep=""),
                   paste("Disp",brsize[,1],"_",brsize[,2],sep=""),
                   paste("PR",brsize[,1],"_",brsize[,2],sep=""))

for(i in 1:nweek) {
  prdf[prdf$WEEK == allweek[i],2+1:nrow(brsize)] = prarray[i,,]
  prdf[prdf$WEEK == allweek[i],2+nrow(brsize)+1:nrow(brsize)] = farray[i,,]
  prdf[prdf$WEEK == allweek[i],2+2*nrow(brsize)+1:nrow(brsize)] = darray[i,,]
  prdf[prdf$WEEK == allweek[i],2+3*nrow(brsize)+1:nrow(brsize)] = tprarray[i,,]
}

prdf = prdf[order(prdf$IRI_KEY,prdf$WEEK),]

# some storeids appear in the household panel that are not in the store panel. Why?
# what should be done about this?
# for now, let's only keep households that shop at the stores in the store panel

goodstore = hhmerge$IRI_KEY %in% allstore

goodhhstore = aggregate(goodstore,by=list(hhmerge$PANID),FUN=min)

update.droplist("hhmerge","bought from store outside the store panel",hhmerge,hhmerge[hhmerge$PANID %in% goodhhstore[goodhhstore$x==1,1],],"PANID")

hhmerge = hhmerge[hhmerge$PANID %in% goodhhstore[goodhhstore$x==1,1],]

hhmerge = merge(hhmerge,prdf)

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]



hhinitmerge$BRAND = hhinitmerge$L5

hhinitmerge$brindex = 0

for(i in 1:nrow(brsize)) {
  hhinitmerge$brindex[hhinitmerge$BRAND == names(brandtab)[brsize[i,1]] & hhinitmerge$VOL_EQ == vols[brsize[i,2]]] = i
}

hhinitmerge$totunits = hhinitmerge$UNITS

hhinitmerge$totdollars = hhinitmerge$DOLLARS

hhinitmerge$ppunit = hhinitmerge$totdollars/hhinitmerge$totunits

#hhinitmerge$L5 <- NULL

#hhinitmerge$UNITS <- NULL

#hhinitmerge$DOLLARS <- NULL

hhinitmerge = hhinitmerge[,c("WEEK",        "IRI_KEY",     "PANID",       "totunits",    "totdollars",
   "HH_INC",      "HH_SIZE",     "HH_RACE",     "HH_AGE",      "HH_EDU",
   "HH_OCC",      "maleedu",     "fedu",        "child",       "brindex",
   "ppunit",      "flag")]

hhinitmerge = hhinitmerge[hhinitmerge$PANID %in% goodhhstore[goodhhstore$x==1,1],]

hhinitmerge = merge(hhinitmerge,prdf)

hhinitmerge = hhinitmerge[order(hhinitmerge$PANID,hhinitmerge$WEEK),]

# fill in missing values of hhmerge

hhmerge$brindex[is.na(hhmerge$brindex)] = 0

hhmerge$totunits[is.na(hhmerge$totunits)] = 0

hhmerge$totdollars[is.na(hhmerge$totdollars)] = 0

hhmerge$ppunit[is.na(hhmerge$ppunit)] = -1

hhmerge$flag[is.na(hhmerge$flag)] = 0

for(d in demovars) {

  hhmerge[[d]][is.na(hhmerge[[d]])] = -1

}


hhmerge$totunits[is.na(hhmerge$totunits)] = 0

hhmerge$totdollars[is.na(hhmerge$totdollars)] = 0

hhinitmerge$ppunit[is.na(hhinitmerge$ppunit)] = -1

hhinitmerge$flag[is.na(hhinitmerge$flag)] = 0

for(d in demovars) {

  hhinitmerge[[d]][is.na(hhinitmerge[[d]])] = -1

}

cat("Price differences before taking min prices.\n")

for(i in 1:nrow(brsize)) {
  rows = hhmerge$brindex == i
  cat("Prices same for product ",i,": ",sum(hhmerge$ppunit[rows]==hhmerge[rows,11+i]),"\n",sep="")
  cat("HH price > for product ",i,": ",sum(hhmerge$ppunit[rows]>hhmerge[rows,11+i]),"\n",sep="")
  cat("Store price > for product ",i,": ",sum(hhmerge$ppunit[rows]<hhmerge[rows,11+i]),"\n",sep="")
  cat("Differences: \n")
  print(summary(hhmerge$ppunit[rows]-hhmerge[rows,11+i]))
}

# compute minimum price hh faces across stores

hhmerge1 = hhmerge
hhinitmerge1 = hhinitmerge
if(pricemethod == 1) {
  for(i in ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)) {
    hhmerge1[hhmerge1[,i]==-1,i] = 1000
  }

  #ncol(hhmerge1)-seq(nrow(brsize)*4-1,0,-1)

  minprice = aggregate(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=min)
} else if (pricemethod == 2) {
  minprice = aggregate(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=max)
} else if (pricemethod == 3) {
  for(i in ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)) {
    hhmerge1[hhmerge1[,i]==-1,i] = 1000
  }

  nobsbyhh = aggregate(hhmerge[,1],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=length)
  colnames(nobsbyhh) = c("PANID","WEEK","nobs")
  nobsbyhh$u = runif(nrow(nobsbyhh))

  hhmerge1 = merge(hhmerge1,nobsbyhh)

  hhmerge1 = hhmerge1[order(hhmerge1$PANID,hhmerge1$WEEK),]

  obsnum1 = 1:nrow(hhmerge1)

  lid1 = c(0,hhmerge1$PANID[1:(nrow(hhmerge1)-1)])
  wid1 = c(0,hhmerge1$WEEK[1:(nrow(hhmerge1)-1)])

  fobs1 = hhmerge1$PANID != lid1 | hhmerge1$WEEK != wid1

  fobsnum1 = obsnum1*as.numeric(fobs1)+(1-as.numeric(fobs1))*(nrow(hhmerge1)+1)

  fobsnum11 = aggregate(fobsnum1,by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=min)

  colnames(fobsnum11) = c("PANID","WEEK","firstidnum")

  hhmerge1 = merge(hhmerge1,fobsnum11)
  hhmerge1 = hhmerge1[order(hhmerge1$PANID,hhmerge1$WEEK),]

  hhmerge1$obsnum = obsnum1 - hhmerge1$firstidnum

  inrange = hhmerge1$u >= hhmerge1$obsnum/hhmerge1$nobs & hhmerge1$u < (hhmerge1$obsnum + 1)/hhmerge1$nobs

  hhmerge1$u <- NULL
  hhmerge1$firstidnum <- NULL
  hhmerge1$obsnum <- NULL
  hhmerge1$nobs <- NULL

  minprice = hhmerge1[inrange,c("PANID","WEEK",colnames(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))]

} else if (pricemethod == 4) {

} else if (pricemethod == 5) {

  # figure out how often different stores are visited

  stvisits = aggregate(hhmerge1$IRI_KEY,by=list(hhmerge1$PANID,hhmerge1$IRI_KEY),FUN=length)

  colnames(stvisits) = c("PANID","IRI_KEY","nvisits")

  hhmerge1 = merge(hhmerge1,stvisits)

  hhmerge1 = hhmerge1[order(hhmerge1$PANID,hhmerge1$WEEK,hhmerge1$nvisits),]

  lid1 = c(0,hhmerge1$PANID[1:(nrow(hhmerge1)-1)])
  wid1 = c(0,hhmerge1$WEEK[1:(nrow(hhmerge1)-1)])

  fobs1 = hhmerge1$PANID != lid1 | hhmerge1$WEEK != wid1

  hhmerge1$nvisits <- NULL

  minprice = hhmerge1[fobs1,c("PANID","WEEK",colnames(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))]

}  else if (pricemethod == 6) {

  # figure out how often different stores are visited

  stvisits = aggregate(hhmerge1$IRI_KEY[hhmerge$brindex > 0],by=list(hhmerge1$PANID[hhmerge$brindex > 0],hhmerge1$IRI_KEY[hhmerge$brindex > 0]),FUN=length)

  colnames(stvisits) = c("PANID","IRI_KEY","nvisits")

  hhmerge1 = merge(hhmerge1,stvisits,all.x=TRUE)

  hhmerge1 = hhmerge1[order(hhmerge1$PANID,hhmerge1$WEEK,hhmerge1$nvisits),]

  lid1 = c(0,hhmerge1$PANID[1:(nrow(hhmerge1)-1)])
  wid1 = c(0,hhmerge1$WEEK[1:(nrow(hhmerge1)-1)])

  fobs1 = hhmerge1$PANID != lid1 | hhmerge1$WEEK != wid1

  hhmerge1$nvisits <- NULL

  minprice = hhmerge1[fobs1,c("PANID","WEEK",colnames(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))]

  irikeymerge = hhmerge1[fobs1,c("PANID","WEEK","IRI_KEY")]

  stvisitsinit = aggregate(hhinitmerge1$IRI_KEY[hhinitmerge$brindex > 0],by=list(hhinitmerge1$PANID[hhinitmerge$brindex > 0],hhinitmerge1$IRI_KEY[hhinitmerge$brindex > 0]),FUN=length)

  colnames(stvisitsinit) = c("PANID","IRI_KEY","nvisits")

  hhinitmerge1 = merge(hhinitmerge1,stvisitsinit,all.x=TRUE)

  hhinitmerge1 = hhinitmerge1[order(hhinitmerge1$PANID,hhinitmerge1$WEEK,hhinitmerge1$nvisits),]

  lid1 = c(0,hhinitmerge1$PANID[1:(nrow(hhinitmerge1)-1)])
  wid1 = c(0,hhinitmerge1$WEEK[1:(nrow(hhinitmerge1)-1)])

  fobs1 = hhinitmerge1$PANID != lid1 | hhinitmerge1$WEEK != wid1

  hhinitmerge1$nvisits <- NULL

  minpriceinit = hhinitmerge1[fobs1,c("PANID","WEEK",colnames(hhinitmerge1[,ncol(hhinitmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))]

  irikeyinitmerge = hhinitmerge1[fobs1,c("PANID","WEEK","IRI_KEY")]

}

colnames(minprice) = c("PANID","WEEK",colnames(hhmerge1[,ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))

# compute max of feature, display and price promotion

maxadv = aggregate(hhmerge1[,ncol(hhmerge1)-seq(nrow(brsize)*3-1,0,-1)],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=max)

colnames(maxadv) = c("PANID","WEEK",colnames(hhmerge1[,ncol(hhmerge1)-seq(nrow(brsize)*3-1,0,-1)]))

othervars = aggregate(hhmerge1[,c("brindex","ppunit")],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=max)
colnames(othervars) = c("PANID","WEEK","brindex","ppunit")

totunits = aggregate(hhmerge1[,c("totunits")],by=list(hhmerge1$PANID,hhmerge1$WEEK),FUN=sum)
colnames(totunits) = c("PANID","WEEK","totunits")

minprice = minprice[order(minprice$PANID,minprice$WEEK),]
othervars = othervars[order(othervars$PANID,othervars$WEEK),]
totunits = totunits[order(totunits$PANID,totunits$WEEK),]

hhmerge1 = cbind(minprice[,c("PANID","WEEK")],othervars[c("brindex","ppunit")],minprice[,(ncol(minprice)-nrow(brsize)+1):ncol(minprice)],
                 maxadv[,(ncol(maxadv)-nrow(brsize)*3+1):ncol(maxadv)],totunits[,"totunits"])

if(pricemethod == 2) {
  for(i in ncol(hhmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)) {
    hhmerge1[hhmerge1[,i]==-1,i] = 1000
  }
}

demos1 = aggregate(hhmerge[,c(demovars,"flag")],by=list(hhmerge$PANID),FUN=max,na.rm=TRUE)

colnames(demos1) = c("PANID",demovars,"flag")

for(d in demovars) {
  hhmerge[[d]] <- NULL
}

hhmerge$flag <- NULL

hhmerge = merge(hhmerge1,demos1)

hhmerge = hhmerge[hhmerge$flag==1,]


colnames(minpriceinit) = c("PANID","WEEK",colnames(hhinitmerge1[,ncol(hhinitmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)]))

# compute max of feature, display and price promotion

maxadvinit = aggregate(hhinitmerge1[,ncol(hhinitmerge1)-seq(nrow(brsize)*3-1,0,-1)],by=list(hhinitmerge1$PANID,hhinitmerge1$WEEK),FUN=max)

colnames(maxadvinit) = c("PANID","WEEK",colnames(hhinitmerge1[,ncol(hhinitmerge1)-seq(nrow(brsize)*3-1,0,-1)]))

othervarsinit = aggregate(hhinitmerge1[,c("brindex","ppunit")],by=list(hhinitmerge1$PANID,hhinitmerge1$WEEK),FUN=max)
colnames(othervarsinit) = c("PANID","WEEK","brindex","ppunit")

totunitsinit = aggregate(hhinitmerge1[,c("totunits")],by=list(hhinitmerge1$PANID,hhinitmerge1$WEEK),FUN=sum)
colnames(totunitsinit) = c("PANID","WEEK","totunits")

minpriceinit = minpriceinit[order(minpriceinit$PANID,minpriceinit$WEEK),]
othervarsinit = othervarsinit[order(othervarsinit$PANID,othervarsinit$WEEK),]
totunitsinit = totunitsinit[order(totunitsinit$PANID,totunitsinit$WEEK),]

hhinitmerge1 = cbind(minpriceinit[,c("PANID","WEEK")],othervarsinit[c("brindex","ppunit")],minpriceinit[,(ncol(minpriceinit)-nrow(brsize)+1):ncol(minpriceinit)],
                 maxadvinit[,(ncol(maxadvinit)-nrow(brsize)*3+1):ncol(maxadvinit)],totunitsinit[,"totunits"])

if(pricemethod == 2) {
  for(i in ncol(hhinitmerge1)-3*nrow(brsize)-seq(nrow(brsize)-1,0,-1)) {
    hhinitmerge1[hhinitmerge1[,i]==-1,i] = 1000
  }
}

demosinit1 = aggregate(hhinitmerge[,c(demovars,"flag")],by=list(hhinitmerge$PANID),FUN=max,na.rm=TRUE)

colnames(demosinit1) = c("PANID",demovars,"flag")

for(d in demovars) {
  hhinitmerge[[d]] <- NULL
}

hhinitmerge$flag <- NULL

hhinitmerge = merge(hhinitmerge1,demosinit1)

hhinitmerge = hhinitmerge[hhinitmerge$flag==1,]

# check on differences between store and panel prices
# also fix up the prices so that the prices the household actually paid match the store price

cat("\n")
cat("Price differences before after min prices.\n")

for(i in 1:nrow(brsize)) {
  rows = hhmerge$brindex == i
  cat("Prices same for product ",i,": ",sum(hhmerge$ppunit[rows]==hhmerge[rows,4+i]),"\n",sep="")
  cat("HH price > for product ",i,": ",sum(hhmerge$ppunit[rows]>hhmerge[rows,4+i]),"\n",sep="")
  cat("Store price > for product ",i,": ",sum(hhmerge$ppunit[rows]<hhmerge[rows,4+i]),"\n",sep="")
  cat("Differences: \n")
  print(summary(hhmerge$ppunit[rows]-hhmerge[rows,4+i]))
  hhmerge[rows,4+i] = hhmerge$ppunit[rows]
}

colnames(hhmerge)[ncol(hhmerge)-length(demovars)-1] = "totunits"

hhinit$BRAND = hhinit$L5

hhinit = hhinit[,c("PANID","WEEK","UNITS","DOLLARS","IRI_KEY","BRAND","VOL_EQ",demovars,"ppunit")]

hhinit$brindex = 0

for(i in 1:nrow(brsize)) {
  hhinit$brindex[hhinit$BRAND == names(brandtab)[brsize[i,1]] & hhinit$VOL_EQ == vols[brsize[i,2]]] = i
}

hhinit = hhinit[hhinit$PANID %in% unique(hhmerge$PANID),]

save(hhmerge,hhinit,brsize,vols,droplist,dropindx,hhinitmerge,prdf,irikeymerge,irikeyinitmerge,file=paste(outputpath,"mergedpanel_",product,"_allyears.RData",sep=""))



