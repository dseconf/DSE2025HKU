# begin data processing from make_ijc_data

# set negative prices to 1000 - this wasn't always being set properly

pricecols = which(substr(colnames(hhmerge),1,7)=="Product")
for(i in pricecols) {
  hhmerge[hhmerge[,i] < 0,i] = 1000
  hhinitmerge[hhinitmerge[,i] < 0,i] = 1000
}

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]
hhinitmerge = hhinitmerge[order(hhinitmerge$PANID,hhinitmerge$WEEK),]

# package size
packsize = vols[1:nsize]

# distribution of units by hh

unitsbyhh = aggregate(hhmerge$totunits,by=list(hhmerge$PANID),FUN=sum)

# number of trips per hh

obsbyhh = aggregate(hhmerge$totunits,by=list(hhmerge$PANID),FUN=length)

hhmerge$flag <- NULL

# remove individuals that don't appear in the initial data
# we could integrate simulate out choices for these guys, but will leave alone for now

hhinit = hhinit[order(hhinit$PANID,hhinit$WEEK),]

initid = unique(hhinit$PANID)

hhmerge = hhmerge[hhmerge$PANID %in% initid,]

hhinitmerge = hhinitmerge[hhinitmerge$PANID %in% initid,]

# construct household level consumption rate

unitsbyhh = aggregate(hhmerge$totunits,by=list(hhmerge$PANID),FUN=sum)
colnames(unitsbyhh) = c("PANID","hhunits")

hhmerge = merge(hhmerge,unitsbyhh)

brsize1 = rbind(c(0,0),brsize)

volbyhh = aggregate(c(0,packsize)[brsize1[hhmerge$brindex+1,2]+1]*hhmerge$totunits,by=list(hhmerge$PANID),FUN=sum)
colnames(volbyhh) = c("PANID","hhvol")

hhmerge = merge(hhmerge,volbyhh)

maxweekhh = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=max)
colnames(maxweekhh) = c("PANID","maxweekhh")
minweekhh = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=min)
colnames(minweekhh) = c("PANID","minweekhh")

mergemax = merge(hhmerge[,c("PANID","WEEK")],maxweekhh)
mergemax = mergemax[order(mergemax$PANID,mergemax$WEEK),]
mergemin = merge(hhmerge[,c("PANID","WEEK")],minweekhh)
mergemin = mergemin[order(mergemin$PANID,mergemin$WEEK),]

nweek = mergemax$maxweekhh-mergemin$minweekhh

hhmerge$crate = hhmerge$hhvol/nweek

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]

# construct interpurchase times

ip1 = data.frame(cbind(hhmerge$WEEK[hhmerge$totunits>0],hhmerge$PANID[hhmerge$totunits>0]))
colnames(ip1) = c("WEEK","PANID")
lpid1 = c(0,ip1$PANID[1:(nrow(ip1)-1)])
lweek1 = c(0,ip1$WEEK[1:(nrow(ip1)-1)])
iptime = data.frame(ip1$PANID[ip1$PANID==lpid1], lweek1[ip1$PANID==lpid1] - ip1$WEEK[ip1$PANID==lpid1])
colnames(iptime) = c("PANID","ip")
iptime$ip = -iptime$ip

maxweek = max(hhmerge$WEEK)
minweekp = rep(maxweek,nrow(hhmerge))

minweekp[hhmerge$totunits > 0] = hhmerge$WEEK[hhmerge$totunits > 0]
minweekphh = aggregate(minweekp,by=list(hhmerge$PANID),FUN=min)

colnames(minweekphh) = c("PANID","minweekphh")

maxweekp = rep(maxweek,nrow(hhmerge))

maxweekp[hhmerge$totunits > 0] = hhmerge$WEEK[hhmerge$totunits > 0]
maxweekphh = aggregate(maxweekp,by=list(hhmerge$PANID),FUN=max)

colnames(maxweekphh) = c("PANID","maxweekphh")

ipfirst = minweekphh$minweekphh - minweekhh$minweekhh
iplast =  maxweekphh$maxweekphh - maxweekhh$maxweekhh

minweek = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=min)

colnames(minweek) = c("PANID","minweek")

hhmerge = merge(hhmerge,minweek)

maxweek = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=max)

colnames(maxweek) = c("PANID","maxweek")

hhmerge = merge(hhmerge,maxweek)

# create imputed inventory measure - for constructing inventory state space

hhmerge$sumq = cumsum(c(0,packsize)[brsize1[hhmerge$brindex+1,2]+1]*hhmerge$totunits)
startq = hhmerge[hhmerge$WEEK==hhmerge$minweek,c("PANID","sumq")]
startq$startq = startq$sumq
startq$sumq <- NULL

hhmerge = merge(hhmerge,startq)

hhmerge$sumq = hhmerge$sumq - hhmerge$startq

hhmerge$inv = hhmerge$sumq - hhmerge$crate*(hhmerge$WEEK-hhmerge$minweek+1)

# scale up inventory by assuming in first purchase inventory is 0

hhmerge = merge(hhmerge,minweekphh)

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]

shifter = hhmerge[hhmerge$WEEK==hhmerge$minweekphh,c("PANID","totunits","inv","brindex")]
shifter$ishift = c(0,packsize)[brsize1[shifter$brindex+1,2]+1]*shifter$totunits - shifter$inv
shifter$totunits <- NULL
shifter$inv <- NULL
shifter$brindex <- NULL

hhmerge = merge(hhmerge,shifter)
hhmerge$invnew = hhmerge$inv + hhmerge$ishift

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]

# drop out observations where it looks like we see many or very few purchases

unitsperyr = 52*unitsbyhh$hhunits/(maxweek$maxweek-minweek$minweek)

qupryr = quantile(unitsperyr,probs=quprobs)

badid = maxweek$PANID[unitsperyr < qupryr[1] | unitsperyr > qupryr[2]]

update.droplist("hhmerge","number of units per year too high or too small",hhmerge,hhmerge[!(hhmerge$PANID %in% badid),],"PANID")

hhmerge = hhmerge[!(hhmerge$PANID %in% badid),]
nweek = nweek[!(hhmerge$PANID %in% badid)]

# drop observations where less than 5 units are ever purchased

badid2 = unique(hhmerge$PANID[hhmerge$hhunits < minunits])

maxip = aggregate(iptime$ip,by=list(iptime$PANID),FUN=max)

update.droplist("hhmerge","number of units below minunits",hhmerge,hhmerge[!(hhmerge$PANID %in% badid2),],"PANID")

hhmerge = hhmerge[!(hhmerge$PANID %in% badid2),]
nweek = nweek[!(hhmerge$PANID %in% badid2)]

# remove individuals who have a break of more than 40 weeks

badid3 = maxip[maxip$x > maxpgap,1]

update.droplist("hhmerge","purchase gap too large",hhmerge,hhmerge[!(hhmerge$PANID %in% badid3),],"PANID")

hhmerge = hhmerge[!(hhmerge$PANID %in% badid3),]
nweek = nweek[!(hhmerge$PANID %in% badid3)]

# remove individuals who purchase more than bigJ units

badid4 = unique(hhmerge$PANID[hhmerge$totunits > bigJ])

update.droplist("hhmerge","totunits above bigJ",hhmerge,hhmerge[!(hhmerge$PANID %in% badid4),],"PANID")

hhmerge = hhmerge[!(hhmerge$PANID %in% badid4),]
nweek = nweek[!(hhmerge$PANID %in% badid4)]

# remove individuals with consumption rates that are too high or too low

obskeep = hhmerge$crate > cratebounds[1] & hhmerge$crate < cratebounds[2]

update.droplist("hhmerge","consumption rates outside of upper or lower bounds",hhmerge,hhmerge[obskeep,],"PANID")

hhmerge = hhmerge[obskeep,]

# some products are never bought by consumers in the estimation sample - need to drop them

goodbrand = which(1:nrow(brsize) %in% as.integer(names(table(hhmerge$brindex))))

badbrand = setdiff(1:nrow(brsize),goodbrand)

xx = data.frame(cbind(1:length(goodbrand),goodbrand))
colnames(xx) = c("newind","oldind")

xy = data.frame(1:max(goodbrand))
colnames(xy) = c("oldind")

xx = merge(xy,xx,all.x=TRUE)

hhmerge$brindex[hhmerge$brindex > 0] = xx$newind[hhmerge$brindex[hhmerge$brindex > 0]]

# remove products that are not in goodbrand

for(i in 1:nrow(brsize)) {
  if(i %in% badbrand) {
    pname = paste("Product",brsize[i,1],"_",brsize[i,2],sep="")
    hhmerge[[pname]] <- NULL

    pname = paste("Feat",brsize[i,1],"_",brsize[i,2],sep="")
    hhmerge[[pname]] <- NULL

    pname = paste("Disp",brsize[i,1],"_",brsize[i,2],sep="")
    hhmerge[[pname]] <- NULL

    pname = paste("PR",brsize[i,1],"_",brsize[i,2],sep="")
    hhmerge[[pname]] <- NULL
  }
}

brsize = brsize[goodbrand,]

# fill in missing prices with modal prices
# I have found that there is usually very little missing on this end

if(fillprices) {
  for(i in 1:nrow(brsize)) {
    pname = paste("Product",brsize[i,1],"_",brsize[i,2],sep="")
    ptable = table(hhmerge[[pname]])
    if(sum(hhmerge[[pname]]==1000)>0) {
      hhmerge[[pname]][hhmerge[[pname]]==1000] = as.numeric(names(ptable)[which(ptable==min(max(ptable)))])
    }
  }
}

# reduce demographic variables to categories we can feed to the estimation
# code

hhmerge = hhmerge[order(hhmerge$PANID,hhmerge$WEEK),]
lpanid = c(0,hhmerge$PANID[1:(nrow(hhmerge)-1)])

first = hhmerge$PANID != lpanid

# compute percentiles of income, age

if("HH_INC" %in% demovars) {
  medincome = quantile(hhmerge$HH_INC[first],probs=0.5)
  hhmerge$highincome = hhmerge$HH_INC >= medincome
}

if("HH_AGE" %in% demovars) {
  medage = quantile(hhmerge$HH_AGE[first],probs=0.5)  # median age is 55+
  hhmerge$agedummy = hhmerge$HH_AGE >= medage
}

if("HH_EDU" %in% demovars) {
  somecollege = hhmerge$HH_EDU%in% c(5,6)
  college = hhmerge$HH_EDU%in% c(7,8)

  hhmerge$somecollege = somecollege
  hhmerge$college = college
}

if("HH_SIZE" %in% demovars) {
  fsize2 = hhmerge$HH_SIZE > 2
  hhmerge$fsize2 = fsize2
}

if("child" %in% demovars) {
  fchild = hhmerge$child < 8
  hhmerge$fchild = fchild
}

# data for estimating brand coefficients

hhmergebr = hhmerge[hhmerge$totunits > 0,]

# data that includes no purchase and no store visit weeks

allid = unique(hhmerge$PANID)
allid = allid[order(allid)]

allweek = unique(hhmerge$WEEK)
allweek = allweek[order(allweek)]

hhbig = data.frame(cbind(kronecker(allid,rep(1,length(allweek))),rep(allweek,length(allid))))
colnames(hhbig) = c("PANID","WEEK")

hhbig = merge(hhbig,hhmerge,all.x=TRUE)

hhbig$brindex[is.na(hhbig$brindex)] = -1
hhbig$totunits[is.na(hhbig$totunits)] = -1
hhbig = hhbig[order(hhbig$PANID,hhbig$WEEK),]

# there are some individuals who drop out of the dataset early, so remove those observations

weekbds = list(c(1114,1165),c(1166,1217),c(1218,1269),c(1270,1321),c(1322,1373),c(1374,1426),c(1427,1478))

maxweekm = aggregate(hhmerge$WEEK,by=list(hhmerge$PANID),FUN=max)
colnames(maxweekm) = c("PANID","maxweek")
hhbig$maxweek <- NULL
hhbig = merge(hhbig,maxweekm)
hhbig = hhbig[order(hhbig$PANID,hhbig$WEEK),]
for(i in 1:7) {
  hhbig$maxweek[hhbig$maxweek >= weekbds[[i]][1] & hhbig$maxweek <= weekbds[[i]][2]] = weekbds[[i]][2]
}

hhbig = hhbig[hhbig$WEEK <= hhbig$maxweek,]
hhbig$maxweek <- NULL

# resize hhinit to match the new data

nobs = nrow(hhbig)
nhhs = length(unique(hhbig$PANID))

bigid = unique(hhbig$PANID)

hhinit = hhinit[hhinit$PANID %in% bigid,]
initid = unique(hhinit$PANID)

hhinitmerge = hhinitmerge[hhinitmerge$PANID %in% bigid,]

if(!setequal(bigid,initid)) {
  stop("different ids in init data and estimation panel")
}

for(j in 1:nsize) {
  hhinit[[paste("Volunits",j,sep="")]] = 0
}
hhinit$VolunitsOther = 0

for(j in 1:nsize) {
  hhinit[[paste("Volunits",j,sep="")]][hhinit$VOL_EQ == vols[j]] = hhinit$UNITS[hhinit$VOL_EQ == vols[j]]
}
hhinit$VolunitsOther[!(hhinit$VOL_EQ %in% vols[1:nsize])] = hhinit$UNITS[!(hhinit$VOL_EQ %in% vols[1:nsize])]

# make hhinitbig
#allidinit = unique(hhminiterge$PANID)
#allidinit = allidinit[order(allidinit)]

#allweekinit = unique(hhinitmerge$WEEK)
#allweekinit = allweekinit[order(allweekinit)]

#hhinitbig = data.frame(cbind(kronecker(allidinit,rep(1,length(allweekinit))),rep(allweekinit,length(allidinit))))
#colnames(hhinitbig) = c("PANID","WEEK")

#hhinitbig = merge(hhinitbig,hhinitmerge,all.x=TRUE)

#hhinitbig$brindex[is.na(hhinitbig$brindex)] = -1
#hhinitbig$totunits[is.na(hhinitbig$totunits)] = -1
#hhinitbig = hhinitbig[order(hhinitbig$PANID,hhinitbig$WEEK),]


# aggregate hhinit to remove duplicate observations
initag = aggregate(hhinit$VOL_EQ*hhinit$UNITS,by=list(hhinit$PANID,hhinit$WEEK),FUN=sum)
colnames(initag) = c("PANID","WEEK","UNITS")

for(j in 1:nsize) {
  initag1 = aggregate(hhinit[[paste("Volunits",j,sep="")]],by=list(hhinit$PANID,hhinit$WEEK),FUN=sum)
  colnames(initag1) = c("PANID","WEEK",paste("Volunits",j,sep=""))
  initag = cbind(initag,initag1[[paste("Volunits",j,sep="")]])
  colnames(initag)[length(colnames(initag))] = paste("Volunits",j,sep="")
}

initag1 = aggregate(hhinit$VolunitsOther,by=list(hhinit$PANID,hhinit$WEEK),FUN=sum)
colnames(initag1) = c("PANID","WEEK","VolunitsOther")
initag = cbind(initag,initag1$VolunitsOther)
colnames(initag)[length(colnames(initag))] = "VolunitsOther"

iweek = min(hhinit$WEEK):max(hhinit$WEEK)
initnew = data.frame(cbind(kronecker(initid,rep(1,length(iweek))),rep(iweek,length(initid))))
colnames(initnew) = c("PANID","WEEK")

hhinitbig = merge(initnew,initag,all.x = TRUE)

hhinitbig = hhinitbig[order(hhinitbig$PANID,hhinitbig$WEEK),]
hhinitbig$UNITS[is.na(hhinitbig$UNITS)] = 0

for(j in 1:nsize) {
  hhinitbig[[paste("Volunits",j,sep="")]][is.na(hhinitbig[[paste("Volunits",j,sep="")]])] = 0
}
hhinitbig$VolunitsOther[is.na(hhinitbig$VolunitsOther)] = 0

# merge in prices - used in simulation from the data

hhinitbig = merge(hhinitbig,hhinitmerge[,c("PANID","WEEK","brindex","ppunit",paste("Product",brsize[,1],"_",brsize[,2],sep=""))],all.x=TRUE)

hhinitbig = hhinitbig[,c("PANID","WEEK","brindex","ppunit",paste("Product",brsize[,1],"_",brsize[,2],sep=""),"UNITS",paste("Volunits",1:nsize,sep=""),"VolunitsOther")]

hhinitbig = hhinitbig[order(hhinitbig$PANID,hhinitbig$WEEK),]

hhinitbig$brindex[is.na(hhinitbig$brindex)] = -1

save(hhmerge,hhmergebr,hhbig,hhinitmerge,hhinitbig,brsize,vols,droplist,dropindx,bigJ,file=paste(outputpath,"estimationdata_",product,"_allyears.RData",sep=""))
