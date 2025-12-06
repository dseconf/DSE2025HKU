# functions that are used by estimation_starter.r
# and artificialdata_ijc.r

# compute purchase hazard
# this is meant to work with hhbig where we aren't missing
# intervening observations

purch.hazard <- function(totunits,data,maxip = 40) {

  #ip1 = data[,c("WEEK","PANID")]
  phazard = rep(0,maxip)
  for(i in 1:maxip) {
    lpid1 = c(rep(0,i),data$PANID[1:(nrow(data)-i)])
    lpurch1 = c(rep(0,i),totunits[1:(nrow(data)-i)])
    s1 = sum( (totunits > 0)*(lpid1 == data$PANID)*(lpurch1 > 0) )
    s2 = sum( (lpurch1 > 0)*(lpid1 == data$PANID) )
    phazard[i] = s1/s2
  }

  return(phazard)

}

# create data on interpurchase times

do.iptime <- function(totunits,data) {

  ip1 = data.frame(cbind(data$WEEK[totunits>0],data$PANID[totunits>0]))
  colnames(ip1) = c("WEEK","PANID")
  lpid1 = c(0,ip1$PANID[1:(nrow(ip1)-1)])
  lweek1 = c(0,ip1$WEEK[1:(nrow(ip1)-1)])
  iptime = data.frame(ip1$PANID[ip1$PANID==lpid1], lweek1[ip1$PANID==lpid1] - ip1$WEEK[ip1$PANID==lpid1])
  colnames(iptime) = c("PANID","ip")
  iptime$ip = -iptime$ip

  return(iptime)

}

# transform a small vector of parameters into a big one with bounds imposed

tform <- function(par,info) {

  np = length(info$xfull)
  nps = nrow(par)

  tf = info$tform

  parbig = matrix(0,nrow=np,ncol=info$nhhs)

  j=0
  for(i in 1:np) {
    if(!info$fixed[i]) {
      j=j+1
      parbig[i,] = par[j,]
      if(abs(tf[i]) == 1) {
        parbig[i,] = sign(tf[i])*exp(parbig[i,])
      } else if (abs(tf[i]) == 2) {
        parbig[i,] = info$lbounds[i] + (info$ubounds[i]-info$lbounds[i])*exp(parbig[i,])/(1+exp(parbig[i,]))
      }
    } else {
      if(info$paramequal[i] > 0) {
        if(info$paramequal[i] >= i) {
          stop("Error - paramequal[i] must be less than i")
        }
        parbig[i,] = parbig[info$paramequal[i],]
      } else {
        if(!is.null(info$paramstart)) {
          parbig[i,] = info$paramstart[i,]
        } else {
          parbig[i,] = info$xfull[i]
        }
      }
    }
  }

  if(info$crfix) {
    parbig[info$nbrand+2,] = parbig[info$nbrand+1,]
  }

  if(info$data.crate) {
    parbig[info$nbrand+1,] = info$crate
    parbig[info$nbrand+2,] = info$crate
  }

  if(info$sizeshifter) {
    for(i in 1:info$nbrand) {
      if(info$brsize[i,2] > 1 & info$sizebrand[i]) {
        parbig[i,] = parbig[i,] + parbig[npbig-info$nsize+info$brsize[i,2],]
      }
    }
  }

  return(parbig)

}

# undo transformation - take bigger matrix & turn to small one

untform <- function(par,info) {

  np = length(info$xfull)
  nps = sum(!info$fixed)

  tf = info$tform

  parsmall = matrix(0,nrow=nps,ncol=info$nhhs)

  j=0
  for(i in 1:np) {
    if(!info$fixed[i]) {
      j=j+1
      parsmall[j,] = par[i,]
      if(abs(tf[i]) == 1) {
        parsmall[j,] = log(sign(tf[i])*par[i,])
      } else if (abs(tf[i]) == 2) {
        a = info$lbounds[i]
        h = info$ubounds[i] - info$lbounds[i]
        parsmall[j,] = log((par[i,]-a)/h) - log(1-(par[i,]-a)/h)
      }
      if(info$sizeshifter) {
        if(i <= info$nbrand & info$sizebrand[i]) {
          if(info$brsize[i,2] > 1) {
            parsmall[j,] = parsmall[j,] - par[np-info$nsize+info$brsize[i,2],]
          }
        }
      }
    }
  }

  return(parsmall)

}

# compute bstate matrix

getbstates <- function(nsize,maxbottles) {

  nbstates = nsize+1
  for(i in 2:maxbottles) {
    nbstates = nbstates + nsize^i
  }

  bstates = matrix(as.integer(0),nrow=maxbottles,ncol=nbstates)
  for(i in 1:nsize) {
    bstates[1,i+1] = i
  }
  # recursion would be better here
  if(maxbottles > 1) {
    sti=nsize+1
    for(i in 1:nsize) {
      for(j in 1:nsize) {
        sti = sti + 1
        bstates[1,sti] = i
        bstates[2,sti] = j
      }
    }
    if(maxbottles > 2) {
      for(i in 1:nsize) {
        for(j in 1:nsize) {
          for(k in 1:nsize) {
            sti = sti + 1
            bstates[1,sti] = i
            bstates[2,sti] = j
            bstates[3,sti] = k
          }
        }
      }
    }
    if(maxbottles > 3) {
      for(i in 1:nsize) {
        for(j in 1:nsize) {
          for(k in 1:nsize) {
            for(l in 1:nsize) {
              sti = sti + 1
              bstates[1,sti] = i
              bstates[2,sti] = j
              bstates[3,sti] = k
              bstates[4,sti] = l
            }
          }
        }
      }
    }
    if(maxbottles > 4) {
      for(i in 1:nsize) {
        for(j in 1:nsize) {
          for(k in 1:nsize) {
            for(l in 1:nsize) {
              for(m in 1:nsize) {
                sti = sti + 1
                bstates[1,sti] = i
                bstates[2,sti] = j
                bstates[3,sti] = k
                bstates[4,sti] = l
                bstates[5,sti] = m
              }
            }
          }
        }
      }
    }
    if(maxbottles > 5) {
      for(i in 1:nsize) {
        for(j in 1:nsize) {
          for(k in 1:nsize) {
            for(l in 1:nsize) {
              for(m in 1:nsize) {
                for(n in 1:nsize) {
                  sti = sti + 1
                  bstates[1,sti] = i
                  bstates[2,sti] = j
                  bstates[3,sti] = k
                  bstates[4,sti] = l
                  bstates[5,sti] = m
                  bstates[6,sti] = n
                }
              }
            }
          }
        }
      }
    }
    if(maxbottles > 6) {
      for(i in 1:nsize) {
        for(j in 1:nsize) {
          for(k in 1:nsize) {
            for(l in 1:nsize) {
              for(m in 1:nsize) {
                for(n in 1:nsize) {
                  for(o in 1:nsize) {
                    sti = sti + 1
                    bstates[1,sti] = i
                    bstates[2,sti] = j
                    bstates[3,sti] = k
                    bstates[4,sti] = l
                    bstates[5,sti] = m
                    bstates[6,sti] = n
                    bstates[7,sti] = o
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  #make an array mapping states to indexes.  Needed for transitions
  revarray = array(as.integer(-1),rep(nsize+1,maxbottles))
  #browser()
  for(i in 1:nbstates) {
    inds = matrix(bstates[,i]+1,nrow=1,ncol=maxbottles,byrow=TRUE)
    revarray[inds]=i
  }
  return(list(as.integer(nbstates),bstates,revarray))

}

# brand coefficient only log-likelihood

llbrand <- function(x,Mcmc,info,hhmergebr,printing=FALSE) {

  if(length(x) != sum(!info$fixed)) {
    stop("parameter length doesn't match number of nonfixed params")
  }
  inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
  inputs$pstart[!info$fixed,] = x

  temp = tform(inputs$pstart[!info$fixed,],info)

  if(printing) {
    cat("Input parameters: ","\n",sep="")
    print(temp[!info$fixed,1])
  }

  ptm = proc.time()[3]
  br.out = brandlik(inputs,Mcmc,info,hhmergebr)

  if(printing) {
    cat("elapsed time:", proc.time()[3]-ptm,"\n")

    cat("log likelihood: ",sum(br.out$llbrand),"\n")

  }

  #return(br.out)
  return(-sum(br.out$llbrand))

}


# compute log-likelihood - for maximum likelihood estimation
# this version solves for the value function

ll <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,indflag=FALSE,hhinclude=1:nhhs) {

  if(length(x) != sum(!info$fixed)) {
    stop("parameter length doesn't match number of nonfixed params")
  }
  if(!is.null(info$paramstart)) {
    stop("paramstart exists - mle will not take inputs")
  }
  inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
  inputs$pstart[!info$fixed,] = x

  temp = tform(matrix(inputs$pstart[!info$fixed,],nrow=sum(!info$fixed),ncol=nhhs),info)
  temp = matrix(temp[!info$fixed,1],nrow=sum(!info$fixed),ncol=1)

  cat("Input parameters: ","\n",sep="")
  if(!is.null(Mcmc$pnames)) {
    rownames(temp) = Mcmc$pnames[!info$fixed]
  }
  print(temp)

  dfind = which(rownames(temp) == "discount factor")
  scind = which(rownames(temp) == "stockout")
  pind = which(rownames(temp) == "price")
  baddf = FALSE
  badsc = FALSE
  badp = FALSE

  if(length(dfind) > 0) {
    baddf = temp[dfind] > 0.999
  }

  if(length(scind) > 0) {
    badsc = temp[scind] > 25
  }

  if(length(pind) > 0) {
    badp = abs(temp[pind]) < 1e-5
  }

  if(max(abs(temp)) < 500 & !baddf & !badsc & !badp & max(abs(temp[1:pind])) < 25) {
    ptm = proc.time()[3]
    ml.out = mlelik(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr)
    cat("elapsed time:", proc.time()[3]-ptm,"\n")

    cat("log likelihood: ",sum(ml.out$llbrand+ml.out$llqhh),"\n")

    if(indflag) {
      return(ml.out)
    } else {
      return(-sum(ml.out$llbrand[hhinclude]+ml.out$llqhh[hhinclude]))
    }

  } else {

    return(1e6)

  }

}

# non transformed for getting hessian

ll.hess <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,indflag=FALSE,hhinclude=1:nhhs) {

  if(length(x) != sum(!info$fixed)) {
    stop("parameter length doesn't match number of nonfixed params")
  }
  if(!is.null(info$paramstart)) {
    stop("paramstart exists - mle will not take inputs")
  }
  inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
  inputs$pstart[!info$fixed,] = x

  temp = tform(matrix(inputs$pstart[!info$fixed,],nrow=sum(!info$fixed),ncol=nhhs),info)

  x1 = temp[!info$fixed,1]

  info1 = info
  info1$tform[!info$fixed] = 0

  ll1 <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,hhinclude) {
    inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
    inputs$pstart[!info$fixed,] = x
    ml.out = mlelik(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr)
    return(-sum(ml.out$llbrand[hhinclude]+ml.out$llqhh[hhinclude]))
  }

  cat("Input parameters: ","\n",sep="")
  print(temp[!info$fixed,1])

  return(hessian(ll1,x1,Mcmc=Mcmc,info=info1,hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr,hhinclude=hhinclude))

}


# compute log likelihood variance matrix

llcov <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,gradtol = 1e-4,hhinclude=1:nhhs) {

  if(length(x) != sum(!info$fixed)) {
    stop("parameter length doesn't match number of nonfixed params")
  }

  inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
  inputs$pstart[!info$fixed,] = x

  #temp = tform(inputs$pstart[!info$fixed,],info)
  temp = tform(matrix(inputs$pstart[!info$fixed,],nrow=sum(!info$fixed),ncol=nhhs),info)

  cat("Input parameters: ","\n",sep="")
  print(temp[!info$fixed,1])

  ml.out = mlelik(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr)
  #cat("abc\n")
  ll1 = ml.out$llbrand+ml.out$llqhh

  llgrad = matrix(0,nrow=length(x),ncol=length(hhinclude))

  for(i in 1:length(x)) {
    #cat("i ",i,"\n")
    x1 = x
    if(x1[i] == 0) {
      h=gradtol
    } else {
      h = gradtol*abs(x1[i])
    }
    x1[i] = x1[i] + h

    inputs = list(pstart = matrix(info$xfull,nrow=length(info$xfull),ncol=nhhs))
    inputs$pstart[!info$fixed,] = x1

    #temp = tform(inputs$pstart[!info$fixed,],info)
    temp = tform(matrix(inputs$pstart[!info$fixed,],nrow=sum(!info$fixed),ncol=nhhs),info)

    o1 = mlelik(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr)

    llgrad[i,] = (o1$llbrand[hhinclude]+o1$llqhh[hhinclude] - ll1[hhinclude])/h

  }

  cat("information matrix rank ", rankMatrix(llgrad%*%t(llgrad)), ", (", length(x),")\n", sep="")

  if(rankMatrix(llgrad%*%t(llgrad)) < length(x)) {

    vcov = ginv(llgrad%*%t(llgrad))

    iflag = FALSE

  } else {

    vcov = chol2inv(chol(llgrad%*%t(llgrad)))

    iflag = TRUE

  }

  return(list(vcov,llgrad,iflag))

}


# function to create simulated data from parameter estimates

simulate.data <- function(xin,hhbig,hhmerge,hhmergebr,info,Mcmc,loadvf=FALSE,
                          vffile = "vfest_save.RData",seed=4321581) {

  set.seed(seed)

  nhhs = info$nhhs

  param = tform(matrix(rep(xin,nhhs),nrow=length(xin),ncol=nhhs),info)

  paramstart = untform(param,info)

  bwmatrix = diag(nrow(param))

  if(loadvf) {
    load(vffile)
  } else {
    out.vf = computevf(paramstart,1+249*(!myopic),info,hhbig,hhmerge,hhmergebr,2,as.integer(gridmethod),as.integer(1),as.integer(1),bwmatrix)
    save(out.vf,file=vffile)
    cat("vf computed\n")
  }

  data2 = hhmerge[,c("PANID","WEEK",colnames(hhmerge)[substr(colnames(hhmerge),1,7)=="Product"])]

  data1 = merge(hhbig[,1:(info$prindex-1)],data2,all.x = TRUE)

  data1$totunits = hhbig$totunits

  nchh = ncol(data1)

  rm(data2)

  done = FALSE

  iter = 0

  out.vf$ivbig = out.vf$ivbig[info$expandbig,]

  bwmatrix = diag(nrow(param))/( 4/(3*Mcmc$nsave) )^0.2

  indexes = 1:10

  simchoices = choicesimijc(param, data1, info, out.vf$ivinfo, out.vf$ivbig,
  out.vf$vfbig, out.vf$oldparam, out.vf$rgridsave,as.integer(info$nrgrid),
  as.integer(indexes), as.integer(c(1,nhhs)), as.integer(info$prindex),
  as.integer(gridmethod),Mcmc$vfinterpmethod,bwmatrix)

  return(list(out.vf = out.vf, sim = simchoices))

}

# compute inventory, for sample selection

getinv <- function(hhinitbig,hhbig,hhcrate,brsize,vols) {

  ninitt = length(unique(hhinitbig$WEEK))

  allids = unique(hhinitbig$PANID)

  allids = allids[order(allids)]

  nt = length(unique(hhbig$WEEK))

  allweeks = unique(hhbig$WEEK)
  allweeks = allweeks[order(allweeks)]

  initx = matrix(hhinitbig$UNITS,nrow=ninitt,ncol=length(allids),byrow=FALSE)  #this is volume, so we are ok

  inds = hhbig$brindex > 0

  volsbig = data.frame(hhbig$PANID[inds],hhbig$WEEK[inds],vols[brsize[hhbig$brindex[inds],2]])
  colnames(volsbig) = c("PANID","WEEK","volume")

  hhbig = merge(hhbig,volsbig,all.x=TRUE)
  hhbig = hhbig[order(hhbig$PANID,hhbig$WEEK),]
  hhbig$volume[is.na(hhbig$volume)] = 0
  hhbig$volume = hhbig$volume*hhbig$totunits

  hhbig2 = data.frame(kronecker(allids,rep(1,nt)),rep(allweeks,length(allids)))
  colnames(hhbig2) = c("PANID","WEEK")

  hhbig2 = merge(hhbig2,hhbig[,c("PANID","WEEK","volume")],all.x=TRUE)

  hhbig2 = hhbig2[order(hhbig2$PANID,hhbig2$WEEK),]

  xbig = matrix(hhbig2$volume,nrow=nt,ncol=length(allids),byrow=FALSE)

  #initinv = NULL
  #inv = rep(0,length(allids))
  for(i in 1:ninitt) {
    if(i==1) {
      initinv = matrix(pmax(initx[i,]-hhcrate,0),nrow=1,ncol=length(allids))
    } else {
      initinv = rbind(initinv,matrix(pmax(inv+initx[i,]-hhcrate,0),nrow=1,ncol=length(allids)))
    }
    inv = initinv[i,]
  }

  for(i in 1:nt) {
    if(i==1) {
      invbig = matrix(pmax(initinv[ninitt,]+xbig[i,]-hhcrate,0),nrow=1,ncol=length(allids))
      nzero = matrix(as.numeric(invbig[i,]==0),nrow=1,ncol=length(allids))
    } else {
      invbig = rbind(invbig,matrix(pmax(inv+xbig[i,]-hhcrate,0),nrow=1,ncol=length(allids)))
      nzero = rbind(nzero,matrix((invbig[i,]==0)*(1+nzero[i-1,]),nrow=1,ncol=length(allids)))
    }
    inv = invbig[i,]
  }

  return(list(initinv=initinv,invbig=invbig,nzero=nzero,xbig=xbig,initx=initx))

}

# functions for making adjustments to inventory/consumption rates

cppFunction('double minhhinv(NumericVector crate, NumericVector initx, NumericVector xbig) {
  int inobs = initx.size();
  int nobs = xbig.size();
  double cr = crate(0);
  double inv = 0, mininv = 0;
  for(int i=0;i<inobs;i++) {
    inv += initx(i) - cr;
    if(inv < 0) {
      inv = 0;
    }
  }
  for(int i=0;i<nobs;i++) {
    inv += xbig(i) - cr;
    if(inv < 0) {
      inv = 0;
    }
    if(i==0) {
      mininv = inv;
    } else {
      if(inv < mininv) {
        mininv = inv;
      }
    }
  }
  return(mininv);
}')

getnewcrate <- function(crate,initx,xbig) {
  nhhs = ncol(initx)
  crnew = rep(0,nhhs)
  for(i in 1:nhhs) {
    if(minhhinv(hhcrate[i],inv1$initx[,i],inv1$xbig[,i]) > 0) {
      delt = 0.01
      crstart = crate[i]
      while(minhhinv(crstart,inv1$initx[,i],inv1$xbig[,i]) > 0) {
        crstart = crstart + delt
      }
      crnew[i] = crstart
    } else {
      crnew[i] = crate[i]
    }
  }
  return(crnew)
}

cppFunction('NumericVector avghhinv(NumericVector crate, NumericVector initx, NumericVector xbig) {
  int inobs = initx.size();
  int nobs = xbig.size();
  double cr = crate(0);
  double inv = 0, mininv = 0;
  double iavgi = 0, avgi = 0;
  for(int i=0;i<inobs;i++) {
    inv += initx(i) - cr;
    if(inv < 0) {
      inv = 0;
    }
    iavgi += inv;
  }
  iavgi /= (double)inobs;
  for(int i=0;i<nobs;i++) {
    inv += xbig(i) - cr;
    if(inv < 0) {
      inv = 0;
    }
    if(i==0) {
      mininv = inv;
    } else {
      if(inv < mininv) {
        mininv = inv;
      }
    }
    avgi += inv;
  }
  avgi /= (double)nobs;
  NumericVector ret(2);
  ret[0] = iavgi;
  ret[1] = avgi;
  return(ret);
}')

getnewcrate2 <- function(initx,xbig,crub) {
  foo <-  function(x,initx,xbig) {
    temp = avghhinv(x,initx,xbig)
    return(abs(temp[1]-temp[2]))
  }
  nhhs = ncol(initx)
  crnew = rep(0,nhhs)
  crbdseq = seq(0.2,crub,0.2)
  for(i in 1:nhhs) {
    done = FALSE
    j = 1
    while(!done) {
      crnew[i] = optimize(foo,c(0,crbdseq[j]),initx=initx[,i],xbig=xbig[!is.na(xbig[,i]),i])$minimum
      j=j+1
      done = j > length(crbdseq) | abs(crnew[i] - crbdseq[j-1]) > 1e-4
      #if(i==2) {browser()}
    }
  }
  return(crnew)
}

# compute empirical purchase probs given imputed inventory

emp.probs <- function(hhbig,hhmerge,hhinitbig,hhcrate,brsize,vols,first,
                      ibds = seq(0,50,0.5)) {

  allids = unique(hhinitbig$PANID)

  allids = allids[order(allids)]

  inv1 = getinv(hhinitbig,hhbig,hhcrate,brsize,vols)

  allweek = unique(hhbig$WEEK)
  allweek = allweek[order(allweek)]
  nt = length(allweek)

  invdata = data.frame(rep(allids[1],nt),allweek,inv1$invbig[,1])
  colnames(invdata) = c("PANID","WEEK","invc")

  for(i in 2:length(allids)) {
    i1 = data.frame(rep(allids[i],nt),allweek,inv1$invbig[,i])
    colnames(i1) = c("PANID","WEEK","invc")
    invdata = rbind(invdata,i1)
  }
  rm(i1)

  hhm1 = hhmerge[,c("PANID","WEEK","brindex")]

  hhm1 = merge(hhm1,invdata)

  hhm1 = hhm1[order(hhm1$PANID,hhm1$WEEK),]

  hhm1$linv = c(NA,hhm1$invc[1:(nrow(hhm1)-1)])

  hhm1 = hhm1[!first,]

  hhm1$madepurchase = hhm1$brindex > 0

  pprobs1 = rep(0,length(ibds))

  pprobs1[1] = sum(hhm1$madepurchase[hhm1$linv==0])/sum(hhm1$linv==0)

  for(i in 2:length(ibds)) {
    inds = hhm1$linv > ibds[i-1] & hhm1$linv <= ibds[i]
    pprobs1[i] = sum(hhm1$madepurchase[inds])/sum(inds)
  }

  return(list(pprobs=pprobs1,hhm=hhm1))

}

# predicted probs as a fn of parameter

pred.probs <- function(xin,Mcmc,info,hhbig,hhmerge,hhmergebr,
                       hhm,first,ibds = seq(0,50,0.5),all=FALSE) {

  temp = ll(xin,Mcmc=Mcmc,info=info,
            hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr,
            indflag=TRUE)

  predprobs = data.frame(cbind(hhm$linv,1-temp$llbig[hhbig$brindex>=0,1,1][!first]))
  colnames(predprobs) = c("linv","predprob")

  pprobs2 = rep(0,length(ibds))

  pprobs2[1] = sum(predprobs$predprob[hhm$linv==0])/sum(hhm$linv==0)

  for(i in 2:length(ibds)) {
    inds = hhm$linv > ibds[i-1] & hhm$linv <= ibds[i]
    pprobs2[i] = sum(predprobs$predprob[inds])/sum(inds)
  }

  if(all) {
    return(list(pprobs2,temp))
  } else {
    return(pprobs2)
  }

}

pprob.fit <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,
                      hhm,first,e.pprobs,ibds,xin1,xf) {

  xin = xin1
  xin[!xf] = x

  pprobs = pred.probs(xin,Mcmc,info,hhbig,hhmerge,hhmergebr,
                      hhm,first,ibds = ibds)
  #browser()
  return(sum( (pprobs-e.pprobs)*(pprobs-e.pprobs), na.rm=TRUE ))

}


imputeinv <- function(hhbig,hhmerge,hhinitbig,hhcrate,brsize,vols) {

  allids = unique(hhinitbig$PANID)

  allids = allids[order(allids)]

  inv1 = getinv(hhinitbig,hhbig,hhcrate,brsize,vols)

  allweek = unique(hhbig$WEEK)
  allweek = allweek[order(allweek)]
  nt = length(allweek)

  allweeki = unique(hhinitbig$WEEK)
  allweeki = allweeki[order(allweeki)]
  nti = length(allweeki)

  invdatai = data.frame(rep(allids[1],nti),allweeki,inv1$initinv[,1])
  colnames(invdatai) = c("PANID","WEEK","invc")

  for(i in 2:length(allids)) {
    i1 = data.frame(rep(allids[i],nti),allweeki,inv1$initinv[,i])
    colnames(i1) = c("PANID","WEEK","invc")
    invdatai = rbind(invdatai,i1)
  }
  rm(i1)

  invdata = data.frame(rep(allids[1],nt),allweek,inv1$invbig[,1])
  colnames(invdata) = c("PANID","WEEK","invc")

  for(i in 2:length(allids)) {
    i1 = data.frame(rep(allids[i],nt),allweek,inv1$invbig[,i])
    colnames(i1) = c("PANID","WEEK","invc")
    invdata = rbind(invdata,i1)
  }
  rm(i1)

  hhbi = hhinitbig[,c("PANID","WEEK","brindex")]

  hhbi = merge(hhbi,invdatai)

  hhbi = hhbi[order(hhbi$PANID,hhbi$WEEK),]

  hhb1 = hhbig[,c("PANID","WEEK","brindex")]

  hhb1 = merge(hhb1,invdata)

  hhb1 = hhb1[order(hhb1$PANID,hhb1$WEEK),]

  return(list(hhbi,hhb1))

}


# compute log-likelihood - for individual-specific parameters
# this is not set up for optimization

ll.heterogeneous <- function(x,Mcmc,info,hhbig,hhmerge,hhmergebr,indflag=FALSE,hhinclude=1:nhhs) {

  inputs = list(pstart = x)

  ptm = proc.time()[3]
  ml.out = mlelik(inputs,Mcmc,info,hhbig,hhmerge,hhmergebr)
  cat("elapsed time:", proc.time()[3]-ptm,"\n")

  cat("log likelihood: ",sum(ml.out$llbrand+ml.out$llqhh),"\n")

  if(indflag) {
    return(ml.out)
  } else {
    return(-sum(ml.out$llbrand[hhinclude]+ml.out$llqhh[hhinclude]))
  }

}

# pred probs on heterogeneous estimates

pred.probs.heterogeneous <- function(xin,Mcmc,info,hhbig,hhmerge,hhmergebr,
                       hhm,first,ibds = seq(0,50,0.5),all=FALSE) {

  temp = ll.heterogeneous(xin,Mcmc=Mcmc,info=info,
            hhbig=hhbig,hhmerge=hhmerge,hhmergebr=hhmergebr,
            indflag=TRUE)
  #browser()
  predprobs = data.frame(cbind(hhm$linv,1-temp$llbig[hhbig$brindex>=0,1,1][!first]))
  colnames(predprobs) = c("linv","predprob")

  pprobs2 = rep(0,length(ibds))

  pprobs2[1] = sum(predprobs$predprob[hhm$linv==0])/sum(hhm$linv==0)

  for(i in 2:length(ibds)) {
    inds = hhm$linv > ibds[i-1] & hhm$linv <= ibds[i]
    pprobs2[i] = sum(predprobs$predprob[inds])/sum(inds)
  }

  if(all) {
    return(list(pprobs2,temp))
  } else {
    return(pprobs2)
  }

}
