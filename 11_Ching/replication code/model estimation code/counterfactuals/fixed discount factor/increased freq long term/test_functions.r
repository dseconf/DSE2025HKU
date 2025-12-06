# functions to test the IJC code

# spline functions


# b spline phi

phi.b <- function(t) {
  at = abs(t)
  return(as.numeric(at > 1 & at <= 2)*(2-at)^3 + as.numeric(at <= 1)*(4-6*at^2+3*at^3))
}

# evaluate b spline

bspline <- function(x,cuts) {

  n = length(cuts)-1
  res = matrix(0,nrow=length(x),ncol=n+3)
  h = cuts[2]-cuts[1]
  if(h > 0) {
    for(k in 1:(n+3)) {
      res[,k] = phi.b( (x-cuts[1])/h - (k-2) )
    }
  } else {
    res[,1] = cuts[1]
  }

  return(res)
}

# code to test fitiv C code

# fit inclusive value regression for heterogeneous households

fitiv.hh.r <- function(param,data,info,index,nsave,nhhtodo,compare=FALSE,testiv=list(0)) {

  idlist = unique(data$PANID)

  prcols = which(substr(colnames(data),1,7)=="Product")

  nsize = info$nsize
  bigJ = info$bigJ
  nbrand = info$nbrand

  ncutiv = info$ncutiv

  nivdeprow = 2*ncutiv+5

  nrowiv = ncutiv+4

  compvec = vector("list",length=nhhtodo)

  totobs = nrow(data)

  obsoffset = 0

  dropmat = matrix(FALSE,nrow=nhhs,ncol=nsize)

  cnum = rep(0,nhhs)

  iv.all = NULL
  bigx.all = NULL

  for(hh in 1:nhhtodo) {

    cutmat = matrix(0,nrow=ncutiv+1,ncol=nsize)

    if(compare) {

      ccutmat = matrix(testiv[[2]][1:((ncutiv+1)*nsize)+((ncutiv+1)*nsize)*(index-1)+(ncutiv+1)*nsize*nsave*(hh-1)],nrow=ncutiv+1,ncol=nsize)

    }

    co = param[,hh]

    data1 = data[data$PANID==idlist[hh],]

    nobs = nrow(data1)

    iv = matrix(0,nrow=nobs,ncol=nsize)

    a = rep(0,nsize)
    b = rep(0,nsize)
    h = rep(0,nsize)

    dropsize = rep(FALSE,nsize)

    for(j in 1:nsize) {

      ubig = data1[,prcols]*co[nbrand+3] + t(matrix(co[1:nbrand],nrow=nbrand,ncol=nobs))
      ubig = exp(ubig)*( data1[,prcols] <= 999 & t(matrix(brsize[,2]==j,nrow=nbrand,ncol=nobs)) )

      iv[,j] = log(apply(ubig,1,sum))

      a[j] = min(iv[,j])
      b[j] = max(iv[,j])
      h[j] = (b[j]-a[j])/(ncutiv+1)
      cutmat[,j] = a[j] + (0:ncutiv)*h[j]

      if(a[j] == b[j] | length(table(iv[,j])) <= 5 ) {dropsize[j] = TRUE}

    }

    dropmat[hh,] = dropsize

    assign("cutmat",cutmat,envir=globalenv())

    # construct b spline basis

    bigx = NULL

    for(j in 1:nsize) {
      #if(dropsize[j]) {
      #  bigx = cbind(bigx,rep(cutmat[1,j],nrow(ubig)),matrix(0,nrow=nrow(ubig),ncol=ncutiv+3))
      #} else {
      bigx = cbind(bigx,bspline(iv[,j],a[j] + (0:(ncutiv+1))*h[j]))
      #}
    }

    assign("iv",iv,envir=globalenv())
    assign("bigx",bigx,envir=globalenv())

    iv.all = rbind(iv.all,iv)
    bigx.all = rbind(bigx.all,bigx)

    # do autocorrelated regression

    bigx.l = bigx[1:(nrow(bigx)-1),]

    ivcoef = matrix(0,nrow=nrowiv*nsize,ncol=nsize)
    residmat = matrix(0,nrow=nobs-1,ncol=nsize)
    civcoef = matrix(0,nrow=nrowiv*nsize,ncol=nsize)

    ivblocksize = nrowiv*nsize*nsize*nsave

    for(j in 1:nsize) {
      if(dropsize[j]) {
        ivcoef[nrowiv*(j-1)+1,j] = 1
      } else {
        if(info$splinetype == 3) {
          xcols = nrowiv*(j-1)+1:nrowiv
        } else if(info$splinetype == 2) {
          xcols = NULL
          for(k in 1:nsize) {
            if(k==j) {
              xcols = c(xcols,nrowiv*(k-1)+1:nrowiv)
            } else if(!dropsize[k]) {
              xcols = c(xcols,nrowiv*(k-1)+1:(nrowiv-1))
            }
          }

        } else {
          xcols = 1:ncol(bigx)
        }

        cnum[hh] = kappa(t(bigx.l[,xcols])%*%bigx.l[,xcols])

        res = lm(iv[2:nrow(bigx),j]~bigx.l[,xcols]-1)
        ivcoef[xcols,j] = res$coefficients

        residmat[,j] = res$residuals
        if(compare) {
          offset = nsize*nrowiv*(j-1)+nsize*nrowiv*nsize*(index-1)+ivblocksize*(hh-1)
          civcoef[,j] = testiv[[1]][offset+1:(nrowiv*nsize)]
        }
      }
    }
    #print(hh)
    #if(hh==166) {
    #  browser()
    #}
    rvar = var(residmat)
    if(sum(dropsize) > 0) {
      vmean = mean(diag(rvar[!dropsize,!dropsize]))
      for(i in 1:nsize) {
        if(dropsize[i]) {
          rvar[i,] = 0
          rvar[,i] = 0
          rvar[i,i] = vmean
        }
      }
    }
    #if(hh==298) {
    #  browser()
    #}
    ivvari = chol2inv(chol(rvar))

    if(compare) {
      civvari = matrix(testiv[[3]][1:(nsize*nsize)+(nsize*nsize)*(index-1)+nsave*nsize*nsize*(hh-1)],nrow=nsize,ncol=nsize)
    }

    detiv = sqrt(det(rvar))

    if(compare) {
      cdetiv = testiv[[4]][index+nsave*(hh-1)]

      civ = matrix(0,nrow=nobs,ncol=nsize)

      for(j in 1:nsize) {
        civ[,j] = testiv[[6]][1:nobs + obsoffset + totobs*(j-1)]
      }

      # extract the C versions for comparison

      compvec[[hh]] = vector("list",length=6)

      compvec[[hh]][[1]] = list(ivcoef,civcoef)

      compvec[[hh]][[2]] = list(cutmat,ccutmat)

      compvec[[hh]][[3]] = list(ivvari,civvari)

      compvec[[hh]][[4]] = list(detiv,cdetiv)

      compvec[[hh]][[5]] = list(iv,civ)

    } else {

      compvec[[hh]] = vector("list",length=6)

      compvec[[hh]][[1]] = ivcoef

      compvec[[hh]][[2]] = cutmat

      compvec[[hh]][[3]] = ivvari

      compvec[[hh]][[4]] = detiv

      compvec[[hh]][[5]] = iv

    }

    obsoffset = obsoffset + nrow(data1)

  }

  assign("iv.all",iv.all,envir=globalenv())
  assign("bigx.all",bigx.all,envir=globalenv())

  return(list(compvec,dropmat,cnum))

}

# version assuming all households are homogeneous

fitiv.r <- function(param,data,info,index,nsave,nhhtodo,testiv,ivbig) {

  idlist = unique(data$PANID)

  prcols = which(substr(colnames(data),1,7)=="Product")

  nsize = info$nsize
  bigJ = info$bigJ
  nbrand = info$nbrand

  ncutiv = info$ncutiv

  nrowiv = ncutiv+4

  compvec = vector("list",length=nhhtodo)

  totobs = nrow(data)

  obsoffset = 0

  # compute inclusive value

  cutmat = matrix(0,nrow=ncutiv+1,ncol=nsize)

  co = param[,1]

  iv = matrix(0,nrow=totobs,ncol=nsize)

  a = rep(0,nsize)
  b = rep(0,nsize)
  h = rep(0,nsize)

  for(j in 1:nsize) {

    ubig = data[,prcols]*co[nbrand+3] + t(matrix(co[1:nbrand],nrow=nbrand,ncol=totobs))
    ubig = exp(ubig)*( data[,prcols] <= 999 & t(matrix(info$brsize[,2]==j,nrow=nbrand,ncol=totobs)) )

    iv[,j] = log(apply(ubig,1,sum))

    a[j] = min(iv[,j])
    b[j] = max(iv[,j])
    h[j] = (b[j]-a[j])/(ncutiv+1)
    cutmat[,j] = a[j] + (0:ncutiv)*h[j]

  }

  assign("cutmat",cutmat,envir=globalenv())

  # construct b spline basis

  bigx = NULL

  for(j in 1:nsize) {
    bigx = cbind(bigx,bspline(iv[,j],a[j] + (0:(ncutiv+1))*h[j]))
  }

  assign("iv",iv,envir=globalenv())
  assign("bigx",bigx,envir=globalenv())

  # do autocorrelated regression

  lpanid = c(0,data$PANID[1:(nrow(data)-1)])

  first = data$PANID != lpanid

  nregobs = sum(!first)

  bigx.l = bigx[c(!first[2:length(first)],FALSE),]

  ivcoef = matrix(0,nrow=nrowiv*nsize,ncol=nsize)
  residmat = matrix(0,nrow=nregobs,ncol=nsize)
  civcoef = matrix(0,nrow=nrowiv*nsize,ncol=nsize)

  ivblocksize = nrowiv*nsize*nsize*nsave

  ivres = vector("list",length=nsize)

  for(j in 1:nsize) {
    #browser()
    if(info$splinetype == 3) {
      res = lm(iv[!first,j]~bigx.l[,nrowiv*(j-1)+1:nrowiv]-1)
      ivcoef[nrowiv*(j-1)+1:nrowiv,j] = res$coefficients
    } else if(info$splinetype == 2) {
      xcols = NULL
      for(k in 1:nsize) {
        if(k==j) {
          xcols = c(xcols,nrowiv*(k-1)+1:nrowiv)
        } else {
          xcols = c(xcols,nrowiv*(k-1)+1:(nrowiv-1))
        }
      }

      res = lm(iv[!first,j]~bigx.l[,xcols]-1)
      ivcoef[xcols,j] = res$coefficients
    } else {
      res = lm(iv[!first,j]~bigx.l-1)
      ivcoef[,j] = res$coefficients
    }
    residmat[,j] = res$residuals
    ivres[[j]] = res

  }

  rvar = var(residmat)
  ivvari = chol2inv(chol(rvar))

  detiv = sqrt(det(rvar))

  # extract the C versions for comparison

  civ = ivbig

  for(hh in 1:nhhtodo) {

    ccutmat = matrix(testiv[[2]][1:((ncutiv+1)*nsize)+((ncutiv+1)*nsize)*(index-1)+(ncutiv+1)*nsize*nsave*(hh-1)],nrow=ncutiv+1,ncol=nsize)

    for(j in 1:nsize) {
      offset = nsize*nrowiv*(j-1)+nsize*nrowiv*nsize*(index-1)+ivblocksize*(hh-1)
      civcoef[,j] = testiv[[1]][offset+1:(nrowiv*nsize)]
    }

    civvari = matrix(testiv[[3]][1:(nsize*nsize)+(nsize*nsize)*(index-1)+nsave*nsize*nsize*(hh-1)],nrow=nsize,ncol=nsize)

    cdetiv = testiv[[4]][index+nsave*(hh-1)]

    #civ = matrix(0,nrow=totobs,ncol=nsize*bigJ)

    #for(j in 1:nsize) {
    #  civ[,j] = testiv[[6]][1:totobs + obsoffset + totobs*(j-1)]
    #}

    compvec[[hh]] = vector("list",length=6)

    compvec[[hh]][[1]] = list(ivcoef,civcoef)

    compvec[[hh]][[2]] = list(cutmat,ccutmat)

    compvec[[hh]][[3]] = list(ivvari,civvari)

    compvec[[hh]][[4]] = list(detiv,cdetiv)

    #compvec[[hh]][[5]] = list(ivdepcoef,civdepcoef)

  }

  return(list(compvec,list(iv,civ),ivres,bigx.l,first))

}
