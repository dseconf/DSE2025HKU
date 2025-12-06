#define ARMA_DONT_USE_CXX11
#ifdef defined __unix__
  #include <sys/resource.h>
#endif
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#include <cmath>
#include <Rcpp/Benchmark/Timer.h>
#include <iostream>
#include <fstream>
namespace RcPar = RcppParallel;


int debugprint=0;

double pi = 3.14159265359;

struct llinputs {
  int nco;
  int nobs;
  int nbrand;
  int nsize;
  int nhhs;
  int sptype;
  int ncutiv;
  int inttype;
  int nsave;
  int ntilde;
  int nrgrid;
  int necoef;
  int capj;
  int nsim;
  int ninitt;
  int retinv;
  int myopic;
  int cmodel;
  int debug;
  int usevf;
  int idrawtype;
  int first20;
  int hinvbound;
  int genrnd;
  int nd;
  int nq;
  int ncut;
  int pflag;
  int hhprint;
  int repprint;
  int ncpgrid;
  int rcpgrid;
  int ncrate;
  double dg;
  double maxinv;
  int hmodel;
  int usecdraws;
  int initvf;
  int invmodel;
  int nb;
  int nbsmall;
  int maxbottles;
  int lomega;
  int nobsinit;
  int lencdraws;
  int lenicdraws;
  int cfrepkeep;
};

// simple max function

double dmax(double a,double b)
{
	if(a > b) {
		return(a);
	} else {
		return(b);
	}
}

double dmin(double a,double b)
{
	if(a < b) {
		return(a);
	} else {
		return(b);
	}
}

double dsign(double x) {
  return (x > 0) - (x < 0);
}

// Index checking function

void indcheck(int ind, int lb, int ub, const char * msg,...) {

  if(ind < lb || ind > ub-1) {
    va_list args;
    va_start(args,msg);

    char buf[100];
    vsnprintf(buf,100,msg,args);
    Rcpp::Rcout << buf << ". ";
    Rcpp::Rcout << "Index value: " << ind << "; Bounds: [" << lb << ", " << ub << "] .";
    throw std::range_error("Index out of bounds.");

    //error("%s. Index out of bounds.  Index value %d: Bounds: [%d,%d] \n",buf,ind,lb,ub);

    va_end(args);
  }

}

// compute min or max of vector

void maxvec(const double * x, int nx, double * mx, int * indm, int getmin) {

  *mx = x[0];
  *indm = 0;

  if(getmin) {
    for(int i=1;i<nx;i++) {
      if(x[i] < *mx) {
        *mx = x[i];
        *indm = i;
      }
    }
  } else {
    for(int i=1;i<nx;i++) {
      if(x[i] > *mx) {
        *mx = x[i];
        *indm = i;
      }
    }
  }

}

// translate inventory to inventory index

int itoiind(const double invprime, const double * istates, const int nistates, const int first20) {

  int iind = 0;
  double dif = 0;
  if(invprime <= istates[0]) {
    iind = 1;
  } else if (invprime >= istates[nistates-1]) {
    iind = nistates;
  } else if(invprime > istates[first20-1]) {
    dif = istates[first20] - istates[first20-1];
    iind = first20 + (int)((invprime-istates[first20-1])/dif);
  } else {
    dif = istates[1] - istates[0];
    iind = 1 + (int)((invprime-istates[0])/dif);
  }
  if(iind > nistates) {
    Rcpp::Rcout << "bug in inventory indexing:" << std::endl;
    Rcpp::Rcout << "invprime: " << invprime << std::endl;
    Rcpp::Rcout << "dif: " << dif << std::endl;
    Rcpp::Rcout << "first20: " << first20 << std::endl;
    Rcpp::Rcout << "istates[first20]" << istates[first20] << "; istates[first20-1]: " << istates[first20-1] << std::endl;
    throw std::range_error("vf interpolation problem");
  }
  return(iind);

}

// b spline phi function

double phifn(double t) {
  double at = fabs(t);

  if(at > 1 && at <= 2) {
    return( (2.0-at)*(2.0-at)*(2.0-at) );
  } else if (at <= 1) {
    return( 4.0-6.0*at*at+3.0*at*at*at );
  } else {
    return(0.0);
  }

}


// evaluate b spline

double bspline(double x, const double * spcoef, const double * cuts, int lcuts) {

  int n = lcuts-1;
  double res = 0.0;
  double h = cuts[1]-cuts[0];
  if(h==0) {
    res = spcoef[0]*x;
  } else {
    for(int k=0;k<n+3;k++) {
      //Rprintf("spcoef[%d]: %f; phifn: %f\n",k,spcoef[k],phifn( (x-cuts[0])/h - ((double)(k+1)-2.0) ));
      res += spcoef[k]*phifn( (x-cuts[0])/h - ((double)(k+1)-2.0) );
    }
  }

  return(res);
}

// compute b spline basis points

void bsplinei(double x, double * cuts, int lcuts, double * y) {

  int n = lcuts-1;
  double h = cuts[1]-cuts[0];
  if(h==0) {
    y[0] = cuts[0];
    for(int k=1;k<n+3;k++) {
      y[k] = 0;
    }
  } else {
    for(int k=0;k<n+3;k++) {
      y[k] = phifn( (x-cuts[0])/h - ((double)(k+1)-2.0) );
    }
  }

}

// ols wrapper for arma

void ols(const arma::mat &X, const arma::colvec &y, arma::colvec &coef, arma::colvec &resid, 
	 double * cnum, const int debugprint) {

  // this might be a little faster
  //coef = inv_sympd(X.t() * X) * X.t() * y;

  if(debugprint) {
    arma::mat xx = X.t()* X;
    *cnum = arma::cond(xx);
    //Rcpp::Rcout << "x'x: " << std::endl;
    //xx.print();
  }

  coef = arma::solve(X, y);      // fit model y ~ X
  resid = y - X*coef;            // residuals

}

// covariance matrix calculation

void covc(const arma::mat &x, const int n1, const int n2, arma::mat &cov) {

  double xmean[n2];

  for(int j=0;j<n2;j++) {
    xmean[j] = 0;
    for(int i=0;i<n1;i++) {
      xmean[j] += x(i,j);
    }
    xmean[j] /= (double)n1;
  }

  for(int i=0;i<n2;i++) {
    for(int j=i;j<n2;j++) {
      cov(i,j) = 0;
      for(int k=0;k<n1;k++) {
        cov(i,j) += (x(k,i)-xmean[i])*(x(k,j)-xmean[j]);
      }
      cov(i,j) /= (double)n1 - 1.0;
      cov(j,i) = cov(i,j);
    }
  }

}


// quantile function

double quantile(const Rcpp::NumericVector x, const double q) {

  Rcpp::NumericVector y = clone(x);
  std::sort(y.begin(),y.end());

  return y[x.size()*(q - 0.000000001)];

}

void tform(const Rcpp::NumericMatrix &par, const Rcpp::NumericVector &xfull,
           const Rcpp::NumericVector &tf, const Rcpp::IntegerVector &fixed,
           const Rcpp::IntegerVector &paramequal,
           const Rcpp::NumericMatrix &paramstart,
           const Rcpp::NumericVector &lbounds, const Rcpp::NumericVector &ubounds,
           const Rcpp::NumericVector &crate, const int useparamstart, const int nbrand,
           const int npsmall, const int npbig, const int nhhs, const int crfix, const int datacrate,
	   const int crhhlb, const int sizeshifter, const int nsize, const int * brsize,
	   const int * sizebrand, Rcpp::NumericMatrix &parbig) {

  for(int k=0;k<nhhs;k++) {
    int j=-1;
    for(int i=0;i<npbig;i++) {

      double lb=lbounds[i];
      if(crhhlb && (i == nbrand || i == nbrand+1)) {
	lb = crate[k];
      }
      if(!fixed[i]) {
        j++;
        parbig[i+npbig*k] = par[j+npsmall*k];
        if(abs(tf[i]) == 1) {
          parbig[i+npbig*k] = dsign(tf[i])*exp(parbig[i+npbig*k]);
        } else if (abs(tf[i]) == 2) {
          parbig[i+npbig*k] = lb+(ubounds[i]-lb)*exp(parbig[i+npbig*k])/
            (1+exp(parbig[i+npbig*k]));
        }
      } else {
        if(paramequal[i] > 0) {
          if(paramequal[i]-1 >= i)
            throw std::runtime_error("Error (C++) paramequal[i] must be less than i");
          parbig[i+npbig*k] = parbig[paramequal[i]-1+npbig*k];
        } else {
          if(useparamstart) {
            parbig[i+npbig*k] = paramstart[i+npbig*k];
          } else {
            parbig[i+npbig*k] = xfull[i];
          }
        }
      }
    }
    if(crfix) {
      parbig[nbrand+1+npbig*k] = parbig[nbrand+npbig*k];
    }
    if(datacrate) {
      parbig[nbrand+1+npbig*k] = crate[k];
      parbig[nbrand+npbig*k] = crate[k];
    }
    if(sizeshifter) {
      for(int kk=0;kk<nbrand;kk++) {
	if(brsize[kk+nbrand] > 1 && sizebrand[kk]) {
	  parbig[kk+npbig*k] += parbig[npbig-1-nsize+brsize[kk+nbrand]+npbig*k];
	}
      }
    }
  }

}


// log likelihood of brand choices

void brchoicellhh(const Rcpp::NumericMatrix &co, const Rcpp::NumericVector &tunits,
                  const Rcpp::NumericVector &panid, const Rcpp::NumericMatrix &pricemat,
                  const Rcpp::IntegerVector &brindex, const Rcpp::IntegerVector &brsize, const int nobs,
                  const int nbrand, const int pflag, const int retindll,
                  const int nhhs, const int varflag, const int nco, Rcpp::NumericVector &rres) {


  double ut[nbrand];

  double ll = 0;

  int drop[nbrand];

  int hh = -1;

  for(int n=0;n<nobs;n++) {
    if(n==0 || panid[n] != panid[n-1]) {
      hh++;
      ll=0;
      if(!varflag)
        rres[hh]=0;
    }
    if(tunits[n] > 0) {
      double maxinv = -1000000;
      int first = 1;
      int printobs = n < -1;
      for(int j=0;j<nbrand;j++) {
        double price = pricemat[n+nobs*j];
        if(hh < 0)
          Rcpp::Rcout << " price [" << j << "]: " << price;
        if(price >= 999 || brsize[((int)brindex[n])-1+nbrand] != brsize[j+nbrand]) {
          drop[j] = 1;
        } else {
          //ut[j] = co[j+nco*hh] + co[nbrand+2+nco*hh]*price*tunits[n];
          // needs to be this with the new model
          ut[j] = co[j+nco*hh] + co[nbrand+2+nco*hh]*price;
          drop[j] = 0;
          if(first) {
            first = 0;
            maxinv = ut[j];
          } else {
            maxinv = std::max(maxinv,ut[j]);
          }
          if(printobs) {
            Rcpp::Rcout << "obs " << n << ", ut[" << j << "]: " << ut[j] << std::endl;
          }
        }
      }// my "clever" idea here doesn't work properly if some parameters get huge
      double s = 0;
      for(int j=0;j<nbrand;j++) {
        if(!drop[j]) {
          //ut[j] -= maxinv;
          s += exp(ut[j]);
        }
      }

      //ll += ut[((int)brindex[n])-1];

      //ll -= log(s);
      ll += log(exp(ut[((int)brindex[n])-1])/s);
      if(varflag) {
        rres[n] = ut[((int)brindex[n])-1]-log(s);
      } else {
        rres[hh] = ll;
      }

      if(hh < 0) {
        Rcpp::Rcout << std::endl << "br ll info, " << hh << ", " << n << ", " << s << ", " << ll  << std::endl;
        Rcpp::Rcout << "drop";
        for(int j=0;j<nbrand;j++)
          Rcpp::Rcout << " " << drop[j];
        Rcpp::Rcout << std::endl;
      }

    }
  }

}

// fitiv for homogeneous model (not parallel, so hhinds should cover all households in sample)

void fitiv(const Rcpp::NumericMatrix &co,  const Rcpp::NumericVector &panid,
           const Rcpp::NumericMatrix &pricemat,
           const Rcpp::IntegerVector &hhinds, const Rcpp::IntegerVector &brsize,
           const Rcpp::IntegerVector &obshhinds, const int nco, const int nobs,
           const int nbrand, const int nsize, const int nhhs, const int sptype,
           const int ncutiv, Rcpp::NumericVector &riv,
           Rcpp::NumericVector &rcutmat, Rcpp::NumericVector &rivcoef,
           Rcpp::NumericVector &rivvari, Rcpp::NumericVector &rdetiv,
           const int debug = 0) {

  int hh = hhinds[0]-2;

  int obsend;
  if(hhinds[1]==nhhs) {
    obsend = nobs;
  } else {
    obsend = obshhinds[hhinds[1]];
  }

  // first, compute inclusive values

  int totobs = obsend-(obshhinds[hhinds[0]-1]-1);

  int liv = riv.size();

  double ut[nbrand];

  int drop[nbrand];

  for(int n=obshhinds[hhinds[0]-1]-1;n<obsend;n++) {
    if(n==obshhinds[hhinds[0]-1]-1 || panid[n] != panid[n-1]) {
      hh++;
    }

    for(int size=0;size<nsize;size++) {
      double umax = -10000000;
      for(int j=0;j<nbrand;j++) {
        double price = pricemat[n+nobs*j];
        if(price >= 999 || brsize[j+nbrand] != size+1) {
          drop[j] = 1;
        } else {
          ut[j] = co[j+nco*hh] + co[nbrand+2+nco*hh]*price;
          drop[j] = 0;
          if(ut[j] > umax)
            umax = ut[j];
        }
      }
      double s = 0;
      for(int j=0;j<nbrand;j++) {
        if(!drop[j]) {
          s += exp(ut[j]-umax);
        }
      }
      s=log(s)+umax;
      int indx = n-(obshhinds[hhinds[0]-1]-1)+totobs*size;
      indcheck(indx,0,liv,"iv size bug: n: %d",n);
      riv[indx] = s;
    }
  }

  int hhlen = hhinds[1]-(hhinds[0]-1);

  int nrowiv = 0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  //int ivblocksize = sptype == 2 && ncutiv > 0 ? nrowiv*nsize*nsize : nrowiv*nsize;

  int ivarriblocksize = nsize*nsize;

  int cutblocksize = (ncutiv+1)*nsize;

  int bobs = (obshhinds[hhinds[0]-1]-1);

  double tempb[nrowiv*nsize];

  //int begin1 = (int)begin;
  //int end1 = (int)end;

  int cutmsize = ncutiv == 0 ? 1 : nsize*(ncutiv+1);

  double cutpt[cutmsize];

  // spline cut points

  if(ncutiv > 0) {

    for(int k=0;k<nsize;k++) {
      double miniv, maxiv;
      int di;
      int indx = totobs*k;
      maxvec(riv.begin()+indx,totobs,&miniv,&di,1);
      maxvec(riv.begin()+indx,totobs,&maxiv,&di,0);
      double h = (maxiv-miniv)/((double)ncutiv+1.0);
      for(int hh=hhinds[0]-1;hh<hhinds[1];hh++) {
        for(int c=0;c<ncutiv+1;c++) {
          cutpt[c+(ncutiv+1)*k] = miniv + ((double)c)*h;
          indcheck(c+k*(ncutiv+1)+cutblocksize*(hh-(hhinds[0]-1)),0,hhlen*cutblocksize,"cutmat index bug.");
          rcutmat[c+k*(ncutiv+1)+cutblocksize*(hh-(hhinds[0]-1))] = cutpt[c+(ncutiv+1)*k];
        }
      }
    }
  }

  // construct matrices of x and y vars for regression

  int bigxcols = ncutiv == 0 ? nsize+1 : nrowiv*nsize;

  double bigx[bigxcols*totobs];
  //Rcpp::Rcout << "totobs-hhlen: " << totobs-hhlen << ", bigxcols: " << bigxcols << std::endl;
  arma::mat bigx1(totobs-hhlen,bigxcols);
  bigx1.zeros();

  arma::mat residmat(totobs-hhlen,nsize); residmat.zeros();

  arma::mat ymat(totobs-hhlen,nsize); ymat.zeros();

  int indx = 0;
  int indx1 = 0;

  for(int hh=hhinds[0]-1;hh<hhinds[1];hh++) {

    int obsstart = obshhinds[hh]-1;
    int obsend = (hh < nhhs-1) ? obshhinds[hh+1]-1 : nobs;

    int lobs = obsend-obsstart;

    // b spline basis

    for(int n=obsstart;n<obsend;n++) {
      if(ncutiv == 0) {
        bigx[indx+totobs*0] = 1.0;
        if(n < obsend-1)
          bigx1(indx1,0) = 1.0;
      }
      for(int k=0;k<nsize;k++) {
        indcheck(indx+totobs*k,0,liv,"iv bspline index bug.");
        if(n < obsend-1) {
          ymat(indx1,k) = riv[indx+1+totobs*k];
        }
        if(ncutiv == 0) {
          bigx[indx-obsstart+lobs*(1+k)] = riv[indx+totobs*k];
          if(n < obsend-1) {
            bigx1(indx1,1+k) = riv[indx+totobs*k];
          }
        } else {
          bsplinei(riv[indx+totobs*k], cutpt+(ncutiv+1)*k, ncutiv+2, tempb+k*nrowiv);
          for(int j=0;j<nrowiv;j++) {
            indcheck(indx+totobs*(j+k*nrowiv),0,bigxcols*totobs,"bigx bspline index bug.");
            bigx[indx+totobs*(j+k*nrowiv)] = tempb[j+k*nrowiv];
            if(n < obsend-1) {
              bigx1(indx1,j+k*nrowiv) = tempb[j+k*nrowiv];
              if(debug && n < 10) {
                Rcpp::Rcout << "in if (" << indx1 << ", " << j+k*nrowiv << "): " << bigx1(indx1,j+k*nrowiv) << std::endl;}
            }
          }
        }
      }
      indx++;
      if(n < obsend-1)
        indx1++;
    }

  }

  if(debug) {
    Rcpp::Rcout << "bigx els:" << std::endl;
    for(int i=0;i<10;i++) {
      Rcpp::Rcout << "[" << i+1 << "]:";
      for(int j=0;j<bigxcols;j++) {
        Rcpp::Rcout << " " << bigx[i+totobs*j];
      }
      Rcpp::Rcout << std::endl;
    }
    Rcpp::Rcout << "bigx1 els:" << std::endl;
    for(int i=0;i<10;i++) {
      Rcpp::Rcout << "[" << i+1 << "]:";
      for(int j=0;j<bigxcols;j++) {
        Rcpp::Rcout << " " << bigx1(i,j);
      }
      Rcpp::Rcout << std::endl;
    }
  }

  for(int k=0;k<nsize;k++) {
    //indcheck(obsstart+1-bobs+totobs*k,0,liv,"iv ols coeff index bug (lb).");
    //indcheck(obsstart+1-bobs+totobs*k+lobs-2,0,liv,"iv ols coeff index bug (ub).");
    arma::colvec riv1(totobs-hhlen);
    for(int j=0;j<totobs-hhlen;j++)
      riv1(j) = ymat(j,k);
    arma::colvec residout(totobs-hhlen); residout.zeros();
    //int ncoef = sptype == 3 ? nrowiv : bigxcols;
    int ncoef = bigxcols;
    if(sptype == 3) {
      ncoef = nrowiv;
    } else if(sptype == 2) {
      ncoef = nrowiv + (nsize-1)*(nrowiv-1);
    }
    arma::colvec coeftemp(ncoef); coeftemp.zeros();
    double cnum=0;
    if(sptype == 3) {
      ols(bigx1.cols(k*nrowiv,(k+1)*nrowiv-1), riv1, coeftemp, residout, &cnum, 0);
    } else {
      // including everything leads to collinearity problems
      arma::uvec xcols(ncoef);
      int cindx = 0;
      for(int kk=0;kk<nsize;kk++) {
	if(kk==k) {
	  for(int ii=0;ii<nrowiv;ii++) {
	    xcols(cindx) = ii+kk*nrowiv;
	    cindx++;
	  }
	} else {
	  for(int ii=0;ii<nrowiv-1;ii++) {
	    xcols(cindx) = ii+kk*nrowiv;
	    cindx++;
	  }
	}
      }
      ols(bigx1.cols(xcols), riv1, coeftemp, residout, &cnum, 0);
    }

    // might be slow?
    for(int j=0;j<totobs-hhlen;j++) {
      residmat(j,k) = residout(j);
    }

    if(ncutiv == 0) {

      for(int hh=hhinds[0]-1;hh<hhinds[1];hh++) {
        for(int j=0;j<bigxcols;j++) {
          indcheck(j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivblocksize,"iv coeff index bug (coeff assignment).");
          rivcoef[j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j);  //could use memcpy??
        }
      }

    } else {

      for(int hh=hhinds[0]-1;hh<hhinds[1];hh++) {
	int cindx = 0;
	for(int j=0;j<bigxcols;j++) {
	  indcheck(j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivblocksize,"iv coeff index bug (coeff assignment).");
	  if(sptype == 3) {
	    if(!(j >= k*nrowiv && j <= (k+1)*nrowiv-1)) {
	      rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = 0;  //could use memcpy??
	    } else {
	      rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j-k*nrowiv);  //could use memcpy??
	    }
	  } else if (sptype == 2) {

	    if( (j >= k*nrowiv && j <= (k+1)*nrowiv-1) | ( (j+1)%nrowiv != 0) ) {
	      rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(cindx);
	      cindx++;
	    } else {
	      rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = 0;
	    }
	    
	  } else {
	    rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j);  //could use memcpy??
	  }
	}
      }

    }

  }

  // compute inverse variance matrix and determinant for variance

  arma::mat covmat(nsize,nsize); covmat.zeros();
  arma::mat covi(nsize,nsize); covi.zeros();

  covc(residmat, totobs-hhlen, nsize, covmat);

  // check for zeros - occasionally can get zero on diagonal if no price variation

  for(int j=0;j<nsize;j++) {
    if(covmat(j,j) == 0)
      covmat(j,j) = 0.1;
  }
  //if(arma::rank( covmat ) < nsize) {
  //  covmat.print();
  //  throw std::range_error("covariance matrix in ols nonsingular");
  //}
  covi = arma::inv_sympd( covmat );

  double deti = sqrt(arma::det(covmat)); // double check we don't want inverse determinant here

  for(int hh=hhinds[0]-1;hh<hhinds[1];hh++) {

    indcheck(ivarriblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivarriblocksize,"ivvari index bug (lb).");
    indcheck(ivarriblocksize*(hh-(hhinds[0]-1))+nsize*nsize-1,0,hhlen*ivarriblocksize,"ivvari index bug (ub).");

    for(int j=0;j<nsize;j++) {
      for(int k=0;k<nsize;k++) {
        rivvari[ivarriblocksize*(hh-(hhinds[0]-1))+j+nsize*k] = covi(j,k);
      }
    }

    indcheck(hh-(hhinds[0]-1),0,hhlen,"detiv index bug.");

    rdetiv[hh-(hhinds[0]-1)] = deti;

  }

}



// parallel worker for fitivhh

struct fitivpar : public RcPar::Worker
{
  // source variables - Rcpp passes pointers so & isn't really needed
  RcPar::RVector<double> riv;  // this is an input but making it a const
  // screws up with arma
  const RcPar::RVector<int> obshhinds;
  const RcPar::RVector<int> hhinds;
  const int nobs;
  const int nbrand;
  const int nsize;
  const int nhhs;
  const int sptype;
  const int ncutiv;
  const int nrowiv;
  const int bobs;
  const int totobs;
  const int liv;
  const int ivblocksize;
  const int ivarriblocksize;
  const int cutblocksize;
  const int hhlen;
  const RcPar::RMatrix<int> dropmat;
  const RcPar::RVector<int> ivdrop;

  // destination vectors
  RcPar::RVector<double> rcutmat;
  RcPar::RVector<double> rivcoef;
  RcPar::RVector<double> rivvari;
  RcPar::RVector<double> rdetiv;
  RcPar::RVector<int> ivdrop2;

  // initialize with source and destination
  fitivpar(Rcpp::NumericVector riv, const Rcpp::IntegerVector obshhinds,
           const Rcpp::IntegerVector hhinds,
           const int nobs,
           const int nbrand, const int nsize, const int nhhs,
           const int sptype, const int ncutiv, const int nrowiv,
           const int bobs, const int totobs, const int liv,
           const int ivblocksize, const int ivarriblocksize,
           const int cutblocksize, const int hhlen,
	   const Rcpp::IntegerMatrix dropmat, const Rcpp::IntegerVector ivdrop,
           Rcpp::NumericVector rcutmat, Rcpp::NumericVector rivcoef,
           Rcpp::NumericVector rivvari, Rcpp::NumericVector rdetiv,
	   Rcpp::IntegerVector ivdrop2)
    : riv(riv), obshhinds(obshhinds), hhinds(hhinds),
    nobs(nobs), nbrand(nbrand),
    nsize(nsize), nhhs(nhhs), sptype(sptype), ncutiv(ncutiv), nrowiv(nrowiv),
    bobs(bobs), totobs(totobs), liv(liv),
    ivblocksize(ivblocksize), ivarriblocksize(ivarriblocksize),
    cutblocksize(cutblocksize), hhlen(hhlen), dropmat(dropmat), ivdrop(ivdrop),
    rcutmat(rcutmat), rivcoef(rivcoef), rivvari(rivvari),
      rdetiv(rdetiv), ivdrop2(ivdrop2)  {}

   // compute inclusive values for specified households
   void operator()(std::size_t begin, std::size_t end) {

     double tempb[nrowiv*nsize];

     int begin1 = (int)begin;
     int end1 = (int)end;

     int cutmsize = ncutiv == 0 ? 1 : nsize*(ncutiv+1);

     for(int hh=begin1;hh<end1;hh++) {

       int obsstart = obshhinds[hh]-1;
       int obsend = (hh < nhhs-1) ? obshhinds[hh+1]-1 : nobs;

       int redo = 1;
       int ivdrop1 = ivdrop[hh];
       ivdrop2[hh] = ivdrop[hh];
       
       int lobs = obsend-obsstart;

       double maxcnum = 0;
       
       arma::mat covmat(nsize,nsize); covmat.zeros();
       arma::mat covi(nsize,nsize); covi.zeros();

       while(redo) {

	 int bigxcols = (ncutiv == 0 || ivdrop1) ? nsize+1 : nrowiv*nsize;

	 double bigx[bigxcols*lobs];

	 double cutpt[cutmsize];

	 // compute break points

	 if(ncutiv > 0 && !ivdrop1) {

	   for(int k=0;k<nsize;k++) {
	     double miniv, maxiv;
	     int di;
	     int indx = obsstart-bobs+totobs*k;
	     indcheck(indx,0,liv,"iv size bug (maxvec lb): hh: %d, size: %d",hh,k);
	     indcheck(indx+lobs-1,0,liv,"iv size bug (maxvec ub): hh: %d, size: %d",hh,k);
	     maxvec(riv.begin()+indx,lobs,&miniv,&di,1);
	     maxvec(riv.begin()+indx,lobs,&maxiv,&di,0);
	     double h = (maxiv-miniv)/((double)ncutiv+1.0);
	     for(int c=0;c<ncutiv+1;c++) {
	       cutpt[c+(ncutiv+1)*k] = miniv + ((double)c)*h;
	       indcheck(c+k*(ncutiv+1)+cutblocksize*(hh-(hhinds[0]-1)),0,hhlen*cutblocksize,"cutmat index bug.");
	       rcutmat[c+k*(ncutiv+1)+cutblocksize*(hh-(hhinds[0]-1))] = cutpt[c+(ncutiv+1)*k];
	     }
	   }

	 }

	 // b spline basis

	 for(int n=obsstart;n<obsend;n++) {
	   if(ncutiv == 0 || ivdrop1)
	     bigx[n-obsstart+lobs*0] = 1.0;
	   for(int k=0;k<nsize;k++) {
	     indcheck(n-bobs+totobs*k,0,liv,"iv bspline index bug.");
	     if(ncutiv == 0 || ivdrop1) {
	       bigx[n-obsstart+lobs*(1+k)] = riv[n-bobs+totobs*k];
	     } else {
	       bsplinei(riv[n-bobs+totobs*k], cutpt+(ncutiv+1)*k, ncutiv+2, tempb+k*nrowiv);
	       for(int j=0;j<nrowiv;j++) {
		 indcheck(n-obsstart+lobs*(j+k*nrowiv),0,nrowiv*nsize*lobs,"bigx bspline index bug.");
		 bigx[n-obsstart+lobs*(j+k*nrowiv)] = tempb[j+k*nrowiv];
	       }
	     }
	   }
	 }

	 arma::mat bigx1(lobs-1,bigxcols);
	 bigx1.zeros();

	 for(int n=0;n<lobs-1;n++) {
	   for(int j=0;j<bigxcols;j++) {
	     bigx1(n,j) = bigx[n+lobs*j];
	   }
	 }

	 arma::mat residmat(lobs-1,nsize); residmat.zeros();
	 //Rcpp::Rcout << "aaa" << std::endl;
	 for(int k=0;k<nsize;k++) {
	   indcheck(obsstart+1-bobs+totobs*k,0,liv,"iv ols coeff index bug (lb).");
	   indcheck(obsstart+1-bobs+totobs*k+lobs-2,0,liv,"iv ols coeff index bug (ub).");
	   arma::colvec riv1(riv.begin()+obsstart+1-bobs+totobs*k,lobs-1,false);
	   arma::colvec residout(lobs-1); residout.zeros();
	   int ncoef = bigxcols;
	   if(sptype == 3) {
	     ncoef = nrowiv;
	   } else if(sptype == 2) {
	     int ndrop = 0;
	     for(int kk=0;kk<nsize;kk++) {
	       if(kk!=k) {
		 ndrop += dropmat(hh,kk);
	       }
	     }
	     if(ncutiv == 0 || ivdrop1) {
	       ncoef = nsize + 1 - ndrop;
	     } else {
	       ncoef = nrowiv + (nsize-1-ndrop)*(nrowiv-1);
	     }
	   }
	   arma::colvec coeftemp(ncoef); coeftemp.zeros();
	   if(dropmat(hh,k)) {
	     coeftemp(0) = 1.0;
	   } else {
	     double cnum = 0;
	     if(sptype == 3) {
	       ols(bigx1.cols(k*nrowiv,(k+1)*nrowiv-1), riv1, coeftemp, residout, &cnum, 0);
	     } else {
	       // including everything leads to collinearity problems
	       arma::uvec xcols(ncoef);
	       if(ncutiv == 0 || ivdrop1) {
		 xcols(0) = 0;
		 int cindx = 1;
		 for(int kk=0;kk<nsize;kk++) {
		   if(!dropmat(hh,kk)) {
		     xcols(cindx) = kk+1;
		     cindx++;
		   }
		 }
	       } else {
		 int cindx = 0;
		 for(int kk=0;kk<nsize;kk++) {
		   if(kk==k) {
		     for(int ii=0;ii<nrowiv;ii++) {
		       xcols(cindx) = ii+kk*nrowiv;
		       cindx++;
		     }
		   } else {
		     if(!dropmat(hh,kk)) {
		       for(int ii=0;ii<nrowiv-1;ii++) {
			 xcols(cindx) = ii+kk*nrowiv;
			 cindx++;
		       }
		     }
		   }
		 }
	       }
	       ols(bigx1.cols(xcols), riv1, coeftemp, residout, &cnum, 1);
	       if(hh==-1) {
		 Rcpp::Rcout << "hh: " << hh << "; size " << k << "; condition number: " << cnum << std::endl;
		 /*coeftemp.print();
		   bigx1.cols(xcols).print();
		   for(int n=obsstart;n<obsend;n++) {
		   Rcpp::Rcout << "iv[" << n << "]: ";
		   for(int k=0;k<nsize;k++) {
		   Rcpp::Rcout << riv[n-bobs+totobs*k];
		   }
		   Rcpp::Rcout << std::endl;
		   }
		   throw std::runtime_error("aaa");*/
	       }
	     }
	     if(cnum > maxcnum) {
	       maxcnum = cnum;
	     }
	     // might be slow?
	     for(int j=0;j<lobs-1;j++) {
	       residmat(j,k) = residout(j);
	     }
	   }

	   if(ncutiv == 0 || ivdrop1) {

	     //Rcpp::Rcout << "coefficients[" << k << "]: ";

	     for(int j=0;j<bigxcols;j++) {
	       indcheck(j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivblocksize,"iv coeff index bug (coeff assignment).");
	       //rivcoef[j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j);
	       int cindx = 0;
	       if(j==0) {
		 rivcoef[j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(cindx);  //could use memcpy??
	       } else if(!dropmat(hh,j-1)) {
		 cindx++;
		 rivcoef[j+bigxcols*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(cindx);
	       }
	       //Rcpp::Rcout << coeftemp(j) << ", ";
	     }
	     //Rcpp::Rcout << std::endl;

	   } else {
	     //Rcpp::Rcout << "ccc " << hh << std::endl;
	     int cindx = 0;
	     for(int j=0;j<bigxcols;j++) {
	       indcheck(j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivblocksize,"iv coeff index bug (coeff assignment).");
	       if(sptype == 3) {
		 if(!(j >= k*nrowiv && j <= (k+1)*nrowiv-1)) {
		   rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = 0;  //could use memcpy??
		 } else {
		   rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j-k*nrowiv);  //could use memcpy??
		 }
	       } else if (sptype == 2) {
		 int sindx = j/nrowiv+1;
		 //if(hh == 251) {
		 //Rcpp::Rcout << "j: " << j << ", sindx: " << sindx << ", cindx: " << cindx << ", ncoef: " << ncoef << std::endl;
		 //}
		 if( (j >= k*nrowiv && j <= (k+1)*nrowiv-1) | (( (j+1)%nrowiv != 0) & !dropmat(hh,sindx-1)) ) {
		   rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(cindx);
		   cindx++;
		 } else {
		   rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = 0;
		 }
	       } else {
		 rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))] = coeftemp(j);  //could use memcpy??
		 
	       }
	     }
	     //Rcpp::Rcout << "ddd" << std::endl;

	   }

	 }
	 //Rcpp::Rcout << "bbb" << std::endl;
	 // compute inverse variance matrix and determinant for variance


	 covc(residmat, lobs-1, nsize, covmat);

	 if(ivdrop1) {
	   redo = 0;
	 } else {
	   redo = maxcnum > 1e12;
	   ivdrop1 = redo;
	 }
	 
       }
       ivdrop2[hh] = ivdrop1;
       //for(int j=0;j<10;j++)
       //  Rcpp::Rcout << riv[j] << ", ";
       //Rcpp::Rcout << std::endl;
       

       //Rcpp::Rcout << residmat.row(1) << std::endl;

       //Rcpp::Rcout << covmat << std::endl;

       // check for zeros - occasionally can get zero on diagonal if no price variation

       double vmean = 0;
       int nkeep = 0;
       for(int j=0;j<nsize;j++) {
	 nkeep += !dropmat(hh,j);
	 if(!dropmat(hh,j)) {
	   vmean += covmat(j,j);
	 }
       }
       vmean /= (double)nkeep;

       for(int j=0;j<nsize;j++) {
	 if(dropmat(hh,j)) {
	   for(int k=0;k<nsize;k++) {
	     covmat(k,j) = 0;
	     covmat(j,k) = 0;
	   }
	   covmat(j,j) = vmean;
	 }
	 //if(covmat(j,j) == 0)
	 //  covmat(j,j) = 0.1;
       }
       /*if(arma::rank( covmat ) < nsize) {
	 
	 for(int j=0;j<bigxcols;j++) {
	   Rcpp::Rcout << "ivcoef(" << j << ", ): ";
	   for(int k=0;k<nsize;k++) {
	     Rcpp::Rcout << "  " << rivcoef[j+nsize*nrowiv*k+ivblocksize*(hh-(hhinds[0]-1))];
	   }
	   Rcpp::Rcout << std::endl;
	 }
         covmat.print();
         throw std::range_error("covariance matrix in ols singular");
	 }*/
       covi = arma::inv_sympd( covmat );

       double deti = sqrt(arma::det(covmat)); // double check we don't want inverse determinant here

       indcheck(ivarriblocksize*(hh-(hhinds[0]-1)),0,hhlen*ivarriblocksize,"ivvari index bug (lb).");
       indcheck(ivarriblocksize*(hh-(hhinds[0]-1))+nsize*nsize-1,0,hhlen*ivarriblocksize,"ivvari index bug (ub).");

       for(int j=0;j<nsize;j++) {
         for(int k=0;k<nsize;k++) {
           rivvari[ivarriblocksize*(hh-(hhinds[0]-1))+j+nsize*k] = covi(j,k);
         }
       }

       indcheck(hh-(hhinds[0]-1),0,hhlen,"detiv index bug.");

       rdetiv[hh-(hhinds[0]-1)] = deti;

     }

   }

};

// fit inclusive values

void fitivhh(const Rcpp::NumericMatrix &co,  const Rcpp::NumericVector &panid,
             const Rcpp::NumericMatrix &pricemat,
             const Rcpp::IntegerVector &hhinds, const Rcpp::IntegerVector &brsize,
             const Rcpp::IntegerVector &obshhinds, const int nco, const int nobs,
             const int nbrand, const int nsize, const int nhhs, const int sptype,
             const int ncutiv, const Rcpp::IntegerMatrix &dropmat, 
	     const Rcpp::IntegerVector &ivdrop, Rcpp::NumericVector &riv,
             Rcpp::NumericVector &rcutmat, Rcpp::NumericVector &rivcoef,
             Rcpp::NumericVector &rivvari, Rcpp::NumericVector &rdetiv,
	     Rcpp::IntegerVector &ivdrop2) {

  int hh = hhinds[0]-2;

  int obsend;
  if(hhinds[1]==nhhs) {
    obsend = nobs;
  } else {
    obsend = obshhinds[hhinds[1]];
  }

  // first, compute inclusive values

  //double sumutil;

  int totobs = obsend-(obshhinds[hhinds[0]-1]-1);

  int liv = riv.size();

  double ut[nbrand];

  int drop[nbrand];

  //std::ofstream outfile;

  //outfile.open("inclusivevalue.csv");

  for(int n=obshhinds[hhinds[0]-1]-1;n<obsend;n++) {
    if(n==obshhinds[hhinds[0]-1]-1 || panid[n] != panid[n-1]) {
      hh++;
    }

    //outfile << hh;

    for(int size=0;size<nsize;size++) {
      double umax = -10000000;
      for(int j=0;j<nbrand;j++) {
        double price = pricemat[n+nobs*j];
        //if(n==0) {
        //  Rcpp::Rcout << "price[" << n << ", " << j << "]: " << pricemat[n+nobs*j] << std::endl;
        //  Rcpp::Rcout << "brsize[" << n << ", " << j << "]: " << brsize[j+nbrand] << std::endl;
        //}
        if(price >= 999 || brsize[j+nbrand] != size+1) {
          drop[j] = 1;
        } else {
          ut[j] = co[j+nco*hh] + co[nbrand+2+nco*hh]*price;
          drop[j] = 0;
          if(ut[j] > umax)
            umax = ut[j];
          //if(n == 0)
          //  Rcpp::Rcout << "size " << size << " ut[" << j << "]: " << ut[j] << std::endl;
        }
      }
      double s = 0;
      for(int j=0;j<nbrand;j++) {
        if(!drop[j]) {
          s += exp(ut[j]-umax);
        }
      }
      s=log(s)+umax;
      int indx = n-(obshhinds[hhinds[0]-1]-1)+totobs*size;
      indcheck(indx,0,liv,"iv size bug: n: %d",n);
      riv[indx] = s;
      //outfile << ", " << riv[indx];
      //if(n==0)
      //  Rcpp::Rcout << "iv[" << size << "]: " << riv[indx] << std::endl;
    }
    //outfile << std::endl;
  }

  //outfile.close();

  int hhlen = hhinds[1]-(hhinds[0]-1);

  int nrowiv = 0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;

  int ivarriblocksize = nsize*nsize;

  int cutblocksize = (ncutiv+1)*nsize;

  int bobs = (obshhinds[hhinds[0]-1]-1);

  fitivpar ivpar(riv,obshhinds,hhinds,nobs,nbrand,nsize,nhhs,sptype,
                 ncutiv,nrowiv,bobs,totobs,liv,ivblocksize,
                 ivarriblocksize,cutblocksize,hhlen,dropmat,ivdrop,
                 rcutmat,rivcoef,rivvari,rdetiv,ivdrop2);

  parallelFor(hhinds[0]-1,hhinds[1],ivpar);

  //ivpar(hhinds[0]-1,hhinds[1]);

}

// compute euclidean distance

double edist(const double * x, const double * y, const int n) {

  double d = 0;

  for(int i=0;i<n;i++) {

    d += (x[i]-y[i])*(x[i]-y[i]);

  }

  return(sqrt(d));

}


// compute vf indexes where parameters are closest.

void getclosestind(const double * co, const double * cosave, const int * indexes, const int nco, const int nrep, const int nsave, const int rep, int * coorder) {

  double distvec[nrep];

  for(int i=0;i<nrep;i++) {
    distvec[i] = edist(co,cosave+(indexes[i]-1)*nco,nco);
  }

  if(nrep > 1 && rep > nrep) {
    //sort(distvec,nrep,coorder);
    //different sort, but keep the order stored in case of ties
    for(int i=0;i<nrep;i++) {
      int index = (rep-i-2)%nsave;
      indcheck(index,0,nsave,"coorder index error: %d, %d, %d",rep,index,i);
      if(i == 0) {
        coorder[0] = index;
      } else {
        coorder[i] = index;
        for(int j=0;j<i;j++) {
          indcheck(i-j-1,0,nsave,"coorder (i-j-1) index error: %d, %d",index,i);
          indcheck(i-j,0,nsave,"coorder (i-j) index error: %d, %d",index,i);
          if(distvec[index] < distvec[coorder[i-j-1]]) {
            int itemp = coorder[i-j-1];
            coorder[i-j-1] = index;
            coorder[i-j] = itemp;
          } else {
            break;
          }
        }
      }
    }

  } else {
    coorder[0] = 0;
  }
}

// multivariate normal kernel

double mvnormk(const double * d, const double * iv, const int nd) {

  double mm = 0, r=0;

  for(int i=0;i<nd;i++) {
    mm = 0;
    for(int j=0;j<nd;j++) {
      mm += d[j]*iv[j+nd*i];
    }
    r += mm*d[i];
  }

  return(r);

}

// compute kernel weights if we use standard IJC

void getkweights(const double * co, const double * cosave, const int * indexes,
                 const double * bwmatrix, const int nco, const int nrep, const int nsave, const int rep, int * coorder, double * kweights, int print) {

  double d[nco];

  if(nrep > 1 && rep > nrep) {
    //go back through indexes & get weights
    //I'm rescaling the weights by maxweight so that we
    // don't get errors with exp being zero
    double maxweight = -10000000;
    for(int i=0;i<nrep;i++) {
      int index = (rep-i-2)%nsave;
      indcheck(index,0,nsave,"coorder index error: %d, %d, %d",rep,index,i);
      coorder[i] = index;

      for(int j=0;j<nco;j++) {
        d[j] = co[j] - cosave[j+index*nco];  //this was indexes[i] but I think needs to be index??
      }

      kweights[i] = -0.5*mvnormk(d,bwmatrix,nco);
      maxweight = dmax(maxweight,kweights[i]);
      if(print) {
        for(int j=0;j<nco;j++) {
          Rcpp::Rcout << "index " << index << ", co[" << j << "]: " << co[j] << ", " << cosave[j+index*nco] << std::endl;
        }
        Rcpp::Rcout << "kernel kweights[" << i << "]: " << kweights[i] << std::endl;
      }

    }

    double ksum = 0.0;

    for(int i=0;i<nrep;i++) {
      kweights[i] = exp(kweights[i]-maxweight);
      ksum += kweights[i];
    }
    if(ksum == 0) {
      for(int i=0;i<nrep;i++) {
        kweights[i] = 0;
      }
    } else {
      for(int i=0;i<nrep;i++) {
        kweights[i] /= ksum;
      }
    }

  } else {
    coorder[0] = 0;
    kweights[0] = 1;
  }

}

// compute average predicted iv tomorrow given today's iv

void ivpred(const double * ivgrid, const double * ivcoef, const double * cutmat, int nd, int ncut, int sptype, int inttype, int capj, double * x) {

  if(inttype == 1) {
    for(int j=0;j<nd-1;j++) {
      if(ncut == 0) {
        x[j] = ivcoef[0+(3*ncut+3)*j] + ivcoef[1+(3*ncut+3)*j]*ivgrid[j];  //out of date
      } else {
        if(sptype == 1) {
          // need to update this
          int istart = 0;
          if(ivgrid[j] < ivcoef[0+(3*ncut+3)*j]) {
            istart = ncut;
          } else if (ivgrid[j] >= ivcoef[ncut-1+(3*ncut+3)*j]) {
            istart = 3*ncut;
          } else {
            int k = 0;
            while(ivgrid[j] >= ivcoef[k+(3*ncut+3)*j])
              k++;
            istart = ncut+2*k;
          }
          x[j] = ivcoef[istart+0+(3*ncut+3)*j] + ivcoef[istart+1+(3*ncut+3)*j]*ivgrid[j];
        } else {
          // this is up to date
          x[j] = 0;
          for(int k=0;k<nd-1;k++) {
            x[j] += bspline(ivgrid[k],ivcoef+(ncut+4)*k+(nd-1)*(ncut+4)*j,cutmat+k*(ncut+1),ncut+2);
          }
        }
      }
    }
  } else if (inttype == 2) {
    // monte carlo integration
    // here set qpts rows to number of sim draws and cols to num of dims
    for(int j=0;j<nd-1;j++) {
      if(ncut == 0) {
        x[j] = ivcoef[0+nd*j];
        for(int k=0;k<nd-1;k++)
          x[j] += ivcoef[1+k+nd*j]*ivgrid[k];
      } else {
        if(sptype == 1) {
          // need to update this
          int istart = 0;
          if(ivgrid[j] < ivcoef[0+(3*ncut+3)*j]) {
            istart = ncut;
          } else if (ivgrid[j] >= ivcoef[ncut-1+(3*ncut+3)*j]) {
            istart = 3*ncut;
          } else {
            int k = 0;
            while(ivgrid[j] >= ivcoef[k+(3*ncut+3)*j])
              k++;
            istart = ncut+2*k;
          }
          x[j] = ivcoef[istart+0+(3*ncut+3)*j] + ivcoef[istart+1+(3*ncut+3)*j]*ivgrid[j];
        } else {
          // this is up to date
          x[j] = 0;
          for(int k=0;k<nd-1;k++) {
            x[j] += bspline(ivgrid[k],ivcoef+(ncut+4)*k+(nd-1)*(ncut+4)*j,cutmat+k*(ncut+1),ncut+2);
          }
        }
      }
    }
  }

}



// compute importance sample function difference
// the detiv and dg should be square root of the determinant (they are
// defined this way in the fitiv procedure)

double impfn(double * ivpoint, const double * rgrid, int niv, const double * ivvari, const double * gvari, double detiv, double dg, const double * gmean) {

  double d1[niv], d2[niv];

  for(int i=0;i<niv;i++) {
    d1[i] = ivpoint[i] - rgrid[i];
    d2[i] = rgrid[i] - gmean[i];
  }

  return(exp(-0.5*(mvnormk(d1,ivvari,niv)-mvnormk(d2,gvari,niv)))*dg/detiv);

}

// importance sampler for gridmethod == 2

double impfn2(double * ivpoint, const double * rgrid, int niv, const double * ivvari, double detiv, double g) {

  double d1[niv], d2[niv];

  double denom = 1.0;

  for(int i=0;i<niv;i++) {
    d1[i] = ivpoint[i] - rgrid[i];
    denom *= 2.0*pi;
  }

  return( (exp(-0.5*mvnormk(d1,ivvari,niv))/(sqrt(denom)*detiv))/g );

}


// value function approximation for interpolation method 1 (Nearest neighbor)


// compute average value function at particular state point (I think ivpoint should be a fitted value?? faster to compute outside this loop)

double vfavgijchh(double * ivpoint, int iind, const double * vf, const double * rgrid, const int * indexes, int * coorder,
                  int ntilde, int nsave, int nrgrid, const double * ivvari, const double * gvari, const double detiv, const double dg, const double * gmean,
                  int niv, const double * igrid, int nistates, int first20, int b, const int nb, int print) {

  if(ntilde == 0) {
    return(0);
  } else {
    double avgv = 0.0;
    double wsum = 0.0;
    double w = 0.0;

    for(int i=0;i<ntilde;i++) {
      for(int j=0;j<nrgrid;j++) {
        w = impfn(ivpoint,rgrid+niv*(j+nrgrid*(indexes[coorder[i]]-1)),niv,ivvari,gvari,detiv,dg,gmean);

        indcheck(indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b)),0,nb*nsave*nrgrid*nistates,"vfavgichh index bug; coorder[i]: %d; indexes[i]: %d; j: %d; b: %d",
                 coorder[i],indexes[coorder[i]],j,b);
        avgv += vf[indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b))]*w;
        wsum += w;
      }
    }

    if(wsum == 0) {
      return(0);
    } else {
      return(avgv/wsum);
    }
  }

}


// version for gridmethod = 2.  Here, the weights are computed outside of the
// function. Note that in this version, we can update vf.
// I've also changed the order in which things are done, so that we can easily
// update the VF - the method should be equivalent to the old one though.
// note that if iter > 0, then we pass in the updated vf rather than the big
// vector of saved value functions.

double vfavgijchh2(double * ivtprobs, int iind, const double * vf, const double * rgrid, const int * indexes, int * coorder, double * kweights,
                  int ntilde, int nsave, int nrgrid, const double * ivvari, const double * gvari, const double detiv, const double dg, const double * gmean,
                   int niv, const double * igrid, int nistates, int first20, int b, const int nb, int print, int iter, int vfinterpmethod) {

  //if(print)
  //  Rcpp::Rcout << "aaa" << ntilde << std::endl;
  if(ntilde == 0) {
    return(0);
  } else {
    double avgv = 0.0;
    double wsum = 0.0;
    double w = 0.0;

    if(print)
      Rcpp::Rcout << "in vfavgijchh2, iind: " << iind << ", b: " << b << ", nsave: "<<
	nsave << ", nrgrid: " << nrgrid << ", nistates: " << nistates << std::endl;
    if(iter == 0) {

      for(int j=0;j<nrgrid;j++) {
        w = ivtprobs[j];
        double avgv1 = 0.0;
        //double w1 = 0;
        for(int i=0;i<ntilde;i++) {
	  indcheck(coorder[i],0,nsave,"vfavg coorder bug");
	  indcheck(indexes[coorder[i]]-1,0,nsave,"vfavg indexes bug");
          indcheck(indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b)),0,nb*nsave*nrgrid*nistates,"vfavgichh2 index bug; coorder[i]: %d; indexes[i]: %d; iind %d; nistates %d; j %d; b %d",
                   coorder[i],indexes[coorder[i]],iind,nistates,j,b);
          double ww = vfinterpmethod == 1 ? 1.0 : kweights[i];
          if(ww < 0)
            throw std::runtime_error("Negative weights - indexing issue in vfavgijchh2");
          avgv1 += vf[indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b))]*ww;
          if(print)
            Rcpp::Rcout << "index: "<< indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b)) <<
	      ", draw: " << i << ", index draw: " << indexes[coorder[i]] << ", inputs: vf: " << vf[indexes[coorder[i]]-1+nsave*(j+nrgrid*(iind+nistates*b))] << ", weight: " << ww << " ";
          //w1 += ww;
        }
        if(print)
          Rcpp::Rcout << std::endl;
        if(vfinterpmethod == 1) {
          avgv += avgv1*w/((double)ntilde);
        } else {
          avgv += avgv1*w;
        }
        //if(w1 > 0.0) {
        //  avgv += avgv1*w/w1;
        //}
        wsum += w;
        if(print) {
          Rcpp::Rcout << "grid point " << j << ", weight, vf:" << w << ", " << avgv1 << std::endl;
        }
      }

    } else {

      // for this to make any sense, we need to have a fixed grid
      // then ivtprobs is the same for all rows so we can just index by row 0

      for(int j=0;j<nrgrid;j++) {
        w = ivtprobs[j];
        avgv += vf[j+nrgrid*iind]*w;
        wsum += w;
      }

    }

    if(wsum == 0) {
      return(0);
    } else {
      return(avgv/wsum);
    }
  }

}

// storage cost function

double scost(const int * b, double inv, double packsize, double crate, const double * omega, int cmodel) {

  if(cmodel == 1) {
    double bnew = 0;
    bnew = inv/packsize;
    return( omega[0]*bnew + omega[1]*bnew*bnew );
  } else if (cmodel == 2) {
    int bcut = (int)omega[0];
    int bnew = 0;
    if(inv > 0) {
      bnew = (int)ceil(inv/packsize);
    }
    if(bnew > bcut) {
      return( omega[1]*((double)(bnew-bcut)) );
    } else {
      return(0.0);
    }
  } else if (cmodel == 3) {
    int bcut = (int)omega[0];
    int bnew = 0;
    if(inv > 0) {
      bnew = (int)ceil(inv/packsize);
    }
    if(bnew > bcut) {
      return( exp(omega[1]*log((double)(bnew-bcut))) );
    } else {
      return(0.0);
    }
  } else {
    double bnew = 0;
    if(inv > 0) {
      bnew = ceil(inv/packsize);
    } else {
      bnew = 0;
    }
    return( omega[0]*bnew + omega[1]*bnew*bnew );
  }

}

// storage cost function for multiple bottle size model (should we deprecate scost?)

double scost2(const double totalvol, const double * omega, int cmodel) {

  if(cmodel <= 1) {
    return( omega[0]*totalvol + omega[1]*totalvol*totalvol );
  } else if (cmodel == 2) {
    if(totalvol > omega[0]) {
      return( omega[1]*(totalvol-omega[0]) );
    } else {
      return(0.0);
    }
  } else if (cmodel == 3) {
    if(totalvol > omega[0]) {
      return( exp(omega[1]*log(totalvol-omega[0])) );
    } else {
      return(0.0);
    }
  } else {
    return( omega[0]*totalvol + omega[1]*totalvol*totalvol );
  }

}

// function to compute next period's inventory in the more complicated model
// note that here, inv is the amount left in the currently opened bottle (slot 1)

void nextinv(const int j, const int size, const int b, const double inv, const int * bstates, const int nb,
	     const int maxbottles, const double * vols, const double crate, const int * revarray,
	     const int nsize, const int* badmissable, int * nextbi, double * nexti, int * nextbflag,
	     const int printflag) {

  double binvprime = -crate;

  int currentb[maxbottles];
  int nextb[maxbottles];
  for(int i=0;i<maxbottles;i++) {
    currentb[i] = bstates[i+maxbottles*b];
  }
  
  if(b > 0) {
    binvprime = inv - crate;
  }

  int overfull = 0;
  
  if(b == 0 || (binvprime < 0 && currentb[1] == 0)) {
    // case 1: no inventory or you use up your last bottle
    if(j == 0) {
      for(int i=0;i<maxbottles;i++) {
	nextb[i] = 0;
      }
      *nexti = 0;
    } else {
      for(int i=0;i<maxbottles;i++) {
	if(i < j) {
	  nextb[i] = size+1;
	} else {
	  nextb[i] = 0;
	}
      }
      *nexti = vols[size]+binvprime;
    }
  } else {
    // case 2: you still have inventory after consumption
    if(binvprime > 0) {
      // case 2a: you don't use up your current bottle
      *nexti = binvprime;
      int zind = -1;
      for(int i=0;i<maxbottles;i++) {
	nextb[i] = currentb[i];
	if(currentb[i] == 0 && zind < 0) {
	  zind = i;
	}
      }
      if(zind >= 0 && j > 0) {
	overfull = zind+j-maxbottles;
	int slot = 0;
	for(int i=zind;i<maxbottles;i++) {
	  if(slot < j) {
	    nextb[i] = size+1;
	  } else {
	    nextb[i] = 0;
	  }
	  slot++;
	}
      }
    } else {
      // case 2b: you do use up your current bottle
      int zind = -1;
      for(int i=1;i<maxbottles;i++) {
	nextb[i-1] = currentb[i];
	if(currentb[i] == 0 && zind < 0) {
	  zind = i-1;
	}
      }
      if(zind == -1 && j > 0) {
	nextb[maxbottles-1] = size+1;  //this only happens if you were full
	overfull = j-1;
      } else if (zind >= 0) {
	int slot = 0;
	overfull = zind+j-maxbottles;
	for(int i=zind;i<maxbottles;i++) {
	  if(slot < j) {
	    nextb[i] = size+1;
	  } else {
	    nextb[i] = 0;
	  }
	  slot++;
	}
      } else {
	nextb[maxbottles-1] = 0;
      }
      *nexti = vols[nextb[0]-1]+binvprime;
    }
  }
  // map nextb into an index

  int indx = 0;
  int maxi=1;
  for(int i=maxbottles-1;i>=0;i--) {
    indx*=nsize+1;
    indx+=nextb[i];  //this will be 0,1,...,nsize, so this is right
    maxi*=nsize+1;
  }
  if(indx < 0 || indx > maxi-1 || printflag) {
    Rcpp::Rcout << "j: " << j << "; size: " << size << std::endl;
    Rcpp::Rcout << "inv: " << inv << "; crate: " << crate << std::endl;
    for(int i=0;i<maxbottles;i++) {
      Rcpp::Rcout << "currentb[" << i << "]: " << currentb[i] << "; ";
    }
    Rcpp::Rcout << std::endl;
    for(int i=0;i<maxbottles;i++) {
      Rcpp::Rcout << "nextb[" << i << "]: " << nextb[i] << "; ";
    }
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "overfull: " << overfull << std::endl;
    indcheck(indx,0,maxi,"bug in reverse b state lookup.");
  }
  // check if the next state is admissable
  *nextbflag = badmissable[revarray[indx]-1];
  if(!badmissable[revarray[indx]-1]) {
    // loop backwards until we find an admissable state
    int mbmax = maxbottles-2;
    for(int mbnew=mbmax;mbnew>=0;mbnew--) {
      indx=0;
      for(int i=mbnew;i>=0;i--) {
	indx*=nsize+1;
	indx+=nextb[i];  //this will be 0,1,...,nsize, so this is right
      }
      if(badmissable[revarray[indx]-1]) {
	break;
      }
    }
    // do we need to reset nexti??
  }
  *nextbflag = *nextbflag && overfull <= 0;
  if(indx < 0 || indx > maxi-1) {
    Rcpp::Rcout << "j: " << j << "; size: " << size << std::endl;
    Rcpp::Rcout << "inv: " << inv << "; crate: " << crate << std::endl;
    for(int i=0;i<maxbottles;i++) {
      Rcpp::Rcout << "currentb[" << i << "]: " << currentb[i] << "; ";
    }
    Rcpp::Rcout << std::endl;
    for(int i=0;i<maxbottles;i++) {
      Rcpp::Rcout << "nextb[" << i << "]: " << nextb[i] << "; ";
    }
    Rcpp::Rcout << std::endl;
    indcheck(indx,0,maxi,"bug in reverse b state lookup.");
  }
  *nextbi = revarray[indx]-1;
  // I was occasionally getting errors where i was really tiny at b=0
  if(*nextbi == 0) {
    *nexti = 0;
  }

}


// wrapper for testing out nextinv

// [[Rcpp::export]]
Rcpp::List nextinvwrap(const int j, const int size, const int b, const double inv,
		       Rcpp::IntegerMatrix bstates, const int nb,
		       const int maxbottles, Rcpp::NumericVector vols, const double crate,
		       Rcpp::IntegerVector revarray, const int nsize, Rcpp::IntegerVector badmissable,
		       const int printflag) {

  int nextbi;
  double nexti;
  int nextbflag;
  
  nextinv(j,size-1,b-1,inv,bstates.begin(), nb, maxbottles, vols.begin(), crate,
	  revarray.begin(), nsize, badmissable.begin(), &nextbi, &nexti, &nextbflag,
	  printflag);

  Rcpp::List ret;
  ret["nextb"] = nextbi+1;
  ret["nexti"] = nexti;
  ret["nextbflag"] = nextbflag;

  return(ret);
  
}

// inclusive value utility for not buying

double utiliv(int j, const int * b, const double ivgrid, double inv, double packsize, double crate, double gamma, double nu, double ccost, const double * omega, int cmodel,
              int printflag) {
  double invprime = inv + packsize*((double)j) - crate;
  double u = ivgrid + (j > 0 ? ccost : 0.0);

  //if(j>0)
  //  u += ivgrid[j-1] + ccost;
  if(invprime < 0)
    invprime = 0;
  if( inv + packsize*((double)j) < crate) {
    return(u+gamma+(-gamma-nu)*(crate - (inv + packsize*((double)j)) )/crate - scost(b,invprime,packsize,crate,omega,cmodel));
  } else {
    return(u+gamma-scost(b,invprime,packsize,crate,omega,cmodel));
  }
}

// inclusive value util with new state space

double utiliv2(const int j, const int size, const double ivgrid, const double crate,
	       const double gamma, const double nu, const double ccost, const double * omega,
               const int lomega, const int cmodel,
	       const double inv, const int b, const int nextb,  const int nextbflag, const int * bstates, const int nb,
	       const int maxbottles, const double * vols, const int printflag) {

  int currentb[maxbottles];
  int bprime[maxbottles];

  double currentvol = 0;
  double totalvol = 0;
  
  for(int i=0;i<maxbottles;i++) {
    currentb[i] = bstates[i+maxbottles*b];
    if(currentb[i] > 0) {
      if(i==0) {
	currentvol += inv;
      } else {
	currentvol += vols[currentb[i]-1];
      }
    }
    bprime[i] = bstates[i+maxbottles*nextb];
    if(bprime[i] > 0) {
      totalvol += vols[bprime[i]-1];
    }
  }
  double nextvol = currentvol + ((double)j)*vols[size];
  
  double u = ivgrid + (j > 0 ? ccost : 0.0);

  if(!nextbflag) {
    u += omega[lomega-1];
    //if(printflag) {
    //  Rcpp::Rcout << "omega[2]: " << omega[lomega-1] << std::endl;
    //}
  }

  // we can put in totalvol or nextvol
  // 
  if( nextvol < crate ) {
    return(u+gamma+(-gamma-nu)*(crate - nextvol )/crate - scost2(totalvol,omega,cmodel));
  } else {
    if(printflag) {
      Rcpp::Rcout << "size, j: " << size << ", " << j;
      Rcpp::Rcout << "; currentb[" << b << "]: ";
      for(int i=0;i<maxbottles;i++) {
	Rcpp::Rcout << " " << bstates[i+maxbottles*b];
      }
      Rcpp::Rcout << "; nextb[" << nextb << "]: ";
      for(int i=0;i<maxbottles;i++) {
	Rcpp::Rcout << " " << bstates[i+maxbottles*nextb];
      }
      Rcpp::Rcout << "; crate: " << crate << "; inv: " << inv << std::endl;
      Rcpp::Rcout << "u: " << u << ", scost(nextvol): " << scost2(nextvol,omega,cmodel) <<
	", scost(totalvol)" << scost2(totalvol,omega,cmodel) << std::endl;
    }
    return(u+gamma-scost2(totalvol,omega,cmodel));
  }
}


// quick vf linear interpolation at inventory point

double vfinterplinfast(double inv, const double * vf, int b, int bbig, const double * igrid, int nistates, int first20, 
                       const int * vffilled, const int invmodel, const int * bstates, const int nb,
		       const int maxbottles, const double * vols) {

  if(invmodel == 2) {
    if(b==0) {
      if(vffilled[0] < 0) {
	Rcpp::Rcout << "Error in state 0 vfinterp" << std::endl;
	throw std::range_error("vf interpolation problem");
      }
      return(vf[0]);
    } else {
      if(inv >= vols[bstates[0+maxbottles*bbig]-1]) {
	int vfindx = nistates*b+itoiind(inv,igrid,nistates,first20)-1;
	if(vffilled[vfindx] < 0) {
	  Rcpp::Rcout << "Error in state " << b << ", " << itoiind(inv,igrid,nistates,first20)-1 << " vfinterp" << std::endl;
	  throw std::range_error("vf interpolation problem");
	}
	return(vf[vfindx]);
      } else if (inv <= igrid[0]) {
	return(vf[nistates*b]);
      } else {
	/*double dif = 0;
	  int iind = 0;
	  if(inv > igrid[first20-1]) {
	  dif = igrid[first20] - igrid[first20-1];
	  iind = first20 + (int)((inv-igrid[first20-1])/dif);
	  } else {
	  dif = igrid[1] - igrid[0];
	  iind = 1 + (int)((inv-igrid[0])/dif);
	  }*/
	int iind = itoiind(inv,igrid,nistates,first20);
	if(vffilled[iind-1+nistates*b] < 0 || vffilled[iind+nistates*b] < 0) {
	  Rcpp::Rcout << "b: " << b << ", vol: " << vols[bstates[0+maxbottles*bbig]-1] << std::endl;
	  Rcpp::Rcout << "inv: " << inv << "; iind: " << iind << std::endl;
	  Rcpp::Rcout << "vffilled[" << iind-1+nistates*b << "]: " << vffilled[iind-1+nistates*b]
		      << "; vffilled[" << iind+nistates*b << "]: " << vffilled[iind+nistates*b]
		      << std::endl;
	  throw std::range_error("vf interpolation problem (invmodel2)");
	}
      
	double vfout = vf[iind-1+nistates*b]+(vf[iind+nistates*b]-vf[iind-1+nistates*b])*(igrid[iind] - inv)/(igrid[iind]-igrid[iind-1]);

	return(vfout);

      }

    }
  } else {
    if(inv >= igrid[nistates-1]) {
      return(vf[nistates*b+nistates-1]);
    } else if (inv <= igrid[0]) {
      return(vf[nistates*b]);
    } else {
      /*double dif = 0;
	int iind = 0;
	if(inv > igrid[first20-1]) {
	dif = igrid[first20] - igrid[first20-1];
	iind = first20 + (int)((inv-igrid[first20-1])/dif);
	} else {
	dif = igrid[1] - igrid[0];
	iind = 1 + (int)((inv-igrid[0])/dif);
	}*/
      int iind = itoiind(inv,igrid,nistates,first20);
      if(vffilled[iind-1+nistates*b] < 0 || vffilled[iind+nistates*b] < 0) {
	Rcpp::Rcout << "b: " << b << std::endl;
	Rcpp::Rcout << "inv: " << inv << "; iind: " << iind << std::endl;
	Rcpp::Rcout << "vffilled[" << iind-1+nistates*b << "]: " << vffilled[iind-1+nistates*b]
		    << "; vffilled[" << iind+nistates*b << "]: " << vffilled[iind+nistates*b]
		    << std::endl;
	throw std::range_error("vf interpolation problem");
      }
      
      double vfout = vf[iind-1+nistates*b]+(vf[iind+nistates*b]-vf[iind-1+nistates*b])*(igrid[iind] - inv)/(igrid[iind]-igrid[iind-1]);

      return(vfout);

    }
  }

}


// compute inclusive value at a given price vector for particular parameter
// this can be called when gridmethod == 2 and we are saving prices rather
// than inclusive values

void computeiv(const double * coef, const double * pricemat,
               const int * brsize, const int nsize, const int nbrand,
               double * riv) {

  double ut[nbrand];
  int drop[nbrand];

  for(int size=0;size<nsize;size++) {
    double umax = -10000000;
    for(int j=0;j<nbrand;j++) {
      double price = pricemat[j];
      if(price >= 999 || brsize[j+nbrand] != size+1) {
        drop[j] = 1;
      } else {
        ut[j] = coef[j] + coef[nbrand+2]*price;
        drop[j] = 0;
        if(ut[j] > umax)
          umax = ut[j];
      }
    }
    double s = 0;
    for(int j=0;j<nbrand;j++) {
      if(!drop[j]) {
        s += exp(ut[j]-umax);
      }
    }
    s=log(s)+umax;
    riv[size] = s;
  }

}


// parallel version of llijc
// note that if hmodel is 1, maxpind should be hmaxpind
// and pind should be hpriceindbig

struct llijcpar : public RcPar::Worker
{

  // input variables
  const RcPar::RMatrix<double> co;
  const RcPar::RVector<double> cosave;
  const RcPar::RVector<double> rgridsave;
  const RcPar::RVector<int> indexes;
  const RcPar::RVector<double> tunits;
  const RcPar::RVector<double> panid;
  const RcPar::RVector<int> brindex;
  const RcPar::RVector<int> brsize;
  const RcPar::RVector<int> obshhinds;
  const RcPar::RMatrix<double> iv;
  const RcPar::RVector<double> vfsave;
  const RcPar::RVector<double> cdraws;
  const RcPar::RVector<double> ic;
  const RcPar::RVector<double> iiv;
  const RcPar::RVector<double> ldraws;
  const RcPar::RVector<double> packsize;
  const RcPar::RVector<double> initx;
  const RcPar::RVector<int> ngrid;
  const RcPar::RVector<double> cutmat;
  const RcPar::RVector<double> ivcoef;
  const RcPar::RVector<double> ivvari;
  const RcPar::RVector<double> detiv;
  const RcPar::RVector<double> pgrid;
  const RcPar::RVector<double> qpts;
  const RcPar::RVector<double> gvari;
  const RcPar::RVector<double> gmean;
  const RcPar::RMatrix<double> rgrid;
  const RcPar::RVector<double> impdist;
  const RcPar::RMatrix<double> bwmatrix;
  const RcPar::RVector<int> pinds;
  const RcPar::RVector<int> bstates;
  const RcPar::RVector<int> revarray;
  const RcPar::RVector<int> badmissable;
  const RcPar::RVector<int> bindex;
  const RcPar::RVector<int> ivdrop;
  const int nobs;
  const int nbrand;
  const int nsize;
  const int inttype;
  const int nsave;
  const int ntilde;
  const int nrgrid;
  const int capj;
  const int nhhs;
  const int ninitt;
  const int ncutiv;
  const int retinv;
  const double maxinv;
  const int myopic;
  const int cmodel;
  const int debug;
  const int sptype;
  const int first20;
  const int hinvbound;
  const int nco;
  const int linitx;
  const int nd;
  const int ncut;
  const double dg;
  const int nistates;
  const double * istates;
  const int ivblocksize;
  const int cutblocksize;
  const int nrep;
  const int nreptilde;
  const double * initinv;
  const int * initb;
  const int nb;
  const int rep;
  const int gridmethod;
  const int vfinterpmethod;
  const int maxpind;
  const int hmodel;
  const int invmodel;
  const int maxbottles;
  const int lomega;

  // output variables
  RcPar::RVector<double> vfll;
  RcPar::RVector<double> utll;
  RcPar::RVector<double> rllrun;
  RcPar::RVector<double> risave1;
  RcPar::RVector<double> risave2;
  RcPar::RVector<double> llhh;

  // initialize with source and destination
  llijcpar(const Rcpp::NumericMatrix co, const Rcpp::NumericVector cosave,
           const Rcpp::NumericVector rgridsave, const Rcpp::IntegerVector indexes,
           const Rcpp::NumericVector tunits, const Rcpp::NumericVector panid,
           const Rcpp::IntegerVector brindex, const Rcpp::IntegerVector brsize,
           const Rcpp::IntegerVector obshhinds, const Rcpp::NumericMatrix iv,
           const Rcpp::NumericVector vfsave, const Rcpp::NumericVector cdraws,
           const Rcpp::NumericVector ic, const Rcpp::NumericVector iiv,
           const Rcpp::NumericVector ldraws, const Rcpp::NumericVector packsize,
           const Rcpp::NumericVector initx, const Rcpp::IntegerVector ngrid,
           const Rcpp::NumericVector cutmat, const Rcpp::NumericVector ivcoef,
           const Rcpp::NumericVector ivvari, const Rcpp::NumericVector detiv,
           const Rcpp::NumericVector pgrid, const Rcpp::NumericVector qpts,
           const Rcpp::NumericVector gvari, const Rcpp::NumericVector gmean,
           const Rcpp::NumericMatrix rgrid, const Rcpp::NumericVector impdist,
           const Rcpp::NumericMatrix bwmatrix, const Rcpp::IntegerVector pinds,
	   const Rcpp::IntegerVector bstates, const Rcpp::IntegerVector revarray,
	   const Rcpp::IntegerVector badmissable, const Rcpp::IntegerVector bindex,
	   const Rcpp::IntegerVector ivdrop,
           const int nobs, const int nbrand, const int nsize, const int inttype,
           const int nsave, const int ntilde, const int nrgrid, const int capj,
           const int nhhs, const int ninitt, const int ncutiv, const int retinv,
           const double maxinv, const int myopic, const int cmodel, const int debug,
           const int sptype, const int first20, const int hinvbound, const int nco,
           const int linitx, const int nd, const int ncut, const double dg,
           const int nistates, const double * istates, const int ivblocksize,
           const int cutblocksize, const int nrep, const int nreptilde,
           const double * initinv, const int * initb, const int nb, const int rep,
           const int gridmethod, const int vfinterpmethod, const int maxpind,
           const int hmodel, const int invmodel, const int maxbottles,
           const int lomega,
           Rcpp::NumericVector vfll, Rcpp::NumericVector utll,
           Rcpp::NumericVector rllrun, Rcpp::NumericVector risave1,
           Rcpp::NumericVector risave2, Rcpp::NumericVector llhh) :
  co(co), cosave(cosave), rgridsave(rgridsave), indexes(indexes), tunits(tunits),
    panid(panid), brindex(brindex), brsize(brsize), obshhinds(obshhinds),
    iv(iv), vfsave(vfsave), cdraws(cdraws), ic(ic), iiv(iiv), ldraws(ldraws),
    packsize(packsize), initx(initx), ngrid(ngrid), cutmat(cutmat),
    ivcoef(ivcoef), ivvari(ivvari), detiv(detiv), pgrid(pgrid), qpts(qpts),
    gvari(gvari), gmean(gmean), rgrid(rgrid), impdist(impdist),
    bwmatrix(bwmatrix), pinds(pinds), bstates(bstates), revarray(revarray),
    badmissable(badmissable), bindex(bindex), ivdrop(ivdrop),
    nobs(nobs), nbrand(nbrand), nsize(nsize),
    inttype(inttype), nsave(nsave), ntilde(ntilde), nrgrid(nrgrid), capj(capj),
    nhhs(nhhs), ninitt(ninitt), ncutiv(ncutiv), retinv(retinv), maxinv(maxinv),
    myopic(myopic), cmodel(cmodel), debug(debug), sptype(sptype), first20(first20),
    hinvbound(hinvbound), nco(nco), linitx(linitx), nd(nd), ncut(ncut), dg(dg),
    nistates(nistates), istates(istates), ivblocksize(ivblocksize),
    cutblocksize(cutblocksize), nrep(nrep), nreptilde(nreptilde),
    initinv(initinv), initb(initb), nb(nb), rep(rep),
    gridmethod(gridmethod), vfinterpmethod(vfinterpmethod),
    maxpind(maxpind), hmodel(hmodel), invmodel(invmodel), maxbottles(maxbottles),
    lomega(lomega),
    vfll(vfll), utll(utll),
    rllrun(rllrun), risave1(risave1), risave2(risave2), llhh(llhh) {}

  // compute log likelihood for specific households
  void operator()(std::size_t begin, std::size_t end) {

    //Rcpp::Rcout << "inside ll operator" << std::endl;

    int begin1 = (int)begin;
    int end1 = (int)end;

    double llrun = 0;
    //double ll = 0;

    //Rcpp::Rcout << "1+capj: " << 1+capj << ". nsize: " << nsize << ". nrgrid: " << nrgrid << std::endl;

    double ut[(1+capj)*nsize];
    memset(ut,0,(1+capj)*nsize*sizeof(double));
    int drop[(1+capj)*nsize];
    memset(drop,0,(1+capj)*nsize*sizeof(int));
    double vfhat[(1+capj)*nsize];
    memset(vfhat,0,(1+capj)*nsize*sizeof(double));

    int bchoice[(1+capj)*nsize];
    memset(bchoice,0,(1+capj)*nsize*sizeof(int));
    double ichoice[(1+capj)*nsize];
    memset(ichoice,0,(1+capj)*nsize*sizeof(double));

    double ivout[nsize];
    memset(ivout,0,nsize*sizeof(double));
    double ivgrid[(1+capj)*nsize];
    memset(ivgrid,0,(1+capj)*nsize*sizeof(double));
    double vftemp[nb*nistates*maxpind];
    memset(vftemp,0,nb*nistates*maxpind*sizeof(double));

    double rgridhh[(nd-1)*nrgrid];
    memset(rgridhh,0,(nd-1)*nrgrid*sizeof(double));
    double ivtprobs[nrgrid*maxpind];
    memset(ivtprobs,0,nrgrid*maxpind*sizeof(double));
    int ivtpfilled[maxpind];
    double invprime = 0;
    int coorder[nsave];
    memset(coorder,0,nsave*sizeof(int));

    double kweights[nsave];

    for(int i=0;i<nsave;i++)
      kweights[i] = -1;

    int lencdraws = cdraws.length();
    int leniv = nobs;

    int bprime = 0;
    int bsmall = 0;

    double clb = 0, cub = 0, gamma = 0, nu = 0, beta = 0, ccost = 0;
    double inv = 0;

    //int lomega = 2;  // update this later
    double omega[lomega];

    int vfind = 0;
    int b = 0;

    double invub = maxinv;
    //Rcpp::Rcout << "nb size: " << nb << ". nistates size: " << nistates << ". maxpind size: " << maxpind << std::endl;

    int vffilled[nb*nistates*maxpind];
    //Rcpp::Rcout << "aaa1" << std::endl;

    //for(int i=0;i<nistates;i++) {
    //  for(int j=0;j<nb;j++) {
    //    vffilled[i+nistates*j] = -1;
    //  }
    //}
    memset(vffilled,-1,nb*nistates*maxpind*sizeof(int));
    memset(ivtpfilled,-1,maxpind*sizeof(int));
    //Rcpp::Rcout << "bbb" << std::endl;

    for(int hh=begin1;hh<end1;hh++) {

      if(!hmodel) {
        memset(vffilled,-1,nb*nistates*maxpind*sizeof(int));
        memset(ivtpfilled,-1,maxpind*sizeof(int));
      }

      indcheck(hh,0,nhhs,"initinv. hh %d",hh);
      inv = initinv[hh];
      b = initb[hh];
      llrun = 0;
      llhh[hh] = 0;
      vfind = hh*nrgrid*nb*nistates*nsave;

      clb = co[nbrand+nco*hh];
      cub = co[nbrand+1+nco*hh];

      gamma = co[nbrand+3+nco*hh];
      nu = co[nbrand+4+nco*hh];
      beta = co[nbrand+5+nco*hh];
      ccost = co[nbrand+7+nco*hh];

      for(int j=0;j<lomega;j++) {
        omega[j] = co[nbrand+9+j+nco*hh];
      }

      invub = hinvbound ? omega[0] : maxinv;

      if(!myopic) {
        if(vfinterpmethod == 1) {
          getclosestind(co.begin() + nco*hh, cosave.begin() + nsave*nco*hh, indexes.begin(), nco, nrep, nsave, rep, coorder);
        } else {
          if(!(hmodel && hh > begin1)) {
            getkweights(co.begin() + nco*hh, cosave.begin() + nsave*nco*hh,
                        indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
                        rep, coorder, kweights, 0);
          }
        }
      }

      int obsend;
      if(hh == nhhs-1) {
        obsend = nobs;
      } else {
        obsend = obshhinds[hh+1]-1;
      }

      if(gridmethod == 2) {
        for(int p=0;p<nrgrid;p++) {
          indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);

          computeiv(co.begin()+nco*hh,rgrid.begin()+nbrand*p,brsize.begin(),
                    nsize,nbrand,rgridhh+(nd-1)*p);

        }
      }
      //Rcpp::Rcout << "rcnew: " << rcnew << std::endl;
      int nextbflag=0;

      for(int n=obshhinds[hh]-1;n<obsend;n++) {
        int units = (int)tunits[n];
        double rcnew;
        indcheck(n,0,lencdraws,"cdraws. n: %d",n);
        rcnew = cdraws[n];
        if(rcnew < 0.5) {
          rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
        } else {
          rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
        }
        //Rcpp::Rcout << "hh: " << hh << ", n: " << n << ", inv: " << inv << ", b: " << b << "rcnew: " << rcnew << std::endl;
        if(units < 0) {
          // household didn't go to the store - there's nothing to model here.
          if(retinv)
            risave2[n] = inv;
	  if(invmodel == 1) {
	    inv = dmax(inv - rcnew,0);
	  } else {
	    nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		    rcnew,revarray.begin(),nsize,badmissable.begin(),&b,&inv,
                    &nextbflag,0);
	  }
          if(debug == 1) {
            rllrun[n] = 1;
          } else if(debug == 2) {
            rllrun[n] = -1;
          }
        } else {
          indcheck(n,0,leniv,"iv. n: %d",n);

          double maxutil = -1000000;

          // compute avg predicted inclusive value, needed for importance weighting
          double rgrid1[nsize];
          for(int size=0;size<nsize;size++) {
            rgrid1[size] = iv(n,size);
          }

	  //int printobs = n <= 51141 ? 51141 : 52607; //51141;
	  int printobs = -1; //n <= 158 ? 158 : 17378;
	  //if(n >= 5163 && n <= 5163) {
	  //  printobs = n;
	  //}

          if(!myopic && beta > 0) {
            ivpred(rgrid1, ivcoef.begin() + ivblocksize*hh, cutmat.begin() + cutblocksize*hh, nd, ncut*(ivdrop[hh]==0), sptype, inttype, capj, ivout);
            if(gridmethod == 2) {
	      if(n==printobs) {
		Rcpp::Rcout << "pinds[" << n << "]: " << pinds[n] << "; ivptfilled[" << n << "]: " << ivtpfilled[pinds[n]-1] << std::endl;
		for(int ii=0;ii<nsize;ii++) {
		  Rcpp::Rcout << "rgrid1[" << ii << "]: " << rgrid1[ii] << ", ";
		}
		Rcpp::Rcout << std::endl;
		for(int ii=0;ii<nd-1;ii++) {
		  Rcpp::Rcout << "ivout[" << ii << "]: " << ivout[ii] << ", ";
		}
		Rcpp::Rcout << std::endl;
		Rcpp::Rcout << "ivcoef:" << std::endl;
		for(int ii =0;ii<ivblocksize;ii++) {
		  Rcpp::Rcout << "[" << ii << "]: " << *(ivcoef.begin() + ivblocksize*hh + ii) << "; ";
		}
		Rcpp::Rcout << std::endl;
	      }
              if(ivtpfilled[pinds[n]-1] < 0) {
                ivtpfilled[pinds[n]-1] = n;
                for(int p=0;p<nrgrid;p++) {
                  ivtprobs[p+nrgrid*(pinds[n]-1)] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,ivvari.begin()+(nd-1)*(nd-1)*hh,
                                                           detiv[hh],impdist[p]);

		  if(n==printobs && p < 10) {
		    Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p+nrgrid*(pinds[n]-1)] << ", " << impdist[p] <<
		      ", " << detiv[hh] << std::endl;
		    for(int ii=0;ii<nd-1;ii++) {
		      Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
		    }
		    Rcpp::Rcout << std::endl;
		  }
                }
              }
            }
          }

          
          if(n == printobs)
            Rcpp::Rcout << "kweights[0]: " << kweights[0] << std::endl;
          double iv0[nsize];
          for(int size=0;size<nsize;size++) {
            iv0[size] = iv(n,size);
            for(int j=0;j<capj+1;j++) {
              ivgrid[j+size*(1+capj)] = iv0[size]*((double)j)/((double)capj);
              if(n==printobs) {
                Rcpp::Rcout << "iv(" << n << "," << size << "): " << iv(n,size) <<
                  " iv0[" << size << "]: " << iv0[size] << " ivgrid[" <<
                  j+size*(1+capj) << "]: " << ivgrid[j+size*(1+capj)] << std::endl;
              }
            }
          }

          // compute value function at points that get visited - could fold this into next loop
	  if(n==printobs) {
	    Rcpp::Rcout << "inv: " << inv << ", b: " << b << ", init inv: " << initinv[hh] << ", init b: " << initb[hh] << std::endl;
	  }


          if(!myopic && beta > 0) {
            for(int size=0;size<nsize;size++) {
              for(int j=0;j<capj+1;j++) {
                if( (j == 0 && size == 0) || j > 0 ) {
		  int bprime1 = 0;
		  if(invmodel == 1) {
		    invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		    //drop[j+capj*size] = invprime > maxinv;
		    invprime = invprime > invub ? invub : invprime;
		  } else {
		    nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			    rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                            &nextbflag,0);
		    bprime1 = bprime;
		    bprime = bindex[bprime]-1;
		    indcheck(bprime,0,nb,"bprime (llijc, vf) error");
		  }
                  int iind = itoiind(invprime,istates,nistates,first20);
                  /*double dif = 0;
		  if(invprime <= istates[0]) {
		    iind = 1;
		  } else if (invprime >= istates[nistates-1]) {
		    iind = nistates;
		  } else if(invprime > istates[first20-1]) {
                    dif = istates[first20] - istates[first20-1];
                    iind = first20 + (int)((invprime-istates[first20-1])/dif);
                  } else {
                    dif = istates[1] - istates[0];
                    iind = 1 + (int)((invprime-istates[0])/dif);
		    }*/
                  int vfindx1 = iind + nistates*(bprime+nb*(pinds[n]-1));
		  int ok = iind < nistates; //check to make sure the state above this one is admissable (I think the second part should always be true)
		  if(invmodel == 2) {
		    ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
		  }
		  if(ok) {
		    indcheck(vfindx1,0,nb*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
		    //Rcpp::Rcout << "Inventory problem " << invprime << ", " << iind << std::endl;
		    //throw std::range_error("stopping");
		    if(vffilled[vfindx1] < 0) {
		      vffilled[vfindx1] = n;
		      if(gridmethod == 2) {
			vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						      nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						      gvari.begin(), detiv[hh], dg, gmean.begin(),
						      nd-1, istates, nistates, first20, bprime, nb, 0, 0, vfinterpmethod);
			if(n==printobs) {
			  //  Rcpp::Rcout << "nreptilde " << nreptilde << std::endl;
			  Rcpp::Rcout << "ivtprobs: " << *(ivtprobs+nrgrid*(pinds[n]-1)) << std::endl;
			  Rcpp::Rcout << "kweights: " << kweights[0] << std::endl;
			  Rcpp::Rcout << "vfin: " << vfsave[0] << std::endl;
			  Rcpp::Rcout << "vf output (" << n << ", " << iind << "): " << vftemp[iind+nistates*bprime] << std::endl;
			}
		      } else {
			vftemp[vfindx1] = vfavgijchh(ivout, iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						     nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						     gvari.begin(), detiv[hh], dg, gmean.begin(),
						     nd-1, istates, nistates, first20, bprime, nb, 0);
		      }
		      if(n == printobs) {
			Rcpp::Rcout << " 1 hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
			  " i': " << invprime << "b': " << bprime << " vf: " << vftemp[vfindx1] << "; "
				    << gridmethod << std::endl;
		      }

		    }

		  }

                  // uncomment if not using price states
                  //vfindx1 = iind + nistates*bprime;
                  //if(vffilled[vfindx1] != n) {
		  
		  vfindx1 = iind-1 + nistates*(bprime+nb*(pinds[n]-1));
		  indcheck(vfindx1,0,nb*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
		  //vfindx1 = iind-1 + nistates*bprime;
		  //if(vffilled[vfindx1] != n) {
		  if(vffilled[vfindx1] < 0) {
		    vffilled[vfindx1] = n;

		    if(gridmethod == 2) {
		      vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						    nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						    gvari.begin(), detiv[hh], dg, gmean.begin(),
						    nd-1, istates, nistates, first20, bprime, nb, 0, 0, vfinterpmethod);
		    } else {
		      vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						   nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						   gvari.begin(), detiv[hh], dg, gmean.begin(),
						   nd-1, istates, nistates, first20, bprime, nb, 0);
		    }

		    if(n == printobs) {
		      Rcpp::Rcout << " 2 hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
			  " i': " << invprime  << "b': " << bprime << " vf: " << vftemp[vfindx1] << "; "
				  << gridmethod << std::endl;
		    }

		  }

                }
              }
            }
          }

          int first = 1;

	  int vfoffset = nistates*nb*(pinds[n]-1);

          for(int size=0;size<nsize;size++) {
            for(int j=0;j<capj+1;j++) {
              int indx = j+(1+capj)*size;
              drop[indx] = 1;
              if( (j == 0 && size == 0) || j > 0 ) {
                drop[indx] = 0;
                // change here - we always allow a choice even if above inventory bound
		if(invmodel == 1) {
		  invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		  invprime = invprime > invub ? invub : invprime;
		  ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize[size], rcnew, gamma, nu, ccost, omega, cmodel, 0);
		} else {
		  nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                          &nextbflag,n==printobs && j==1 && size == 1);
		  bchoice[indx] = bprime;
		  ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				     nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), n==printobs);
		  bprime = bindex[bprime]-1;
		  indcheck(bprime,0,nb,"bprime (llijc, util) error");
		  ichoice[indx] = invprime;
		}
                if(!myopic && beta > 0) {
                  //Rcpp::Rcout << "n " << n << "; ";
                  vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20, vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
                } else {
                  vfhat[indx] = 0.0;
                }
                vfll[n+obsend*(indx)] = vfhat[indx];
                utll[n+obsend*(indx)] = ut[indx];
                if(first) {
                  first=0;
                  maxutil = ut[indx]+vfhat[indx];
                } else {
                  maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);
                }
                if(n==printobs) {
                  Rcpp::Rcout << " ivgrid[" << indx << "]: " << ivgrid[indx] << " ut[" << indx << "]: " << ut[indx] << " vfhat[" << indx << "]: " << vfhat[indx] << " flag: " << nextbflag << std::endl;
                }
              }
            }
          }

          double s1 = 0.0; //exp(ut[0]+vfhat[0]-maxutil);

          for(int size=0;size<nsize;size++) {
            for(int j=0;j<capj+1;j++) {
              if(!drop[j+(1+capj)*size]) {
                s1 += exp(ut[j+(1+capj)*size]+vfhat[j+(1+capj)*size]-maxutil);
              }
            }
          }

          if(debug == 2) {
	    if(n==printobs) {
	      Rcpp::Rcout << "Implied Choice Probabilities:" << std::endl;
	    }
            for(int size=0;size<nsize;size++) {
              for(int j=0;j<capj+1;j++) {
                if(!drop[j+(1+capj)*size]) {
                  rllrun[n+nobs*(j+(1+capj)*size)] = exp(ut[j+(1+capj)*size]+vfhat[j+(1+capj)*size]-maxutil)/s1;
		  if(n==printobs) {
		    Rcpp::Rcout << "size = " << size+1 << ", j = " << j << ": " << rllrun[n+nobs*(j+(1+capj)*size)] << std::endl;
		  }
                }
              }
            }
          }

          double obspack = 0;
          int obschoice = 0;

          if(brindex[n] > 0) {
            obspack = packsize[brsize[nbrand + (int)brindex[n] -1]-1];
            obschoice = (int)tunits[n] + (1+capj)*(brsize[nbrand + (int)brindex[n] -1]-1);
          }

          indcheck(obschoice,0,(1+capj)*nsize,"ut. tunits: %d",((int)tunits[n]));

          // this line shouldn't be necessary since we remove inv above the bound
          //if( dmax(inv + obspack*tunits[n] - rcnew,0) > invub )
          // this disallows purchases the individual could never carry
          if(obspack*tunits[n] > invub)
            llrun -= 1000000;

          double cprob = exp(ut[obschoice] + vfhat[obschoice] - maxutil)/s1;

          llrun += log(cprob);

          if(n==printobs)
            Rcpp::Rcout << "log likelihood[" << n << "]: " << log(cprob) << std::endl;

          //}

          if(debug == 1) {
            //if( dmax(inv + obspack*tunits[n] - rcnew,0) > invub ) {
            if(obspack*tunits[n] > invub) {
              rllrun[n] = -1000000;
            } else {
              rllrun[n] = exp(ut[obschoice]+vfhat[obschoice]-maxutil)/s1;
            }
          }
          if(retinv)
            risave2[n] = inv;
	  if(invmodel==1) {
	    inv = dmax(inv + obspack*tunits[n] - rcnew,0);
	    inv = inv > invub ? invub : inv;
	  } else {
	    inv = ichoice[obschoice];
	    b = bchoice[obschoice];
	  }
        }
      }

      llhh[hh] = llrun;

    }

  }

};


// compute (dynamic) part of log likelihood
// I think vf is unneccessary here

void llijc(const Rcpp::NumericMatrix &co, const Rcpp::NumericVector &cosave,
           const Rcpp::NumericVector &rgridsave, const Rcpp::IntegerVector &indexes,
           const Rcpp::NumericVector &tunits,
           const Rcpp::NumericVector &panid, const Rcpp::IntegerVector &brindex,
           const Rcpp::IntegerVector &hhinds, const Rcpp::IntegerVector &brsize,
           const Rcpp::IntegerVector &obshhinds, const Rcpp::NumericMatrix &iv,
           const Rcpp::NumericVector &vfsave, const Rcpp::NumericVector &cdraws,
           const Rcpp::NumericVector &ic, const Rcpp::NumericVector &iiv,
           const Rcpp::NumericVector &ldraws, const Rcpp::NumericVector &packsize,
           const Rcpp::NumericVector &initx, const Rcpp::IntegerVector &initb,
	   const Rcpp::IntegerVector &inits, const Rcpp::IntegerVector &ngrid,
           const Rcpp::NumericVector &cutmat, const Rcpp::NumericVector &ivcoef,
           const Rcpp::NumericVector &ivvari, const Rcpp::NumericVector &detiv,
           const Rcpp::NumericVector &pgrid, const Rcpp::NumericVector &qpts,
           const Rcpp::NumericVector &gvari, const Rcpp::NumericVector &gmean,
           const Rcpp::NumericMatrix &rgrid, const Rcpp::NumericVector &impdist,
           const Rcpp::NumericMatrix &bwmatrix, const Rcpp::IntegerVector &pinds,
	   const Rcpp::IntegerVector &bstates, const Rcpp::IntegerVector &revarray,
	   const Rcpp::IntegerVector &badmissable, const Rcpp::IntegerVector &bindex,
	   const Rcpp::IntegerVector &ivdrop,
           const llinputs llinfo, const int rep, const int gridmethod,
           const int vfinterpmethod, const int maxpind,
           Rcpp::NumericVector &vfll, Rcpp::NumericVector &utll,
           Rcpp::NumericVector &rllrun,
           Rcpp::NumericVector &risave1, Rcpp::NumericVector &risave2,
           Rcpp::NumericVector &llhh) {

  int nobs = llinfo.nobs;

  int nbrand = llinfo.nbrand;

  int nsize = llinfo.nsize;

  int inttype = llinfo.inttype;

  int nsave = llinfo.nsave;

  int ntilde = llinfo.ntilde;

  int nrgrid = llinfo.nrgrid;

  int capj = llinfo.capj;

  int nhhs = llinfo.nhhs;

  int ninitt = llinfo.ninitt;

  int ncutiv = llinfo.ncutiv;

  int retinv = llinfo.retinv;

  double maxinv = llinfo.maxinv;

  int myopic = llinfo.myopic;

  int cmodel = llinfo.cmodel;

  int debug = llinfo.debug;

  int sptype = llinfo.sptype;

  int first20 = llinfo.first20;

  int hinvbound = llinfo.hinvbound;

  int nco = llinfo.nco;

  int hmodel = llinfo.hmodel;

  int usecdraws = llinfo.usecdraws;

  int initvf = llinfo.initvf;

  int lomega = llinfo.lomega;

  int linitx = initx.size();

  int nd = llinfo.nd;
  int ncut = llinfo.ncut;
  double dg = llinfo.dg;

  int invmodel = llinfo.invmodel;
  int maxbottles = llinfo.maxbottles;

  int nistates = ngrid[nd-1];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }

  int ivblocksize;
  if(sptype == 1) {
    ivblocksize = (2*ncutiv+2)*nsize*nsize;
  } else {
    if(ncutiv == 0) {
      ivblocksize = (nsize+1)*nsize;
    } else {
      ivblocksize = (ncutiv+4)*nsize*nsize;
    }
  }

  int cutblocksize = (ncutiv+1)*nsize;

  int nrep, nreptilde;

  int rep1 = rep+1; //conform to the R code

  if(rep1 <= nsave) {
    nrep = rep1-1;
  } else {
    nrep = nsave;
  }

  if(nrep < ntilde) {
    nreptilde = rep1-1;
  } else {
    nreptilde = ntilde;
  }

  if(initvf) {
    // we're supplying vf and parameter draw starting points
    nrep = nsave;
    nreptilde = ntilde;
    rep1 = nrep+1;
  }

  double initinv[nhhs];
  int initbstate[nhhs];

  //int dummyb[1];

  //dummyb[0] = 1;

  int nb = llinfo.nbsmall;

  double clb = 0, cub = 0;

  double inv = 0;

  int b=0;

  double invub = maxinv;

  //int lomega = 2;  // update this later
  double omega[lomega];

  int nextbflag;

  for(int n=hhinds[0]-1;n<hhinds[1];n++) {
    inv = 0;
    b = 0;
    for(int j=0;j<lomega;j++) {
      omega[j] = co[nbrand+9+j+nco*n];
    }
    clb = co[nbrand+nco*n];
    cub = co[nbrand+1+nco*n];

    for(int t=0;t<ninitt;t++) {
      double chosenq = 0;
      double rcnew;
      if(usecdraws) {
        rcnew = ic[n+nhhs*t];
      } else {
        rcnew = Rcpp::runif(1)[0];
      }
      if(rcnew < 0.5) {
        rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
      } else {
        rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
      }
      indcheck(t+ninitt*n,0,linitx,"initx: n: %d; t %d",n,t);
      if(invmodel == 1) {
	chosenq = initx[t+ninitt*n];  // since this vector is ordered by id, then week, this should be the right indexing
	invub = hinvbound ? omega[0] : maxinv;

	if(retinv)
	  risave1[n+nhhs*t] = inv;
	inv = dmax(inv + chosenq - rcnew,0);
	inv = inv > invub ? invub : inv;
      } else {
	nextinv(initb[t+ninitt*n],inits[t+ninitt*n]-1,b,inv,bstates.begin(),
		nb,maxbottles,packsize.begin(),
		rcnew,revarray.begin(),nsize,badmissable.begin(),&b,&inv,
                &nextbflag,0);
      }
      //if(n == 127 && invmodel == 2) {
      //	Rcpp::Rcout << "init for " << n << ", " << t << ": inv: " << inv << ", b " << b <<
      //	  ", j: " << initb[t+ninitt*n] << ", s: " << inits[t+ninitt*n] << std::endl;
      //}
    }
    initinv[n] = inv;
    if(invmodel == 2) {
      initbstate[n] = b;
    }
  }
  //Rcpp::Rcout << "aaa" << std::endl;
  // compute likelihood at data

  llijcpar llpar(co,cosave,rgridsave,indexes,tunits,panid,brindex,brsize,
                 obshhinds,iv,vfsave,cdraws,ic,iiv,ldraws,packsize,initx,
                 ngrid,cutmat,ivcoef,ivvari,detiv,pgrid,qpts,gvari,gmean,
                 rgrid,impdist,bwmatrix,pinds,bstates,revarray,badmissable,
		 bindex,ivdrop,
                 nobs,nbrand,nsize,inttype,nsave,ntilde,nrgrid,capj,nhhs,
                 ninitt,ncutiv,retinv,maxinv,myopic,cmodel,debug,sptype,
                 first20,hinvbound,nco,linitx,nd,ncut,dg,nistates,istates,
                 ivblocksize,cutblocksize,nrep,nreptilde,initinv,initbstate,nb,
                 rep1,gridmethod,vfinterpmethod,maxpind,hmodel,invmodel,
		 maxbottles,lomega,
		 vfll,utll,rllrun,risave1,risave2,llhh);
  //Rcpp::Rcout << "object created" << std::endl;
  parallelFor(hhinds[0]-1,hhinds[1],llpar);

  //llpar(hhinds[0]-1,hhinds[1]);

  //END_RCPP

}

// density for triangular distribution

double tdist(const double x, const double * bds) {
  double midpoint = 0.5*(bds[1]+bds[2]);
  return( ((x < midpoint)*(x-bds[1])+(x >= midpoint)*(bds[2]-x))*4/((bds[2]-bds[1])*(bds[2]-bds[1])) );
}




// parallel vf update

/// new args - brsize, impdist, gridmethod - add to other calls

struct vfupdatepar : public RcPar::Worker
{

  // inputs
  const RcPar::RVector<int> hhinds;
  const RcPar::RVector<int> bstates;
  const RcPar::RVector<double> pgrid;
  const RcPar::RVector<int> ngrid;
  const RcPar::RMatrix<double> co;
  const RcPar::RVector<double> vfsave;
  const RcPar::RVector<double> rgrid;
  const RcPar::RVector<double> cosave;
  const RcPar::RVector<double> rgridsave;
  const RcPar::RVector<double> packsize;
  const RcPar::RVector<double> cutmat;
  const RcPar::RVector<double> ivcoef;
  const RcPar::RVector<double> ivvari;
  const RcPar::RVector<double> detiv;
  const RcPar::RVector<double> gvari;
  const RcPar::RVector<double> gmean;
  const RcPar::RVector<int> indexes;
  const RcPar::RVector<int> brsize;
  const RcPar::RVector<double> impdist;
  const RcPar::RMatrix<double> bwmatrix;
  const RcPar::RVector<int> revarray;
  const RcPar::RVector<int> badmissable;
  const RcPar::RVector<int> bindex;
  const RcPar::RVector<int> ivdrop;
  const int capj;
  const int nbrand;
  const int cmodel;
  const int ncutiv;
  const int ncut;
  const int sptype;
  const int inttype;
  const int first20;
  const int nsave;
  const int ntilde;
  const int ncrate;
  const double dg;
  const int hinvbound;
  const int lrgrid;
  const int livcoef;
  const int lcutmat;
  const int lvf;
  const int ldetiv;
  const int nd;
  const int nistates;
  const double * istates;
  const int nsize;
  const int nrgrid;
  const R_len_t nb;
  const int nbsmall;
  const int nco;
  const int nhhs;
  const int * dummyb;
  const int ivblocksize;
  const int cutblocksize;
  const int nrep;
  const int nreptilde;
  const int rep;
  const int livvari;
  const int niter;
  const int gridmethod;
  const int vfinterpmethod;
  const int hmodel;
  const int invmodel;
  const int maxbottles;
  const int lomega;

  //output
  RcPar::RVector<double> rvfout;

  // initialize with source and destination
  vfupdatepar(const Rcpp::IntegerVector hhinds, const Rcpp::IntegerVector bstates,
              const Rcpp::NumericVector pgrid, const Rcpp::IntegerVector ngrid,
              const Rcpp::NumericMatrix co, const Rcpp::NumericVector vfsave,
              const Rcpp::NumericVector rgrid, const Rcpp::NumericVector cosave,
              const Rcpp::NumericVector rgridsave, const Rcpp::NumericVector packsize,
              const Rcpp::NumericVector cutmat, const Rcpp::NumericVector ivcoef,
              const Rcpp::NumericVector ivvari, const Rcpp::NumericVector detiv,
              const Rcpp::NumericVector gvari, const Rcpp::NumericVector gmean,
              const Rcpp::IntegerVector indexes, const Rcpp::IntegerVector brsize,
              const Rcpp::NumericVector impdist, const Rcpp::NumericMatrix bwmatrix,
	      const Rcpp::IntegerVector revarray,
	      const Rcpp::IntegerVector badmissable, const Rcpp::IntegerVector bindex,
	      const Rcpp::IntegerVector ivdrop,
              const int capj, const int nbrand, const int cmodel, const int ncutiv,
              const int ncut, const int sptype, const int inttype,
              const int first20, const int nsave, const int ntilde,
              const int ncrate, const double dg, const int hinvbound,
              const int lrgrid, const int livcoef, const int lcutmat,
              const int lvf, const int ldetiv, const int nd,
              const int nistates, const double * istates, const int nsize,
              const int nrgrid, const int nb, const int nbsmall, const int nco,
              const int nhhs, const int * dummyb, const int ivblocksize,
              const int cutblocksize, const int nrep, const int nreptilde,
              const int rep, const int livvari, const int niter,
              const int gridmethod, const int vfinterpmethod, const int hmodel,
	      const int invmodel, const int maxbottles, const int lomega,
              Rcpp::NumericVector rvfout) :
  hhinds(hhinds), bstates(bstates), pgrid(pgrid), ngrid(ngrid), co(co),
    vfsave(vfsave), rgrid(rgrid), cosave(cosave), rgridsave(rgridsave),
    packsize(packsize), cutmat(cutmat), ivcoef(ivcoef), ivvari(ivvari),
    detiv(detiv), gvari(gvari), gmean(gmean), indexes(indexes),
    brsize(brsize), impdist(impdist), bwmatrix(bwmatrix), revarray(revarray),
    badmissable(badmissable), bindex(bindex), ivdrop(ivdrop),
    capj(capj), nbrand(nbrand), cmodel(cmodel), ncutiv(ncutiv),
    ncut(ncut), sptype(sptype), inttype(inttype), first20(first20),
    nsave(nsave), ntilde(ntilde), ncrate(ncrate), dg(dg),
    hinvbound(hinvbound), lrgrid(lrgrid), livcoef(livcoef),
    lcutmat(lcutmat), lvf(lvf), ldetiv(ldetiv), nd(nd), nistates(nistates),
    istates(istates), nsize(nsize), nrgrid(nrgrid), nb(nb), nbsmall(nbsmall), 
    nco(nco),
    nhhs(nhhs), dummyb(dummyb), ivblocksize(ivblocksize),
    cutblocksize(cutblocksize), nrep(nrep), nreptilde(nreptilde),
    rep(rep), livvari(livvari), niter(niter), gridmethod(gridmethod),
    vfinterpmethod(vfinterpmethod), hmodel(hmodel), invmodel(invmodel),
    maxbottles(maxbottles), lomega(lomega), rvfout(rvfout) {}

  //compute value function for specified households

  void operator()(std::size_t begin, std::size_t end) {

    int begin1 = (int)begin;
    int end1 = (int)end;

    double cratevec[ncrate];
    double cprobs[ncrate];
    double endpoints[ncrate+1];

    R_len_t c1;

    double ivgrid[(1+capj)*nsize];
    double vftemp[nbsmall*nistates];

    double ivout[(nd-1)*nrgrid];
    double ivtprobs[nrgrid*nrgrid];  // only needed for gridmethod == 2
    int coorder[nsave];

    int drop[nsize*(capj+1)];
    double uu[nsize*(capj+1)];

    double rgridhh[(nd-1)*nrgrid]; // rgrid for a given household

    int maxit = gridmethod == 2 ? niter : 1;

    double kweights[nsave];

    int hhbegin = begin1;
    int hhend = end1;
    //int pbegin = 0;
    //int pend = nrgrid;
    int ibegin = 0;
    int iend = nistates;
    
    if(hmodel) {
      hhbegin=0;
      hhend=1;
      //pbegin=begin1;
      //pend=end1;
      ibegin=begin1;
      iend=end1;
    }

    for(int hh=hhbegin;hh<hhend;hh++) {

      if(hh == -1) {
	for(int i=0;i<nco;i++) {
	  Rcpp::Rcout << " p[" << i << "]: " << co[i+nco*hh];
	}
	Rcpp::Rcout << std::endl;
      }

      // get household coefficient
      // double xi[nbrand];
      double clb = co[nbrand+nco*hh];
      double cub = co[nbrand+1+nco*hh];
      //double alpha = co[nbrand+2+nco*hh];
      double gamma = co[nbrand+3+nco*hh];
      double nu = co[nbrand+4+nco*hh];
      double beta = co[nbrand+5+nco*hh];
      double uf = co[nbrand+6+nco*hh];
      double ccost = co[nbrand+7+nco*hh];
      double stvisitprob = co[nbrand+8+nco*hh];

      int vffilled[nbsmall*nistates];
      for(int b=0;b<nbsmall;b++) {
        for(int i=0;i<nistates;i++) {
          vffilled[i+nistates*b] = 0;
        }
      }

      //for(int j=0;j<nbrand;j++) {
      //  xi[j] = co[j+nco*hh];
      //}

      //int lomega = nco-9-nbrand;
      double omega[lomega];
      for(int j=0;j<lomega;j++) {
        omega[j] = co[9+j+nbrand+nco*hh];
      }

      if(beta == 0) {
        memset(rvfout.begin()+(hh-(hhinds[0]-1))*nrgrid*nbsmall*nistates,0,sizeof(double)*nrgrid*nbsmall*nistates);
      } else {
        double invub = hinvbound ? omega[0] : pgrid[nd-1+nd*(nistates-1)];

        // compute consumption rate draws

        for(int i=0;i<ncrate+1;i++) {
          endpoints[i] = ((double)i)/((double)ncrate)*(cub-clb)+clb;
        }

        for(int i=0;i<ncrate;i++) {
          cratevec[i] = 0.5*(endpoints[i]+endpoints[i+1]);
        }

        if(ncrate == 1) {
          cprobs[0] = 1;
        } else {
          double midpoint = 0.5*(clb+cub);
          for(int i=0;i<ncrate;i++) {
            if( endpoints[i] < midpoint && endpoints[i+1] > midpoint ) {
              double a1 = tdist(endpoints[i],co.begin()+nbrand+nco*hh);
              double a2 = tdist(midpoint,co.begin()+nbrand+nco*hh);
              cprobs[i] = dmin(a1,a2)*(midpoint-endpoints[i])+0.5*(midpoint-endpoints[i])*fabs(a1-a2);
              a1 = tdist(midpoint,co.begin()+nbrand+nco*hh);
              a2 = tdist(endpoints[i+1],co.begin()+nbrand+nco*hh);
              cprobs[i] += dmin(a1,a2)*(endpoints[i+1]-midpoint)+0.5*(endpoints[i+1]-midpoint)*fabs(a1-a2);
            } else {
              double a1 = tdist(endpoints[i],co.begin()+nbrand+nco*hh);
              double a2 = tdist(endpoints[i+1],co.begin()+nbrand+nco*hh);
              cprobs[i] = dmin(a1,a2)*(endpoints[i+1]-endpoints[i])+0.5*(endpoints[i+1]-endpoints[i])*fabs(a1-a2);
            }
          }
        }

        //int indxprint = -1;

        int vfind = hh*nrgrid*nbsmall*nistates*nsave;
        indcheck(nco*hh,0,nco*nhhs,"co offset bug hh: %d",hh);
        indcheck(nsave*nco*hh,0,nsave*nco*nhhs,"cosave offset bug hh: %d",hh);
        if(vfinterpmethod == 1) {
          getclosestind(co.begin() + nco*hh, cosave.begin() + nsave*nco*hh, indexes.begin(), nco, nrep, nsave, rep, coorder);
        } else {
          getkweights(co.begin() + nco*hh, cosave.begin() + nsave*nco*hh,
                      indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
                      rep, coorder, kweights,0);
        }

        if(gridmethod == 1) {
          // average inclusive value tomorrow.
          for(int p=0;p<nrgrid;p++) {
            indcheck((nd-1)*p,0,lrgrid,"rgrid offset bug p: %d",p);
            indcheck(ivblocksize*hh,0,livcoef,"ivcoef offset bug p: %d",p);
            indcheck(cutblocksize*hh,0,lcutmat,"cutmat offset bug p: %d",p);
            ivpred(rgrid.begin()+(nd-1)*p, ivcoef.begin() + ivblocksize*hh, cutmat.begin() + cutblocksize*hh, nd, ncut*(ivdrop[hh]==0), sptype, inttype, capj, ivout+(nd-1)*p);
          }
        } else if(gridmethod == 2) {
          // note we have to compute the importance grid before getting
          // the transition probabilities
          if(hh == -1) {
            int ofset = ivblocksize*hh;
            //for(int ii=0;ii<nd-1;ii++) {
            //  Rcpp::Rcout << "ivcoef[" << ii << "]: " << ivcoef[ofset+0+nd*ii];
            //  for(int kk=0;kk<nd-1;kk++) {
            //    Rcpp::Rcout << ", " << ivcoef[ofset+1+kk+nd*ii];
            //  }
            //  Rcpp::Rcout << std::endl;
	    //}
	    if(ncutiv==0 || ivdrop[hh]) {
	      for(int ii=0;ii<nd-1;ii++) {
		for(int kk=0;kk<nd;kk++) {
		  Rcpp::Rcout << " ivcoef[" << kk + nd*ii << "]: " << ivcoef[ofset+kk + nd*ii];
		}
		Rcpp::Rcout << std::endl;
	      }
	    } else {
	      for(int ii=0;ii<nd-1;ii++) {
		for(int kk=0;kk<nd-1;kk++) {
		  Rcpp::Rcout << " ivcoef[" << (ncut+4)*kk + (nd-1)*(ncut+4)*ii << "]: " << ivcoef[ofset+(ncut+4)*kk + (nd-1)*(ncut+4)*ii];
		}
		Rcpp::Rcout << std::endl;
	      }
	    }
	    ofset = cutblocksize*hh;
	    for(int ii=0;ii<cutblocksize;ii++) {
	      Rcpp::Rcout << "cutmat[" << ii << "]: " << cutmat[ofset+ii];
	      //for(int kk=0;kk<nsize;kk++) {
	      //  Rcpp::Rcout << ", " << ivcoef[ofset+1+kk+nsize*ii];
	      //}
	      Rcpp::Rcout << std::endl;
	    }
	    int p1 = 84;
	    for(int ii=0;ii<nd-1;ii++) {
	      Rcpp::Rcout << "rgrid[" << ii+(nd-1)*p1 << "]: " << rgrid[ii+(nd-1)*p1] << std::endl;
	    }
          }
          for(int p=0;p<nrgrid;p++) {
            indcheck(nbrand*p,0,lrgrid,"rgrid offset bug p: %d",p);
            indcheck(ivblocksize*hh,0,livcoef,"ivcoef offset bug p: %d",p);
            indcheck(cutblocksize*hh,0,lcutmat,"cutmat offset bug p: %d",p);

            computeiv(co.begin()+nco*hh,rgrid.begin()+nbrand*p,brsize.begin(),
                      nsize,nbrand,rgridhh+(nd-1)*p);

            ivpred(rgridhh+(nd-1)*p, ivcoef.begin() + ivblocksize*hh, cutmat.begin() + cutblocksize*hh, nd, ncut*(ivdrop[hh]==0), sptype, inttype, capj, ivout+(nd-1)*p);
          }
          for(int p=0;p<nrgrid;p++) {
            for(int p1=0;p1<nrgrid;p1++) {
              ivtprobs[p1+nrgrid*p] = impfn2(ivout+(nd-1)*p,rgridhh+(nd-1)*p1,nd-1,
                                             ivvari.begin()+(nd-1)*(nd-1)*hh,
                                             detiv[hh],impdist[p1]);
              /*if(p==84 && p1 == 5 && hh == -1) {
                Rcpp::Rcout << "p, p1: " << p << ", " << p1 << ". rgridhh: ";
                for(int ii=0;ii<nd-1;ii++)
                  Rcpp::Rcout << rgridhh[ii+(nd-1)*p1] << " ";
                Rcpp::Rcout << ". ivout: ";
                for(int ii=0;ii<nd-1;ii++)
                  Rcpp::Rcout << ivout[ii+(nd-1)*p1] << " ";
                Rcpp::Rcout << ". ivtprob: " << ivtprobs[p1+nrgrid*p] <<
                  ". impdist: " << impdist[p1] << "." << std::endl;
		  }*/
            }
          }
        }

        for(int iter=0;iter<niter;iter++) {
          for(int p=0;p<nrgrid;p++) {
            // we have to fill in the value function here since it'll get used at all the i states
            // when computing choices.

	    // note - update this to not evaluate at redundant states

            indcheck(vfind,0,lvf,"vf offset bug hh %d",hh);
            indcheck((nd-1)*(nd-1)*hh,0,livvari,"ivvari offset bug hh %d",hh);
            indcheck(hh,0,ldetiv,"detiv offset bug hh %d",hh);
            if(iter == 0) {
              if(gridmethod == 1) {
                for(int b=0;b < nb;b++) {
		  if(badmissable[b]) {
		    int bsmall = bindex[b]-1;
		    for(int i=0;i<nistates;i++) {
		      if(invmodel == 1 || (b==0 && i==0) || (b > 0 && pgrid[nd-1+nd*i] <= packsize[bstates[0+maxbottles*b]-1])) {
			vftemp[i+nistates*bsmall] = vfavgijchh(ivout+(nd-1)*p, i, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
							  nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
							  gvari.begin(), detiv[hh], dg, gmean.begin(),
							  nd-1, istates, nistates, first20, bsmall, nbsmall, 0);
		      }
		    }
                  }
                }
              } else if (gridmethod == 2) {
                for(int b=0;b < nb;b++) {
		  if(badmissable[b]) {
		    int bsmall = bindex[b]-1;
		    for(int i=0;i<nistates;i++) {
		      if(invmodel == 1 || (b==0 && i==0) || (b > 0 && pgrid[nd-1+nd*i] <= packsize[bstates[0+maxbottles*b]-1])) {
			vftemp[i+nistates*bsmall] = vfavgijchh2(ivtprobs+nrgrid*p, i, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
							   nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
							   gvari.begin(), detiv[hh], dg, gmean.begin(),
							   nd-1, istates, nistates, first20, bsmall, nbsmall, hh<0 && i == 0 && p == 0, iter, vfinterpmethod);
		      }
		    }
                  }
                }
              }
            } else {
              for(int b=0;b < nb;b++) {
		if(badmissable[b]) {
		  int bsmall = bindex[b]-1;
		  for(int i=0;i<nistates;i++) {
		    if(invmodel == 1 || (b==0 && i==0) || (b > 0 && pgrid[nd-1+nd*i] <= packsize[bstates[0+maxbottles*b]-1])) {
		      vftemp[i+nistates*bsmall] = vfavgijchh2(ivtprobs+nrgrid*p, i, rvfout.begin()+(hh-(hhinds[0]-1))*nrgrid*nb*nistates,
							 rgridsave.begin(), indexes.begin(), coorder, kweights,
							 nreptilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
							 gvari.begin(), detiv[hh], dg, gmean.begin(),
							 nd-1, istates, nistates, first20, bsmall, nbsmall, 0, iter, vfinterpmethod);
		    }
		  }
                }
              }
            }

            if(gridmethod == 1) {
              for(int i=0;i<nd-1;i++) {
                for(int j=0;j<capj+1;j++) {
                  ivgrid[j+i*(1+capj)] = rgrid[i+(nd-1)*p]*((double)j)/((double)capj);
                }
              }
            } else if(gridmethod == 2) {
              for(int i=0;i<nd-1;i++) {
                for(int j=0;j<capj+1;j++) {
                  ivgrid[j+i*(1+capj)] = rgridhh[i+(nd-1)*p]*((double)j)/((double)capj);
                }
              }
            }

            for(int b=0;b < nb;b++) {
	      if(badmissable[b]) {
		int bsmall = bindex[b]-1;
		for(int i=0;i<nistates;i++) {

		  // only evaluate at admissable states

		  if(invmodel == 1 || (b==0 && i==0) || (b > 0 && pgrid[nd-1+nd*i] <= packsize[bstates[0+maxbottles*b]-1])) {

		    double avgv = 0;
		    for(c1=0;c1 < ncrate;c1++) {
		      for(int ss=0;ss < 2;ss++) {
			double sumeu=0;
			int maxj;
			if(ss==0) {
			  maxj = 1;
			} else {
			  maxj = capj+1;
			}

			double maxutil = -1000000;

                        int nextbflag=0;
                        
			for(int j=0;j < maxj;j++) {
			  for(int sz=0;sz<nsize;sz++) {
			    drop[j+(1+capj)*sz] = 1;
			    if( (j == 0 && sz == 0) || j > 0 ) {
			      drop[j+(1+capj)*sz] = 0;
			      double invprime = 0;
			      double vfprime = 0;
			      int nextb = 0;
			      if(invmodel==1) {
				invprime = dmax(pgrid[nd-1+nd*i]+packsize[sz]*((double)j)-cratevec[c1],0.0);
				invprime = invprime > invub ? invub : invprime;
				double vfprime = vfinterplinfast(invprime, vftemp, b, 0, istates, nistates, first20, vffilled,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
				uu[j+(1+capj)*sz] = (utiliv(j,dummyb,ivgrid[j+(1+capj)*sz],pgrid[nd-1+nd*i], packsize[0], cratevec[c1], gamma, nu, ccost, omega, cmodel, 0) + beta*vfprime)/uf;
				/*if(p+nrgrid*(i+nistates*bsmall)==584 && hh==117) {
				  Rcpp::Rcout << "j,s: " << j << ", " << sz << "; vfprime: " << vfprime << "; u: " << utiliv(j,dummyb,ivgrid[j+(1+capj)*sz],pgrid[nd-1+nd*i], packsize[0], cratevec[c1], gamma, nu, ccost, omega, cmodel, 0) << std::endl;
				  }*/
			      } else {
				nextinv(j,sz,b,pgrid[nd-1+nd*i],bstates.begin(),nb,maxbottles,packsize.begin(),
					cratevec[c1],revarray.begin(),nsize,badmissable.begin(),&nextb,&invprime,
                                        &nextbflag,0);
				int nextb1 = nextb;
				uu[j+(1+capj)*sz] = utiliv2(j,sz,ivgrid[j+(1+capj)*sz],cratevec[c1], gamma, nu, ccost, omega, lomega, cmodel, pgrid[nd-1+nd*i], b, nextb, nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
				nextb = bindex[nextb]-1;
				indcheck(nextb,0,nb,"nextb error");
			    
				uu[j+(1+capj)*sz] = (uu[j+(1+capj)*sz] + beta*vfinterplinfast(invprime, vftemp, nextb, nextb1, istates, nistates, first20, vffilled,invmodel,bstates.begin(),nb,maxbottles,packsize.begin()))/uf;


				
			      /*if(j==0 && ss==0) {
				int currentbb[maxbottles];
				int nextbb[maxbottles];
				Rcpp::Rcout << std::endl;
				Rcpp::Rcout << "b: " << b << "; i: " << i << "; j: " << j << "; size: "
				<< sz << "; c1: " << cratevec[c1] << std::endl;
				Rcpp::Rcout << "************************************" << std::endl;
				for(int k=0;k<maxbottles;k++) {
				Rcpp::Rcout << "currentb[";
				currentbb[k] = bstates[k+maxbottles*b];
				nextbb[k] = bstates[k+maxbottles*nextb];
				Rcpp::Rcout << k << "]: " << currentbb[k] << "; nextb[" << k << "]: " <<
				nextbb[k] << std::endl;
				}
				Rcpp::Rcout << "current i: " << pgrid[nd-1+nd*i] << "; next i: " <<
				invprime << std::endl;
			      
				}*/
			      }
                          
                          
			      maxutil = dmax(maxutil,uu[j+(1+capj)*sz]);

			    }
			  }
			}
			for(int j=0;j < maxj;j++) {
			  for(int sz=0;sz<nsize;sz++) {
			    if( (j == 0 && sz == 0) || j > 0 ) {
			      if(!drop[j+(1+capj)*sz]) {
				uu[j+(1+capj)*sz] -= maxutil;
				sumeu += exp(uu[j+(1+capj)*sz]);
			      }
			    }
			  }
			}
			if(ss==0) {
			  avgv += uf*(log(sumeu)+maxutil)*cprobs[c1]*(1.0-stvisitprob);
			} else {
			  avgv += uf*(log(sumeu)+maxutil)*cprobs[c1]*stvisitprob;
			}
			

		      }
		    }


		    int vfind1 = (hh-(hhinds[0]-1))*nrgrid*nbsmall*nistates;
		    indcheck(vfind1 + p+nrgrid*(i+nistates*bsmall),0,(hhinds[1]-hhinds[0]+1)*nrgrid*nbsmall*nistates,"i: %d; p: %d",i,p);

		    rvfout[vfind1 + p+nrgrid*(i+nistates*bsmall)] = avgv;
		    /*if(p+nrgrid*(i+nistates*bsmall)==584 && hh==117) {
		      Rcpp::Rcout << "rvf: " << avgv << std::endl;
		    }*/
		  }

		}

	      }
	    }
	    //throw std::range_error("aaa");
          }
        }
      }

    }

  }


};

// update value function

// note that when we do updating or evaluation of likelihood - we only need ivdep at the current draw since it only goes in current
// utility
// we also only use cutmat, ivcoef, etc, at past draws since they are used in the likelihood

void vfupdateijc(const Rcpp::IntegerVector &hhinds, const Rcpp::IntegerVector &bstates,
                 const Rcpp::NumericVector &pgrid, const Rcpp::IntegerVector &ngrid,
                 const Rcpp::NumericMatrix &co, const Rcpp::NumericVector &vfsave,
                 const Rcpp::NumericVector &rgrid, const Rcpp::NumericVector &cosave,
                 const Rcpp::NumericVector &rgridsave, const Rcpp::NumericVector &packsize,
                 const Rcpp::NumericVector &cutmat, const Rcpp::NumericVector &ivcoef,
                 const Rcpp::NumericVector &ivvari, const Rcpp::NumericVector &detiv,
                 const Rcpp::NumericVector &gvari, const Rcpp::NumericVector &gmean,
                 const Rcpp::IntegerVector &indexes, const Rcpp::IntegerVector &brsize,
                 const Rcpp::NumericVector &impdist, const Rcpp::NumericMatrix &bwmatrix,
		 const Rcpp::IntegerVector &revarray, const Rcpp::IntegerVector &badmissable,
		 const Rcpp::IntegerVector &bindex, const Rcpp::IntegerVector &ivdrop,
                 const llinputs llinfo, const int rep, const int niter,
                 const int gridmethod, const int vfinterpmethod,
                 Rcpp::NumericVector &rvfout) {

  int capj = llinfo.capj;
  int nbrand = llinfo.nbrand;
  int cmodel = llinfo.cmodel;
  int ncutiv = llinfo.ncutiv;
  int ncut = llinfo.ncutiv;
  int sptype = llinfo.sptype;
  int inttype = llinfo.inttype;
  int first20 = llinfo.first20;
  int nsave = llinfo.nsave;
  int ntilde = llinfo.ntilde;
  int ncrate = llinfo.ncrate; //need to add to the struct
  double dg = llinfo.dg;
  int hinvbound = llinfo.hinvbound;
  int hmodel = llinfo.hmodel;
  int initvf = llinfo.initvf;

  int lrgrid = rgrid.length();
  int livcoef = ivcoef.length();
  int lcutmat = cutmat.length();
  int lvf = vfsave.length();
  int livvari = ivvari.length();
  int ldetiv = detiv.length();

  int invmodel = llinfo.invmodel;
  int maxbottles = llinfo.maxbottles;

  int nd = ngrid.length(), nistates;
  nistates = ngrid[nd-1];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }

  int nsize = nd-1;

  int nrgrid = llinfo.nrgrid;

  int lomega = llinfo.lomega;

  int nb = llinfo.nb;
  int nbsmall = llinfo.nbsmall;
  int nco = llinfo.nco;
  int nhhs = llinfo.nhhs;
  //int *dummyb; // dummy integer pointer, but later this will tell use the composition of package holding
  int dummyb[1];
  *dummyb = 0;  // update this once we've built this in

  int ivblocksize;
  if(sptype == 1) {
    ivblocksize = (2*ncutiv+2)*nsize*nsize;
  } else {
    if(ncutiv==0) {
      ivblocksize = (nsize+1)*nsize;
    } else {
      ivblocksize = (ncutiv+4)*nsize*nsize;
    }
  }

  int cutblocksize = (ncutiv+1)*nsize;

  int nrep, nreptilde;

  int rep1 = rep+1; //conform to the R code

  if(rep1 <= nsave) {
    nrep = rep1-1;
  } else {
    nrep = nsave;
  }

  if(nrep < ntilde) {
    nreptilde = rep1-1;
  } else {
    nreptilde = ntilde;
  }

  if(initvf) {
    // we're supplying vf and parameter draw starting points
    nrep = nsave;
    nreptilde = ntilde;
    rep1 = nrep+1;
  }

  vfupdatepar vfpar(hhinds,bstates,pgrid,ngrid,co,vfsave,rgrid,cosave,
                    rgridsave,packsize,cutmat,ivcoef,ivvari,detiv,gvari,
                    gmean,indexes,brsize,impdist,bwmatrix,revarray,
		    badmissable,bindex,ivdrop,
                    capj,nbrand,cmodel,ncutiv,ncut,
                    sptype,inttype,first20,nsave,ntilde,ncrate,dg,
                    hinvbound,lrgrid,livcoef,lcutmat,lvf,ldetiv,nd,
                    nistates,istates,nsize,nrgrid,nb,nbsmall,nco,nhhs,
                    dummyb,ivblocksize,cutblocksize,nrep,nreptilde,
                    rep1,livvari,niter,gridmethod,vfinterpmethod,hmodel,
		    invmodel,maxbottles,lomega,rvfout);

  //parallelFor(hhinds[0]-1,hhinds[1],vfpar);

  if(hmodel) {
    // for homogeneous model, we can just do household 1, and then copy over
    // the VF for everyone else
    //  I tried parallelizing over inventory states but it is slower.
    vfpar(hhinds[0]-1,hhinds[0]);
    //parallelFor(0,nistates,vfpar);
    for(int hh=hhinds[0];hh<hhinds[1];hh++) {
      int vfind1 = (hh-(hhinds[0]-1))*nrgrid*nbsmall*nistates;
      for(int i=0;i<nrgrid*nistates*nbsmall;i++)
        rvfout[vfind1 + i] = rvfout[i];
    }
  } else {
    //vfpar(hhinds[0]-1,hhinds[1]);
    parallelFor(hhinds[0]-1,hhinds[1],vfpar);
  }

}



// main loop to do MCMC

// [[Rcpp::export]]
Rcpp::List MCMCLoops(Rcpp::List inputs, Rcpp::List Mcmc, Rcpp::List info,
                     Rcpp::List hhbig, Rcpp::List hhmerge, Rcpp::List hhmergebr,
		     Rcpp::IntegerVector hhinclude) {

  BEGIN_RCPP

  int maxrep = Rcpp::as< std::vector<int> >(Mcmc["R"])[0];
  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];
  int crfix = Rcpp::as< std::vector<int> >(info["crfix"])[0];
  int datacrate = Rcpp::as< std::vector<int> >(info["data.crate"])[0];
  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  int usecdraws = Rcpp::as< std::vector<int> >(info["usecdraws"])[0];
  int brandonly = Rcpp::as< std::vector<int> >(info["brandonly"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }
  Rcpp::NumericVector crate(nhhs);
  if(datacrate) {
    SEXP temp(info["crate"]);
    Rcpp::NumericVector temp1(temp);
    crate = clone(temp1);
  }

  SEXP xfull1 = info["xfull"];
  Rcpp::NumericVector xfull(xfull1);

  int nhhinclude = hhinclude.size();
  Rcpp::IntegerVector includeflag(nhhs);
  /*for(int i=0;i<nhhs;i++) {
    int j=0;
    while(j < nhhinclude && hhinclude(j) != i) {
      j++;
    }
    if(j == nhhinclude) {
      includeflag(i) = 0;
    } else {
      includeflag(i) = 1;
    }
    }*/

  for(int i=0;i<nhhinclude;i++) {
    includeflag(hhinclude(i)-1) = 1;
  }

  //Rcpp::Rcout << "maxrep: " << maxrep << std::endl;
  //Rcpp::Rcout << "nhhs: " << nhhs << std::endl;
  //Rcpp::Rcout << "xfull: ";

  int npbig = xfull.size();
  //for(int i=0;i<npbig;i++) {
  //  Rcpp::Rcout << xfull[i] << " ";
  //}
  //Rcpp::Rcout << std::endl;

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdropb(ivdrop1);
    ivdrop = clone(ivdropb);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = clone(ivdrop);
  

  SEXP tf1 = info["tform"];
  Rcpp::NumericVector tf(tf1);
  SEXP fixed1 = info["fixed"];
  Rcpp::IntegerVector fixed(fixed1);
  SEXP pequal1 = info["paramequal"];
  Rcpp::IntegerVector paramequal(pequal1);
  SEXP lbounds1 = info["lbounds"];
  Rcpp::NumericVector lbounds(lbounds1);
  SEXP ubounds1 = info["ubounds"];
  Rcpp::NumericVector ubounds(ubounds1);
  int crhhlb = Rcpp::as< std::vector<int> >(info["crhhlb"])[0];

  if(crhhlb) {
    if(info.containsElementNamed("crate")) {
      SEXP temp(info["crate"]);
      Rcpp::NumericVector temp1(temp);
      crate = clone(temp1);
    } else {
      throw std::range_error("crhhlb is TRUE but crate is undefined");
    }
  }

  if(!inputs.containsElementNamed("pstart"))
    throw std::range_error("pstart is NULL in inputs");

  SEXP inputpstart = inputs["pstart"];
  Rcpp::NumericMatrix pstart1(inputpstart);
  Rcpp::NumericMatrix pstart = Rcpp::clone(inputpstart);
  int npsmall = 0;
  for(int i = 0;i < npbig;i++) {
    npsmall += !fixed(i);
  }

  int useparamstart = 0;
  Rcpp::NumericMatrix paramstart(npbig,nhhs);
  if(info.containsElementNamed("paramstart")) {
    SEXP paramstart1 = info["paramstart"];
    Rcpp::NumericMatrix paramstart1a(paramstart1);
    pstart = Rcpp::clone(paramstart1a);
    paramstart = Rcpp::clone(paramstart1a);
    useparamstart = 1;
  }

  if(datacrate) {
    pstart(nbrand,Rcpp::_) = crate;
  }

  Rcpp::NumericMatrix pold(npsmall,nhhs);
  int j=-1;
  for(int i=0;i<npbig;i++) {
    if(!fixed(i)) {
      j++;
      pold(j,Rcpp::_) = pstart(i,Rcpp::_);
    }
  }

  Rcpp::NumericMatrix poldtf(npbig,nhhs);

  // this will not copy the object

  SEXP tunitstemp(hhmergebr["totunits"]);
  Rcpp::NumericVector tunitsbr(tunitstemp);

  SEXP panidtemp(hhmergebr["PANID"]);
  Rcpp::NumericVector panidbr(panidtemp);

  SEXP brindextemp(hhmergebr["brindex"]);
  Rcpp::IntegerVector brindexbr(brindextemp);

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  int pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  int retindll = Rcpp::as< std::vector<int> >(info["retindll"])[0];
  int varflag = 0;
  if(info.containsElementNamed("varflag"))
    varflag = Rcpp::as< std::vector<int> >(info["varflag"])[0];

  SEXP panidmtemp(hhmerge["PANID"]);
  Rcpp::NumericVector panidmerge(panidmtemp);

  Rcpp::IntegerVector hhinds(2);
  hhinds(0) = 1;
  hhinds(1) = nhhs;

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  SEXP pindsmerge1(info["idpriceindmerge"]);
  Rcpp::IntegerVector pindsmerge(pindsmerge1);

  SEXP pindsbig1(info["idpriceindbig"]);
  Rcpp::IntegerVector pindsbig(pindsbig1);

  int hmaxpind = Rcpp::as< std::vector<int> >(info["hmaxpind"])[0];
  SEXP hpindsbig1(info["hpriceindbig"]);
  Rcpp::IntegerVector hpindsbig(hpindsbig1);

  int maxpind = Rcpp::as< std::vector<int> >(info["maxpind"])[0];
  Rcpp::Rcout << "maxpind: " << maxpind << std::endl;
  int nobsbr = tunitsbr.size();
  int nobsmerge = panidmerge.size();

  int prindex = Rcpp::as< std::vector<int> >(info["prindex"])[0];
  Rcpp::NumericMatrix pricematbr(nobsbr,nbrand);
  SEXP ptemp = hhmergebr[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = hhmergebr[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }

  // this copies - should just work with numericmatrix directly
  //std::vector<double> pricevecbr = Rcpp::as< std::vector<double> >(pricematbr);

  Rcpp::NumericMatrix pricematmerge(nobsmerge,nbrand);
  SEXP ptemp2 = hhmerge[prindex-1];
  Rcpp::NumericVector ptemp3(ptemp2);
  for(int i=0;i<nbrand;i++) {
    ptemp2 = hhmerge[prindex-1+i];
    ptemp3 = ptemp2;
    pricematmerge(Rcpp::_,i) = ptemp3;
  }

  // define inputs for quantity log likelihood
  llinputs llinfo;

  SEXP hhbigcol1 = hhbig[0];
  Rcpp::NumericVector hhbigc1(hhbigcol1);
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];

  llinfo.nb = nb;
  llinfo.nbsmall = nbsmall;
  llinfo.nco = npbig;
  llinfo.nobs = hhbigc1.size();
  llinfo.nbrand = nbrand;
  llinfo.nsize = nsize;
  llinfo.nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  llinfo.sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  llinfo.ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  llinfo.inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];
  int nsave = 10;
  if(Mcmc.containsElementNamed("nsave")) {
    llinfo.nsave = Rcpp::as< std::vector<int> >(Mcmc["nsave"])[0];
    nsave = llinfo.nsave;
  } else {
    llinfo.nsave = 10;
  }
  if(Mcmc.containsElementNamed("ntilde")) {
    llinfo.ntilde = Rcpp::as< std::vector<int> >(Mcmc["ntilde"])[0];
  } else {
    llinfo.ntilde = 3;
  }

  int nitervf = 50; // number of iterations for value function update
  int niterstop = 500; // rep at which we switch to only a single update
  if(Mcmc.containsElementNamed("nitervf")) {
    nitervf = Rcpp::as< std::vector<int> >(Mcmc["nitervf"])[0];
  }

  if(Mcmc.containsElementNamed("niterstop")) {
    niterstop = Rcpp::as< std::vector<int> >(Mcmc["niterstop"])[0];
  }

  llinfo.nrgrid = Rcpp::as< std::vector<int> >(info["nrgrid"])[0];
  llinfo.necoef = Rcpp::as< std::vector<int> >(info["necoef"])[0];
  llinfo.capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];
  llinfo.nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];
  llinfo.ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];
  llinfo.retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];
  llinfo.myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];
  llinfo.cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];
  llinfo.debug = Rcpp::as< std::vector<int> >(info["debug"])[0];
  llinfo.usevf = 1;
  llinfo.idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];
  llinfo.first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];
  llinfo.hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];
  llinfo.genrnd = Rcpp::as< std::vector<int> >(info["genrnd"])[0];
  llinfo.nd = Rf_length(info["ngrid"]);
  llinfo.nq = INTEGER(Rf_getAttrib(info["qpts"], R_DimSymbol))[0];
  llinfo.ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];
  llinfo.pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  llinfo.hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  llinfo.repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0]-1;
  llinfo.ncpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[1];
  llinfo.rcpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[0];
  llinfo.ncrate = Rcpp::as< std::vector<int> >(info["ncost"])[0];
  llinfo.dg = Rcpp::as< std::vector<double> >(info["dg"])[0];
  llinfo.maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];
  llinfo.hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];
  llinfo.invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  llinfo.maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];
  llinfo.usecdraws = usecdraws;
  llinfo.initvf = 0;
  if(llinfo.invmodel == 1) {
    llinfo.lomega = 2;
  } else {
    llinfo.lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }

  int hmodel = llinfo.hmodel;

  Rcpp::NumericVector vfll(llinfo.nobs*nsize*(1+llinfo.capj));
  Rcpp::NumericVector utll(llinfo.nobs*nsize*(1+llinfo.capj));


  if(llinfo.debug && llinfo.retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  Rcpp::NumericVector cosave(npbig*nhhs*llinfo.nsave);
  Rcpp::IntegerVector indexes(llinfo.nsave);

  for(int i=0;i<llinfo.nsave;i++) {
    indexes(i) = i+1;
  }

  SEXP tunitsbigtemp(hhbig["totunits"]);
  Rcpp::NumericVector tunitsbig(tunitsbigtemp);

  SEXP panidbigtemp(hhbig["PANID"]);
  Rcpp::NumericVector panidbig(panidbigtemp);

  SEXP brindexbigtemp(hhbig["brindex"]);
  Rcpp::IntegerVector brindexbig(brindexbigtemp);

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP expandbig1(info["expandbig"]);
  Rcpp::IntegerVector expandbig(expandbig1);

  Rcpp::NumericMatrix ivbig(llinfo.nobs,nsize);
  Rcpp::NumericMatrix ivbig1(llinfo.nobs,nsize);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);

  //int nbstates = Rf_length(info["bstates"]);
  int maxgridlength = ngrid[ngrid.size()-1];

  int vfblocksize1 = llinfo.nsave*llinfo.nrgrid*nbsmall*maxgridlength;
  int vfblocksize2 = llinfo.nrgrid*nbsmall*maxgridlength;

  Rcpp::NumericVector vfsave(vfblocksize1*nhhs);

  if(info.containsElementNamed("vfstart")) {
    SEXP tempvfinfo(info["vfstart"]);
    Rcpp::NumericVector tempvf(tempvfinfo);
    vfsave = Rcpp::clone(tempvf);
    llinfo.initvf = 1;
  }

  // if initvf is true, copy in supplied cosave
  if(llinfo.initvf) {
    SEXP tempcosave1(info["cosave"]);
    Rcpp::NumericVector tempcosave(tempcosave1);
    cosave = Rcpp::clone(tempcosave);
  }

  //Rcpp::Rcout << "aaa" << std::endl;

  Rcpp::NumericVector vf(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  //SEXP cdraws1(info["cdraws"]);
  //Rcpp::NumericVector cdraws(cdraws1);

  // cdraws needs to be reddrawn so we redefine it here
  Rcpp::NumericVector cdraws(llinfo.nobs);

  // icdraws won't be used anymore except for debugging
  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);

  SEXP cdrawtemp(info["cdraws"]);
  Rcpp::NumericVector cdrawinfo(cdrawtemp);

  SEXP iiv1(info["initivdraws"]);
  Rcpp::NumericVector iiv(iiv1);

  SEXP ldraws1(info["logitdraws"]);
  Rcpp::NumericVector ldraws(ldraws1);

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);

  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);

  SEXP initb1(info["initb"]);
  Rcpp::IntegerVector initb(initb1);

  SEXP inits1(info["inits"]);
  Rcpp::IntegerVector inits(inits1);

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  SEXP qpts1(info["qpts"]);
  Rcpp::NumericVector qpts(qpts1);

  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);

  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);

  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);
  
  //std::vector<double> pricevecmerge = Rcpp::as< std::vector<double> >(pricematmerge);

  int llbrandsize = nhhs;
  if(varflag)
    llbrandsize = nobsbr;

  // declarations of variables that get returned/work variables

  Rcpp::NumericVector llbrand(llbrandsize);
  Rcpp::NumericVector llbrand1(llbrandsize);
  Rcpp::NumericVector llbrand1a(llbrandsize);

  Rcpp::NumericMatrix rivout(nobsmerge,nsize);
  Rcpp::NumericMatrix rivout1(nobsmerge,nsize);

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;
  Rcpp::NumericVector rcutmatout(cutblocksize*nhhs);
  Rcpp::NumericVector rcutmatout1(cutblocksize*nhhs);

  int nrowiv=0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivcoefblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  Rcpp::NumericVector rivcoefout(ivcoefblocksize*nhhs);
  Rcpp::NumericVector rivcoefout1(ivcoefblocksize*nhhs);

  int ivvariblocksize = nsize*nsize;
  Rcpp::NumericVector rivvariout(ivvariblocksize*nhhs);
  Rcpp::NumericVector rivvariout1(ivvariblocksize*nhhs);

  Rcpp::NumericVector rdetivout(nhhs);
  Rcpp::NumericVector rdetivout1(nhhs);

  // setup of grids, starting values for B and W,

  Rcpp::IntegerVector fixedpop(npsmall);
  if(Mcmc.containsElementNamed("fixedcoef")) {
    SEXP tempfixed(Mcmc["fixedcoef"]);
    Rcpp::IntegerVector tempfixed1(tempfixed);
    fixedpop = clone(tempfixed1);
  }
  //Rcpp::print(fixedpop);
  //Rcpp::Rcout << "abc" << std::endl;
  int nfixedpop = 0;
  for(int i=0;i<npsmall;i++) {
    nfixedpop += fixedpop(i);
  }
  //Rcpp::Rcout << "def" << std::endl;
  int allfixed = nfixedpop == npsmall;
  int allvarying = nfixedpop == 0;

  if(llinfo.hmodel && !allfixed) {
    Rcpp::Rcout << "Number of population fixed params: "  << nfixedpop << std::endl;
    throw std::runtime_error("if hmodel is true all parameters must be fixed");
  }


  int fixedmbound = nfixedpop;
  if(allvarying)
    fixedmbound = 1;

  // extract and define proposals

  arma::mat proposal(fixedmbound,fixedmbound);
  proposal.eye();
  if(Mcmc.containsElementNamed("proposal")) {
    SEXP tempprop1(Mcmc["proposal"]);
    Rcpp::NumericVector tempprop(tempprop1);
    if(tempprop.length() != nfixedpop*nfixedpop) {
      Rcpp::Rcout << "length of proposal: " << tempprop.length()
                 << "num fixed coefs: " << nfixedpop << std::endl;
      throw std::range_error("Incorrect dimensions for proposal.  Make sure the dimensions correspond to the number of population fixed coefficient.");
    }
    std::copy(tempprop.begin(),tempprop.end(),proposal.begin());
  }

  arma::mat cholprop = chol(proposal);

  // split up the proposal into 2 pieces

  int nbrfixed = 0;
  int ii1 = 0;
  for(int i=0;i<nbrand;i++) {
    if(!fixed(i)) {
      if(fixedpop(ii1)) {
        nbrfixed++;
      }
      ii1++;
    }
  }

  //Rcpp::Rcout << "nbrfixed: " << nbrfixed << " nfixedpop: " << nfixedpop << std::endl;

  int bend1 = nbrfixed;
  int bstart2 = nbrfixed;
  int bend2 = nfixedpop;

  if(nbrfixed == 0 || allvarying)
    bend1 = 1;

  arma::mat cpropbrand = cholprop.submat(0,0,bend1-1,bend1-1);

  if( nfixedpop == nbrfixed ) {
    bstart2 = nfixedpop-1;
    bend2 = nfixedpop;
  }

  if(allvarying) {
    bstart2 = 0;
    bend2 = 1;
  }

  arma::mat cpropother = cholprop.submat(bstart2,bstart2,bend2-1,bend2-1);

  arma::mat gvari1a(gvari.begin(),nsize,nsize,false);
  arma::mat gvar = arma::inv_sympd(gvari1a);
  arma::mat gchol = arma::chol(gvar);

  int nrgrid = llinfo.nrgrid;

  // gridmethod == 1: old random grid of inclusive values
  // gridmethod == 2: fixed grid of price vectors (randomly drawn in advance)
  int gridmethod = 1;
  if(Mcmc.containsElementNamed("gridmethod")) {
    gridmethod = Rcpp::as< std::vector<int> >(Mcmc["gridmethod"])[0];
  }

  int vfinterpmethod = 1;
  if(Mcmc.containsElementNamed("vfinterpmethod")) {
    vfinterpmethod = Rcpp::as< std::vector<int> >(Mcmc["vfinterpmethod"])[0];
  }

  Rcpp::NumericMatrix bwmatrix(npbig,npbig);

  double hthumb = pow(4/(3*( (double)nsave )),0.2);
  for(int i=0;i<npbig;i++) {
    bwmatrix(i,i) = 1/hthumb;
  }

  if(Mcmc.containsElementNamed("bwmatrix")) {
    SEXP bwtemp1(Mcmc["bwmatrix"]);
    Rcpp::NumericMatrix bwtemp(bwtemp1);
    bwmatrix = Rcpp::clone(bwtemp);
  }

  int rgridbdbig = gridmethod == 2 ? nbrand*llinfo.nrgrid*llinfo.nsave : nsize*llinfo.nrgrid*llinfo.nsave;

  int rgridbd = gridmethod == 2 ? nbrand : nsize;

  Rcpp::NumericMatrix rgrid(rgridbd,nrgrid);
  Rcpp::NumericVector rgridsave(rgridbdbig);


  if(gridmethod == 1) {
    for(int i=0;i<nsize;i++) {
      rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
    }
  } else if (gridmethod == 2) {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  Rcpp::NumericVector impdist(nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }

  arma::mat rgrida(rgrid.begin(),rgridbd,nrgrid,false);
  arma::colvec gmeana(gmean.begin(),gmean.size(),false);

  arma::mat tempm(nsize,nsize);
  if(gridmethod != 2) {
    tempm = rgrida.t()*gchol;
    rgrida = tempm.t();
    rgrida.each_col() += gmeana;
  }

  int nistates = ngrid[llinfo.nd-1];

  int nkeep = Rcpp::as< std::vector<int> >(Mcmc["keep"])[0];
  int nwrite = Rcpp::as< std::vector<int> >(Mcmc["nwrite"])[0];

  Rcpp::Environment base("package:base");
  Rcpp::Function save = base["save"];

  int ndrawkeep = maxrep/nkeep;

  // compute number of demographic vars

  SEXP dvarname1(info["dvarnames"]);
  Rcpp::StringVector dvarnames(dvarname1);

  int ndvars = dvarnames.size();

  Rcpp::NumericMatrix bdraw(ndrawkeep,npsmall*ndvars);

  Rcpp::NumericMatrix Wdraw(ndrawkeep,npsmall*npsmall);

  Rcpp::NumericVector llbrsave(ndrawkeep);
  Rcpp::NumericVector llqsave(ndrawkeep);

  // make these TRANSFORMED beta draws
  Rcpp::NumericVector betadraw(ndrawkeep*nhhs*npsmall);

  //Rcpp::List saveoutput;

  //saveoutput["bdraw"] = bdraw;
  //saveoutput["Wdraw"] = Wdraw;
  //saveoutput["llbrsave"] = llbrsave;
  //saveoutput["llqsave"] = llqsave;
  //saveoutput["betadraw"] = betadraw;

  Rcpp::NumericMatrix W(npsmall,npsmall);


  for(int i=0;i<npsmall;i++) {
    if(!fixedpop[i])
      W(i,i) = 1;
  }
  //Rcpp::Rcout << "ndvars: " << ndvars << std::endl;

  Rcpp::NumericMatrix b(npsmall,ndvars);

  for(int i=0;i<npsmall;i++) {
    b(i,0) = Rcpp::mean(pold(i,Rcpp::_));
  }

  // for speed, make a matrix of the demographic variables

  SEXP dvar1(info["dvars"]);
  Rcpp::NumericMatrix dvars(dvar1); //check that dims are right
  Rcpp::Rcout << "dvar rows: " << dvars.nrow() << ". cols: " << dvars.ncol() << std::endl;

  arma::mat W1(npsmall-nfixedpop,npsmall-nfixedpop);
  W1.zeros();

  arma::mat W1chol(npsmall-nfixedpop,npsmall-nfixedpop);
  W1chol.zeros();

  arma::mat W1inv(npsmall-nfixedpop,npsmall-nfixedpop);
  W1inv.zeros();

  arma::mat bnotfixed(npsmall-nfixedpop,nhhs);
  bnotfixed.zeros();

  arma::mat temppnew(npsmall-nfixedpop,nhhs);
  temppnew.zeros();

  arma::mat temppold(npsmall-nfixedpop,nhhinclude);
  temppold.zeros();

  // this is only needed for conformability
  arma::mat temppold1(npsmall-nfixedpop,nhhs);
  temppold1.zeros();

  Rcpp::NumericVector bdiff0(nhhs);
  Rcpp::NumericVector bdiff1(nhhs);

  arma::mat bdiff0a(bdiff0.begin(),bdiff0.size(),false);
  arma::mat bdiff1a(bdiff1.begin(),bdiff1.size(),false);

  arma::mat betamean(npsmall-nfixedpop,1);
  betamean.zeros();
  Rcpp::Rcout << "ndvars: " << ndvars << std::endl;
  arma::mat b1(npsmall-nfixedpop,ndvars);
  b1.zeros();
  Rcpp::Rcout << b1.n_rows << ", " << b1.n_cols << std::endl;
  arma::mat b2(nfixedpop,1);
  b2.zeros();

  arma::mat b2a(nbrfixed,1);
  b2a.zeros();

  arma::mat b2b(nfixedpop-nbrfixed,1);
  b2b.zeros();
  Rcpp::Rcout << "b2b size " << nfixedpop-nbrfixed << std::endl;

  arma::mat btempdiff(npsmall-nfixedpop,1);
  btempdiff.zeros();

  arma::mat bmeandiff(npsmall-nfixedpop,1);
  bmeandiff.zeros();

  Rcpp::NumericVector pnewdrawvec((npsmall-nfixedpop)*nhhs);
  arma::mat pnewdrawmat(pnewdrawvec.begin(),nhhs,npsmall-nfixedpop,false);

  Rcpp::NumericVector llqhh(nhhs);
  Rcpp::NumericVector llqhh1(nhhs);
  Rcpp::NumericVector llqhh1a(nhhs);
  Rcpp::NumericVector rllrun(llinfo.nsim*llinfo.nobs);
  Rcpp::NumericVector risave1(llinfo.nsim*nhhs*llinfo.ninitt);
  Rcpp::NumericVector risave2(llinfo.nsim*llinfo.nobs);

  Rcpp::NumericMatrix pnew(npsmall,nhhs);
  Rcpp::NumericMatrix pnewtf(npbig,nhhs);

  Rcpp::NumericVector r(nhhs);
  Rcpp::NumericVector u(nhhs);

  arma::mat npsmallnf(npsmall-nfixedpop,npsmall-nfixedpop);
  npsmallnf.eye();
  npsmallnf *= npsmall-nfixedpop;

  arma::mat S(npsmall-nfixedpop,npsmall-nfixedpop);
  S.zeros();

  arma::mat Schol(npsmall-nfixedpop,npsmall-nfixedpop);
  Schol.zeros();

  arma::mat XX(npsmall-nfixedpop,npsmall-nfixedpop);
  XX.zeros();

  arma::mat tempb1(1,npsmall - nfixedpop);
  tempb1.zeros();

  arma::mat tempb2(1,nfixedpop);
  tempb2.zeros();

  arma::mat tempb2a(1,nbrfixed);
  tempb2a.zeros();

  arma::mat tempb2b(1,nfixedpop-nbrfixed);
  tempb2b.zeros();

  // varying coefficient rho (1st el is brands, and second is other params)
  Rcpp::NumericVector rho(2);
  if(Mcmc.containsElementNamed("rho")) {
    SEXP temprho(Mcmc["rho"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho in Mcmc must be 2");
    rho = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho(i) = 0.1;
  }

  // fixed coefficients rho1
  Rcpp::NumericVector rho1(2);
  if(Mcmc.containsElementNamed("rho1")) {
    SEXP temprho(Mcmc["rho1"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho1 in Mcmc must be 2");
    rho1 = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho1(i) = 0.1;
  }

  // saves of rho1 - so we can see when it stabilizes
  Rcpp::NumericMatrix rho1save(ndrawkeep,2);

  int splitfixed = 1;
  if(Mcmc.containsElementNamed("splitfixed")) {
    splitfixed = Rcpp::as< std::vector<int> >(Mcmc["splitfixed"])[0];
  }


  int npropsave = 100;
  if(Mcmc.containsElementNamed("npropsave")) {
    npropsave = Rcpp::as< std::vector<int> >(Mcmc["npropsave"])[0];
  }

  int nprint = 10;
  if(Mcmc.containsElementNamed("nprint")) {
    nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  }

  std::vector<std::string> pnamevec(npbig,"V");
  for(int i=0;i<npbig;i++) {
    std::stringstream ss;
    ss << i+1;
    pnamevec[i] += ss.str();
  }

  Rcpp::CharacterVector pnames = Rcpp::wrap(pnamevec);

  if(Mcmc.containsElementNamed("pnames")) {
    SEXP pnames1(Mcmc["pnames"]);
    Rcpp::CharacterVector pnames2(pnames1);
    pnames = clone(pnames2);
  }
  //if(Mcmc.containsElementNamed("pnames")) {
  //  nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  //}

  Rcpp::IntegerMatrix xaccept(npropsave,2);
  Rcpp::IntegerVector acceptflag(nhhs);

  int wmethod = 2;  //inverse wishart drawing method - all seem about the same
  // 1: Rossi code (with demographics)
  // 2: my code (derived from ken train's R code)
  // 3: inverse gamma

  // priors

  // variance matrix prior

  arma::mat S0 = arma::eye(npsmall-nfixedpop,npsmall-nfixedpop);
  int nu = npsmall-nfixedpop+3;
  if(Mcmc.containsElementNamed("S0")) {
    SEXP temp(Mcmc["S0"]);
    Rcpp::NumericMatrix tempS(temp);
    for(int i=0;i<npsmall-nfixedpop;i++) {
      for(int j=0;j<npsmall-nfixedpop;j++) {
        S0(i,j) = tempS(i,j);
      }
    }
  }
  if(Mcmc.containsElementNamed("nu")) {
    nu = Rcpp::as< std::vector<int> >(Mcmc["nu"])[0];
  }
  S0 = S0*( (double)nu );

  //arma::mat Stemp(npsmall-nfixedpop,npsmall-nfixedpop+nhhs);
  arma::mat Stemp(npsmall-nfixedpop,nu+nhhinclude);
  Stemp.zeros();
  
  // mean parameter priors

  arma::mat bA = arma::eye(ndvars,ndvars);

  bA = bA*1000;
  if(Mcmc.containsElementNamed("bA")) {
    SEXP temp(Mcmc["bA"]);
    Rcpp::NumericMatrix tempbA(temp);
    bA = Rcpp::as<arma::mat>(tempbA);  // will clone
  }

  // precomputing to save time
  arma::mat RA = arma::chol(bA);
  arma::mat dvarm = Rcpp::as<arma::mat>(dvars);
  arma::mat dvarminc(nhhinclude,ndvars);
  ii1 = 0;
  for(int i=0;i<nhhs;i++) {
    if(includeflag(i)) {
      dvarminc.row(ii1) = dvarm.row(i);
      ii1++;
    }
  }
  arma::mat WW = arma::join_cols(dvarminc,RA);

  arma::mat betabar(npsmall-nfixedpop,ndvars);
  betabar.zeros();

  if(Mcmc.containsElementNamed("betabar")) {
    SEXP temp(Mcmc["betabar"]);
    Rcpp::NumericMatrix tempbb(temp);
    betabar = Rcpp::as<arma::mat>(tempbb);
  }
  arma::mat betabarvec = betabar;
  betabarvec.reshape(ndvars*(npsmall-nfixedpop),1);

  //this is a fudge when we have no demographic vars
  //We need to do it right with demographics
  double betaVA = 10.0;
  arma::mat betaA(npsmall-nfixedpop,npsmall-nfixedpop);
  betaA.eye();
  betaA = betaA*betaVA;
  if(Mcmc.containsElementNamed("betaA")) {
    SEXP temp(Mcmc["betaA"]);
    Rcpp::NumericMatrix tempbetaA(temp);
    betaA = Rcpp::as<arma::mat>(tempbetaA);  // will clone
  }
  arma::mat betaS0I = arma::inv_sympd(betaA);

  // fixed coefficient priors

  arma::mat bbarfixed(nfixedpop,1);
  bbarfixed.zeros();

  if(Mcmc.containsElementNamed("bbarfixed")) {
    SEXP temp(Mcmc["bbarfixed"]);
    Rcpp::NumericMatrix tempbbf(temp);
    bbarfixed = Rcpp::as<arma::mat>(tempbbf);
  }

  arma::mat bbarAfixed(nfixedpop,nfixedpop);
  bbarAfixed.eye();
  bbarAfixed = bbarAfixed*10.0;
  
  if(Mcmc.containsElementNamed("bbarAfixed")) {
    SEXP temp(Mcmc["bbarAfixed"]);
    Rcpp::NumericMatrix tempbbAf(temp);
    bbarAfixed = Rcpp::as<arma::mat>(tempbbAf);
  }
  arma::mat bbarAfI = arma::inv_sympd(bbarAfixed);

  //nu += nhhs;
  nu += nhhinclude;

  int naccept = 0;

  Rcpp::Rcout << "starting loop" <<std::endl;


  int keep=0;

  for(int rep=0;rep<maxrep;rep++) {

    int prcheck = 0;

    //debugprint = rep == 1334;
    //if(rep==0) {
    if(usecdraws) {
      cdraws = Rcpp::clone(cdrawinfo);
    } else {
      cdraws = Rcpp::runif(llinfo.nobs);
    }
    //}

    Rcpp::Timer timer;

    // initial likelihood calculations
    if(prcheck)
      Rcpp::Rcout << "rep " << rep << std::endl;

    tform(pold,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
          nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,
          brsize.begin(),sizebrand.begin(),poldtf);

    if(prcheck)
      Rcpp::Rcout << "aaa" << std::endl;

    //Rcpp::Rcout << "avg cons rate(1): " << Rcpp::mean(poldtf(nbrand,Rcpp::_)) << std::endl;

    //Rcpp::Rcout << "avg cons rate(2): " << Rcpp::mean(poldtf(nbrand+1,Rcpp::_)) << std::endl;

    brchoicellhh(poldtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
                 nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand);

    //Rcpp::Rcout << "temp ll: " << llbrand << std::endl;
    //throw std::range_error("stopping");

    if(prcheck)
      Rcpp::Rcout << "bbb" << std::endl;

    if(!brandonly) {
    
      if(hmodel) {
        fitiv(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
              obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
              ncutiv, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout);
      } else {
        fitivhh(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
                obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                ncutiv, dropmat, ivdrop, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout,
		ivdrop2);
      }

      for(int i=0;i<llinfo.nobs;i++) {
        for(int j=0;j<nsize;j++) {
          ivbig(i,j) = rivout(expandbig[i]-1,j);
        }
      }

      if(prcheck)
        Rcpp::Rcout << "ccc" << std::endl;

      if(hmodel) {
        llijc(poldtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
              hhinds,brsize,obshhinds,ivbig,vfsave,cdraws,ic,iiv,ldraws,
              packsize,initx,initb,inits,ngrid,rcutmatout,rivcoefout,rivvariout,rdetivout,
              pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,hpindsbig,
	      bstates,revarray,badmissable,bindex,ivdrop,
              llinfo,rep,gridmethod,vfinterpmethod,hmaxpind,
              vfll,utll,rllrun,risave1,risave2,llqhh);
      } else {
        llijc(poldtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
              hhinds,brsize,obshhinds,ivbig,vfsave,cdraws,ic,iiv,ldraws,
              packsize,initx,initb,inits,ngrid,rcutmatout,rivcoefout,rivvariout,rdetivout,
              pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,pindsbig,
	      bstates,revarray,badmissable,bindex,ivdrop2,
              llinfo,rep,gridmethod,vfinterpmethod,maxpind,
              vfll,utll,rllrun,risave1,risave2,llqhh);
      }
      //Rcpp::Rcout << "llqhh[0]:" << llqhh[0] << std::endl;

      //Rcpp::Rcout << "test ll sum" << Rcpp::sum(llqhh) << std::endl;
    }
    if(prcheck)
      Rcpp::Rcout << "ddd" << std::endl;

    int i1=0;

    if(!allfixed) {

      // draw population varying coefficients

      for(int i=0;i<npsmall;i++) {
        if(!fixedpop[i]) {
          int j1=0;
          for(int j=0;j<npsmall;j++) {
            if(!fixedpop[j]) {
              W1(i1,j1) = W(i,j);
              j1++;
            }
          }
          i1++;
        }
      }

      std::copy(pold.begin(),pold.end(),pnew.begin());

      i1=0;
      for(int i=0;i<npsmall;i++) {
        if(!fixedpop[i]) {
          for(int j=0;j<nhhs;j++) {
            temppnew(i1,j) = pold(i,j);
	    temppold1(i1,j) = pold(i,j);
          }
	  for(int j=0;j<nhhinclude;j++) {
            temppold(i1,j) = pold(i,hhinclude[j]-1);
	  }
          i1++;
        }
      }

      W1chol = arma::chol(W1);

      Rcpp::NumericVector tempdraws = Rcpp::rnorm((npsmall-nfixedpop)*nhhs);
      
      
      std::copy(tempdraws.begin(),tempdraws.end(),pnewdrawvec.begin());  // could copy direct to pnewdrawmat

      /*Rcpp::Rcout << "pnewdraws: ";
      for(int i=0;i<npsmall-nfixedpop;i++) {
	Rcpp::Rcout << "[0," << i << "]: " << pnewdrawmat(0,i) << "; ";
      }
      Rcpp::Rcout << std::endl;*/

      temppnew += arma::trans(pnewdrawmat*W1chol)*sqrt(rho(0));
     
      i1=0;
      for(int i=0;i<npsmall;i++) {
        if(!fixedpop[i]) {
          for(int j=0;j<nhhs;j++) {
            bnotfixed(i1,j) = b(i,0);
            pnew(i,j) = temppnew(i1,j);
            // add in demographics - can we vectorize this part?
            for(int k=1;k<ndvars;k++) {
              bnotfixed(i1,j) += dvars(j,k)*b(i,k);
            }
          }
          i1++;
        }
      }

      /*Rcpp::Rcout << "rho:" << sqrt(rho(0)) << std::endl;

      Rcpp::Rcout << W1chol;

      Rcpp::Rcout << "pold: ";
      for(int i=0;i<npsmall;i++) {
	Rcpp::Rcout << pold(i,0) << " ";
      }
      Rcpp::Rcout << std::endl;
      
      Rcpp::Rcout << "pnew: ";
      for(int i=0;i<npsmall;i++) {
	Rcpp::Rcout << pnew(i,0) << " ";
      }
      Rcpp::Rcout << std::endl;*/

      
      if(prcheck)
        Rcpp::Rcout << "eee" << std::endl;

      tform(pnew,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
            nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,
            brsize.begin(),sizebrand.begin(),pnewtf);

      if(!brandonly) {
      
        if(prcheck)
          Rcpp::Rcout << "fff" << std::endl;

        if(hmodel) {
          fitiv(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                ncutiv, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1);
        } else {
          fitivhh(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                  obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                  ncutiv, dropmat, ivdrop, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1,
		  ivdrop2);
        }

        for(int i=0;i<llinfo.nobs;i++) {
          for(int j=0;j<nsize;j++) {
            ivbig1(i,j) = rivout1(expandbig[i]-1,j);
          }
        }
      }
      
      if(prcheck)
        Rcpp::Rcout << "ggg" << std::endl;

      brchoicellhh(pnewtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
                   nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand1);

      if(prcheck)
        Rcpp::Rcout << "hhh" << std::endl;
      if(!brandonly) {
        if(hmodel) {
          llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,hpindsbig,
		bstates,revarray,badmissable,bindex,ivdrop,
                llinfo,rep,gridmethod,vfinterpmethod,hmaxpind,
                vfll,utll,rllrun,risave1,risave2,llqhh1);
        } else {
          llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,pindsbig,
		bstates,revarray,badmissable,bindex,ivdrop2,
                llinfo,rep,gridmethod,vfinterpmethod,maxpind,
                vfll,utll,rllrun,risave1,risave2,llqhh1);
        }
      }
      
      if(prcheck)
        Rcpp::Rcout << "iii" << std::endl;

      std::copy(llbrand1.begin(),llbrand1.end(),llbrand1a.begin());
      std::copy(llqhh1.begin(),llqhh1.end(),llqhh1a.begin());

      if(arma::rank( W1 ) < npsmall-nfixedpop ) {
        W1.print();
        Rcpp::Rcout << "reciprocal condition number: " << rcond(W1) << std::endl;
        throw std::range_error("W1 matrix in b draw singular");
      }

      W1inv = arma::inv_sympd(W1);
      //Rcpp::Rcout << "aaa" << std::endl;
      for(int i=0;i<nhhs;i++) {
        btempdiff = temppnew.col(i);
        bmeandiff = bnotfixed.col(i);
        bdiff1[i] = arma::as_scalar((btempdiff.t() - bmeandiff.t())*W1inv*(btempdiff - bmeandiff));

      }
      //Rcpp::Rcout << "abc" << std::endl;
      //Rcpp::Rcout << W1inv << std::endl;
      for(int i=0;i<nhhs;i++) {
	btempdiff = temppold1.col(i);
        bmeandiff = bnotfixed.col(i);
        bdiff0[i] = arma::as_scalar((btempdiff.t() - bmeandiff.t())*W1inv*(btempdiff - bmeandiff));
	//if(i<=1) {
	//  Rcpp::Rcout << btempdiff << std::endl;
	//  Rcpp::Rcout << bmeandiff << std::endl;
	//}
      }
      //Rcpp::Rcout << "bbb" << std::endl;
      r = exp(llbrand1 + llqhh1 - llbrand - llqhh - 0.5*bdiff1 + 0.5*bdiff0 );
      //Rcpp::Rcout << Rcpp::sum(llbrand) << " "<< Rcpp::sum(llbrand1) << " "
      //	  << Rcpp::sum(bdiff0) << " " << Rcpp::sum(bdiff1) << std::endl;
      //Rcpp::Rcout << "sum r: " << Rcpp::sum(r) << std::endl;
      //throw std::range_error("stopping");
      if(rep == -1) {
	Rcpp::Rcout << "household 522: " << r[521] << std::endl;
	Rcpp::Rcout << "llqhh: " << llqhh1[521] << ", " << llqhh[521] << ". " << std::endl;
	Rcpp::Rcout << "llbrand: " << llbrand1[521] << ", " << llbrand[521] << ". " << std::endl;
	Rcpp::Rcout << "bdiff: " << bdiff0[521] << ", " << bdiff1[521] << ". " << std::endl;
	arma::mat temp1 = temppnew.col(521);
	temp1.print();
	arma::mat temp2 = temppold1.col(521);
	temp2.print();
	Rcpp::Rcout << "pnewtf: " << pnewtf(41,521) << std::endl;
        //throw std::range_error("stopping");
      }
      
      u = Rcpp::runif(nhhs);

      naccept = 0;

      int i2=-1;
      for(int i=0;i<nhhs;i++) {

        acceptflag[i] = u[i] < dmin(1,r[i]) && !std::isnan(r[i]);
        if(includeflag[i])
          i2++;

        // double check isnan here
        // also, check that infinity is properly dealt with, truncate if not working
        if(u[i] < dmin(1,r[i]) && !std::isnan(r[i])) {
          // accept the draw if included hh
	  if(includeflag[i]) {
	    naccept++;
            i1=0;
            for(int j=0;j<npsmall;j++) {
              if(!fixedpop[j]) {
                temppold(i1,i2) = pnew(j,i);
                i1++;
              }
            }
          }
          pold(Rcpp::_,i) = pnew(Rcpp::_,i);
          poldtf(Rcpp::_,i) = pnewtf(Rcpp::_,i);
          ivbig(i,Rcpp::_) = ivbig1(i,Rcpp::_);
          llqhh[i] = llqhh1[i];
          llbrand[i] = llbrand1[i];

          std::copy(rcutmatout1.begin()+cutblocksize*i,rcutmatout1.begin()+cutblocksize*(i+1),rcutmatout.begin()+cutblocksize*i);
          std::copy(rivcoefout1.begin()+ivcoefblocksize*i,rivcoefout1.begin()+ivcoefblocksize*(i+1),rivcoefout.begin()+ivcoefblocksize*i);
          std::copy(rivvariout1.begin()+ivvariblocksize*i,rivvariout1.begin()+ivvariblocksize*(i+1),rivvariout.begin()+ivvariblocksize*i);
          rdetivout[i] = rdetivout1[i];

        }

      }

      //Rcpp::Rcout << "Number of acceptance: " << naccept << std::endl;
      

      if( ((double)naccept)/((double)nhhinclude) > 0.3) {
        rho(0) *= 1.1;
      } else {
        rho(0) *= 0.9;
      }

      // draw hyperparameters b and W
      // here I'm following code from Rossi Allenby and McCulloch's rmultireg
      // note # of equations = # of nonfixed utility coefs
      // # of vars is number of demographic vars (including the constant)

      //temppold: numpars by nhhinclude
      //betabar: ndvars by numpars
      //RA: ndvars by ndvars
      //ZZ should be nhhinclude + ndvars by numpars
      //WW is nhhinclude + ndvars by ndvars
      //IR is ndvars by ndvars
      //btilde is ndvars by numpars
      //E is nhhinclude + ndvars by numpars

      //arma::mat ZZ = arma::join_cols(temppold.t(), RA*betabar.t());

      //arma::mat IR = arma::solve(arma::trimatu(arma::chol(WW.t()*WW)),
      //                           arma::eye(ndvars,ndvars));

      //arma::mat btilde = (IR*IR.t())*(WW.t()*ZZ);

      //Rcpp::Rcout << "ZZ rows and cols" << ZZ.n_rows << ", " << ZZ.n_cols << std::endl;
      //Rcpp::Rcout << "WW rows and cols" << WW.n_rows << ", " << WW.n_cols << std::endl;
      //Rcpp::Rcout << "btilde length" << btilde.n_rows << ", " << btilde.n_cols << std::endl;

      // old
      //arma::mat E = ZZ - WW*btilde;

      //int nre = ndvars == 1 ? nhhinclude : nhhinclude + ndvars;
      int nre = nhhinclude;
      
      arma::mat E(nre,npsmall-nfixedpop);
      if(ndvars == 1) {
	// here define to be beta_n - betabar
	for(int i=0;i<nhhinclude;i++) {
	  i1=0;
	  for(int j=0;j<npsmall;j++) {
	    if(!fixedpop(j)) {
	      E(i,i1) = temppold(i1,i) - b(j,0);
	      i1++;
	    }
	  }
	  //for(int j=0;j<npsmall-nfixedpop;j++) {
	  //  E(i,j) = temppold(j,i) - b(,0);
	  //}
	}
      } else {
	//E = ZZ - WW*btilde; //I think something's wrong here, check this
	i1=0;
	for(int i=0;i<npsmall;i++) {
	  if(!fixedpop[i]) {
	    for(int j=0;j<nhhs;j++) {
	      bnotfixed(i1,j) = b(i,0);
	      for(int k=1;k<ndvars;k++) {
		bnotfixed(i1,j) += dvars(j,k)*b(i,k);
	      }
	    }
	    i1++;
	  }
	}

	//E.fill(-1);

	int i=0;
	for(int ii=0;ii<nhhs;ii++) {
	  if(includeflag(ii)) {
	    i1=0;
	    for(int j=0;j<npsmall;j++) {
	      if(!fixedpop(j)) {
		E(i,i1) = temppold(i1,i) - bnotfixed(i1,ii);
		i1++;
	      }
	    }
	    i++;
	  }
	}

	//Rcpp::Rcout << includeflag << std::endl;

	//E.print();
	//throw std::range_error("stopping");

      }

      // I've moved the calculation where we draw b to after W, since it
      // is a little easier

      if(wmethod == 1) {

        // this code is from Bayesm R package Rwishart

	// this should be robust to hhinclude

        // priors

        arma::mat V1 = E.t()*E + S0;

        if(arma::rank( V1 ) < npsmall-nfixedpop ) {
          W1.print();
          Rcpp::Rcout << "reciprocal condition number: " << rcond(W1) << std::endl;
          throw std::range_error("posterior matrix in W draw singular");
        }

        arma::mat V = arma::inv_sympd(V1);

        int m = npsmall-nfixedpop;
        arma::mat T = arma::zeros(m,m);

        for(int i = 0; i < m; i++) {
          T(i,i) = sqrt(Rcpp::rchisq(1,nu-i)[0]); //rchisq returns a vectorized object, so using [0] allows for the conversion to double
        }

        for(int j = 0; j < m; j++) {
          for(int i = j+1; i < m; i++) {
            T(i,j) = Rcpp::rnorm(1)[0]; //rnorm returns a NumericVector, so using [0] allows for conversion to double
          }
        }

        arma::mat C = arma::trans(T)*arma::chol(V);
        arma::mat CI = arma::solve(arma::trimatu(C),arma::eye(m,m));

        W1 = CI*arma::trans(CI);

      } else if (wmethod == 2) {

        // my old code based on Ken Train's R code

        //S = E.t()*E + npsmallnf;

	S = E.t()*E + S0;
	//Rcpp::Rcout << "E'E:" << E.t()*E << std::endl;
	//Rcpp::Rcout << "v0S0:" << S0 << std::endl;

	//Rcpp::Rcout << "S:" << S << std::endl;
        //Rcpp::NumericVector Stemp1 = Rcpp::rnorm( (npsmall-nfixedpop)*(npsmall-nfixedpop+nhhs) );
	//Rcpp::NumericVector Stemp1 = Rcpp::rnorm( (npsmall-nfixedpop)*(nu+nhhinclude) );
        //std::copy(Stemp1.begin(),Stemp1.end(),Stemp.begin());
	for(int j=0;j<nu;j++) {
	  for(int i=0;i<npsmall-nfixedpop;i++) {
	    Stemp(i,j) = Rcpp::rnorm(1)[0];
	  }
	}

	//Rcpp::Rcout << "Stemp(_,0): " << Stemp.col(0) << std::endl;

        if(arma::rank( S ) < npsmall-nfixedpop ) {
          S.print();
          Rcpp::Rcout << "reciprocal condition number: " << rcond(S) << std::endl;
	  arma::mat temp1 = E.t()*E;
	  temp1.print();
	  S0.print();
	  temppold.print();
	  bnotfixed.print();
	  Rcpp::Rcout << b << std::endl;
	  //E.print();
          throw std::range_error("S matrix in W draw singular");
        }

        Schol = arma::chol(arma::inv_sympd(S));

        XX = Schol.t()*Stemp;
	//arma::mat XX1 = Stemp.t()*Schol;
	//Schol.print();
        arma::mat X1 = XX*XX.t();
	//arma::mat X1 = XX1.t()*XX1;
	//X1.print();
	
        if(arma::rank( X1 ) < npsmall-nfixedpop ) {
          X1.print();
          Rcpp::Rcout << "reciprocal condition number: " << rcond(X1) << std::endl;
          throw std::range_error("X1 matrix in W draw singular");
        }
	// taking this inverse messes things up.  Why??
        W1 = arma::inv_sympd(X1);
	//W1.print();
	//Rcpp::Rcout << "aaa" << std::endl;
	//if(rep >= 100) {throw std::range_error("stopping");}
	//W1 = X1;
      } else if (wmethod == 3) {

        // Inverted gamma parameter by parameter - sometimes the covariances
        // seem hard to identify and this can fix the issue
        // here I follow my fortran code

        for(int i=0;i< npsmall-nfixedpop;i++) {

          Rcpp::NumericVector uu = Rcpp::rnorm(nu);

          double s = S0(i,i);
          for(int j=0;j<nhhinclude;j++)
            s += E(j,i)*E(j,i);
	  
          double r=0;
          for(int j=0;j<nu;j++)
            r+=uu(j)*uu(j);

          W1(i,i) = s/r;

        }

        //W1.print();

      }



      i1=0;
      for(int i=0;i<npsmall;i++) {
        if(!fixedpop[i]) {
          int j1=0;
          for(int j=0;j<npsmall;j++) {
            if(!fixedpop[j]) {
              W(i,j) = W1(i1,j1);
              j1++;
            }
          }
          i1++;
        }
      }

      arma::mat btildedraw = arma::mat(Rcpp::rnorm(ndvars*(npsmall-nfixedpop)));
      btildedraw.reshape(1,ndvars*(npsmall-nfixedpop));
      //Rcpp::Rcout << "btildedraw: " << btildedraw << std::endl;

      W1chol = arma::chol(W1);
      W1inv = arma::inv_sympd(W1);
      //Rcpp::Rcout << b1.n_rows << ", " << b1.n_cols << std::endl;
      if(ndvars == 1) {
	//Rcpp::Rcout << "abc" << std::endl;
	arma::mat tempbmean(npsmall-nfixedpop,1);
	tempbmean.zeros();
	for(int i=0;i<npsmall-nfixedpop;i++) {
	  for(int j=0;j<nhhinclude;j++) {
	    tempbmean(i,0) += temppold(i,j)/((double)nhhinclude);
	  }
	}

	arma::mat betaS1 = arma::inv_sympd(betaS0I + ((double)nhhinclude)*W1inv);
	tempbmean = betaS1*( betaS0I*betabar + ((double)nhhinclude)*W1inv*tempbmean);

	arma::mat betaS1C = arma::chol(betaS1);
	b1 = tempbmean + arma::trans(btildedraw*betaS1C);
	// if no priors used
	//b1 = tempbmean + arma::trans(btildedraw*W1chol)/sqrt((double)nhhinclude);
	//Rcpp::Rcout << "tempbmean " << tempbmean << std::endl;
      } else {
	//b1 = btilde + IR*btildedraw*W1chol;
	//b1 = b1.t();

	arma::mat xtx = arma::kron(arma::trans(dvarminc)*dvarminc,W1inv);

	arma::mat xty = W1inv*(temppold*dvarminc);
	
	int bigrow = ndvars*(npsmall-nfixedpop);
	xty.reshape(xty.n_rows*xty.n_cols,1);
	arma::mat ucholinv = arma::solve(arma::trimatu(arma::chol(xtx+betaS0I)), arma::eye(bigrow,bigrow));
	//Rcpp::Rcout << "abc" << std::endl;
	arma::mat betaS1 = ucholinv*arma::trans(ucholinv);
	arma::mat btilde = betaS1*(xty+betaS0I*betabarvec);
	//Rcpp::Rcout << "abcd" << std::endl;
	arma::mat bb = btilde + arma::trans(arma::chol(betaS1))*arma::trans(btildedraw);
	bb.reshape(npsmall-nfixedpop,ndvars);
	b1=bb;
      }

      //Rcpp::Rcout << "b and W draws" << std::endl;
      //Rcpp::Rcout << b1 << std::endl;
      //Rcpp::Rcout << W1 << std::endl;
      //if(rep == 2) {
      //throw std::range_error("stopping");
      //}
      //Rcpp::Rcout << "zzz" << std::endl;
      //Rcpp::Rcout << b.n_rows << ", " << b.n_cols << std::endl;
      //Rcpp::Rcout << b1.n_rows << ", " << b1.n_cols << std::endl;
      //Rcpp::Rcout << "ndvars: " << ndvars << std::endl;
      i1=0;
      for(int i=0;i<npsmall;i++) {
        if(!fixedpop[i]) {
          for(int j=0;j<ndvars;j++) {
            //Rcpp::Rcout << i << ", " << i1 << ", " << j << std::endl;
            b(i,j) = b1(i1,j);
          }
          i1++;
        }
      }
      //Rcpp::Rcout << "b1:" << std::endl;
      //b1.print();
      //Rcpp::Rcout << "b:" << std::endl;
      //Rcpp::Rcout << b << std::endl;
      //Rcpp::Rcout << "www" << std::endl;
    }

    if(!allvarying) {

      // draw population fixed coefficients

      if(splitfixed) {
	// separate fixed coefs for brands and dynamic pars

	for(int iis=0;iis<2;iis++) {

	  if(!( (iis == 0 && nbrfixed == 0) || (iis == 1 && nbrfixed == nfixedpop) ) ) {

	    i1=0;
	    for(int i=0;i<npsmall;i++) {
	      if(fixedpop(i)) {
		b2(i1,0) = b(i,0);
		i1++;
	      }
	    }

            //Rcpp::Rcout << "initial values for b and b2:" << std::endl;
            //Rcpp::Rcout << b << std::endl;
            //Rcpp::Rcout << b2 << std::endl;

	    if(iis == 0) {

	      Rcpp::NumericVector tempb12 = Rcpp::rnorm( nbrfixed );

	      std::copy(tempb12.begin(),tempb12.end(),tempb2a.begin());

	      i1=0;
              int i2=0;
	      for(int i=0;i<npsmall;i++) {
                if(fixedpop(i)) {
                  if(i1 < nbrfixed) {
                    b2a(i2,0) = b2(i1,0);
                    i2++;
                  }
                  i1++;
                }
              }

	      b2a += sqrt(rho1(0))*arma::trans(tempb2a*cpropbrand);

	      i1=0;
              i2=0;
	      for(int i=0;i<npsmall;i++) {
                if(fixedpop(i)) {
                  if(i1 < nbrfixed) {
                    b2(i1,0) = b2a(i2,0);
                    i2++;
                  }
                  i1++;
                }
              }

              //Rcpp::Rcout << "part a values for b2 and b2a:" << std::endl;
              //Rcpp::Rcout << b2 << std::endl;
              //Rcpp::Rcout << b2a << std::endl;

	    } else {

	      Rcpp::NumericVector tempb12 = Rcpp::rnorm( nfixedpop - nbrfixed );

	      std::copy(tempb12.begin(),tempb12.end(),tempb2b.begin());
              //Rcpp::Rcout << "abc" << std::endl;
	      i1=0;
              int i2=0;
	      for(int i=0;i<npsmall;i++) {
                if(fixedpop(i)) {
                  if(i1 >= nbrfixed) {
                    b2b(i2,0) = b2(i1,0);
                    i2++;
                  }
                  i1++;
                }
              }

              //Rcpp::Rcout << "b2b iteration (BEFORE) " << rep << std::endl;
              //b2b.print();
              //Rcpp::Rcout << "b2" << std::endl;
              //b2.print();

	      b2b += sqrt(rho1(1))*arma::trans(tempb2b*cpropother);
              //Rcpp::Rcout << "bcd" << std::endl;
	      i1=0;
              i2=0;
	      for(int i=0;i<npsmall;i++) {
                //Rcpp::Rcout << "i: " << i << ", i2: " << i2 << ", i3: " << i3 << std::endl;
                if(fixedpop(i)) {
                  if(i1 >= nbrfixed) {
                    b2(i1,0) = b2b(i2,0);
                    i2++;
                  }
                  i1++;
                }
	      }

              //Rcpp::Rcout << "b2b iteration " << rep << std::endl;
              //b2b.print();
              //Rcpp::Rcout << "b2" << std::endl;
              //b2.print();
              //Rcpp::Rcout << "b" << std::endl;
              //Rcpp::Rcout << b << std::endl;
              //throw std::runtime_error("stopping");

	    }

	    std::copy(pold.begin(),pold.end(),pnew.begin());

	    i1=0;
	    for(int i=0;i<npsmall;i++) {
	      if(fixedpop[i]) {
		for(int j=0;j<nhhs;j++)
		  pnew(i,j) = b2(i1,0);
		i1++;
	      }
	    }

	    tform(pnew,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
		  nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
                  sizebrand.begin(),pnewtf);
            if(!brandonly) {
              if(hmodel) {
                fitiv(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                      obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                      ncutiv, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1);
              } else {
                fitivhh(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                        obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                        ncutiv, dropmat, ivdrop, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1,
			ivdrop2);
              }

              for(int i=0;i<llinfo.nobs;i++) {
                for(int j=0;j<nsize;j++) {
                  ivbig1(i,j) = rivout1(expandbig[i]-1,j);
                }
              }
            }
            
	    brchoicellhh(pnewtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
			 nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand1);

            if(!brandonly) {
              if(hmodel) {
                llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                      hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                      packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                      pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,hpindsbig,
		      bstates,revarray,badmissable,bindex,ivdrop,
                      llinfo,rep,gridmethod,vfinterpmethod,hmaxpind,
                      vfll,utll,rllrun,risave1,risave2,llqhh1);
              } else {
                llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                      hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                      packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                      pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,pindsbig,
		      bstates,revarray,badmissable,bindex,ivdrop2,
                      llinfo,rep,gridmethod,vfinterpmethod,maxpind,
                      vfll,utll,rllrun,risave1,risave2,llqhh1);
              }
            }
	    double r1 = 0;

	    for(int i=0;i<nhhinclude;i++) {

	      r1 += llqhh1(hhinclude[i]-1) + llbrand1(hhinclude[i]-1)
		- ( llqhh(hhinclude[i]-1) + llbrand(hhinclude[i]-1) );

	    }

	    r1 = exp(r1);

	    double u1 = R::runif(0,1);

	    if(u1 < dmin(1,r1) && !std::isnan(r1)) {

	      i1=0;
	      for(int i=0;i<npsmall;i++) {
		if(fixedpop(i)) {
		  b(i,0) = b2(i1,0);
		  i1++;
		}
	      }

	      std::copy(pnew.begin(),pnew.end(),pold.begin());
	      std::copy(pnewtf.begin(),pnewtf.end(),poldtf.begin());

	      xaccept(rep%npropsave,iis) = 1;

	      std::copy(ivbig1.begin(),ivbig1.end(),ivbig.begin());

	      std::copy(rcutmatout1.begin(),rcutmatout1.end(),rcutmatout.begin());
	      std::copy(rivcoefout1.begin(),rivcoefout1.begin(),rivcoefout.begin());
	      std::copy(rivvariout1.begin(),rivvariout1.begin(),rivvariout.begin());
	      std::copy(rdetivout1.begin(),rdetivout1.end(),rdetivout.begin());

	      std::copy(llqhh1.begin(),llqhh1.end(),llqhh.begin());
	      std::copy(llbrand1.begin(),llbrand1.end(),llbrand.begin());

	    } else {

	      xaccept(rep%npropsave,iis) = 0;

	    }

	    if(rep%npropsave == 0 && rep > 0) {

	      double facceptfixed = 0.0;

	      for(int k=0;k<npropsave;k++)
		facceptfixed += xaccept(k,iis)/dmin(rep+1,npropsave);

	      if(facceptfixed < 0.3) {
		rho1(iis) *= 0.9;
	      } else {
		rho1(iis) *= 1.1;
	      }

	    }

	  }

	}

      } else {
	// draw all fixed coefs in one block

	i1=0;
	for(int i=0;i<npsmall;i++) {
	  if(fixedpop[i]) {
	    b2[i1] = b(i,0);
	    i1++;
	  }
	}

	arma::mat b2old = b2;

	/*Rcpp::NumericVector tempb12 = Rcpp::rnorm( nfixedpop );

	std::copy(tempb12.begin(),tempb12.end(),tempb2.begin());
	Rcpp::Rcout << "Fixed errors(1): " << tempb12 << std::endl;*/
	for(int i=0;i<nfixedpop;i++) {
	  tempb2(0,i) = Rcpp::rnorm(1)[0];
	}
	//Rcpp::Rcout << "Fixed errors: " << tempb2 << std::endl;

	//Rcpp::Rcout << "fixed rho: " << rho1(0) << std::endl;
	
	b2 += sqrt(rho1(0))*arma::trans(tempb2*cholprop);

	std::copy(pold.begin(),pold.end(),pnew.begin());

	i1=0;
	for(int i=0;i<npsmall;i++) {
	  if(fixedpop[i]) {
	    for(int j=0;j<nhhs;j++)
	      pnew(i,j) = b2[i1];
	    i1++;
	  }
	}

	tform(pnew,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
	      nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
              sizebrand.begin(),pnewtf);
        if(!brandonly) {
          if(hmodel) {
            fitiv(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                  obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                  ncutiv, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1);
          } else {
            fitivhh(pnewtf,  panidmerge, pricematmerge, hhinds, brsize,
                    obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
                    ncutiv, dropmat, ivdrop, rivout1, rcutmatout1, rivcoefout1, rivvariout1, rdetivout1,
		    ivdrop2);
          }

          for(int i=0;i<llinfo.nobs;i++) {
            for(int j=0;j<nsize;j++) {
              ivbig1(i,j) = rivout1(expandbig[i]-1,j);
            }
          }
        }
        
	brchoicellhh(pnewtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
		     nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand1);
        if(!brandonly) {
          if(hmodel) {
            llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                  hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                  packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                  pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,hpindsbig,
		  bstates,revarray,badmissable,bindex,ivdrop,
                  llinfo,rep,gridmethod,vfinterpmethod,hmaxpind,
                  vfll,utll,rllrun,risave1,risave2,llqhh1);
          } else {
            llijc(pnewtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
                  hhinds,brsize,obshhinds,ivbig1,vfsave,cdraws,ic,iiv,ldraws,
                  packsize,initx,initb,inits,ngrid,rcutmatout1,rivcoefout1,rivvariout1,rdetivout1,
                  pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,pindsbig,
		  bstates,revarray,badmissable,bindex,ivdrop2,
                  llinfo,rep,gridmethod,vfinterpmethod,maxpind,
                  vfll,utll,rllrun,risave1,risave2,llqhh1);
          }
        }
	//Rcpp::Rcout << "llqhh1[0]:" << llqhh1[0] << std::endl;

	double r1 = 0;

	for(int i=0;i<nhhinclude;i++) {

	  r1 += llqhh1[hhinclude[i]-1] + llbrand1[hhinclude[i]-1] -
                ( llqhh[hhinclude[i]-1] + llbrand[hhinclude[i]-1] );

	}

	double bdiff2 = arma::as_scalar( (b2.t() - bbarfixed.t())*bbarAfI*(b2 - bbarfixed) );
	double bdiff2old = arma::as_scalar( (b2old.t() - bbarfixed.t())*bbarAfI*(b2old - bbarfixed) );

	r1 += -0.5*bdiff2 + 0.5*bdiff2old;
	
	/*if(((rep+1)%nprint == 0 && nprint > 1) || nprint == 1) {

	  double r2 = 0;
	  double r3 = 0;
	  double r4 = 0;
	  double r5 = 0;
	  for(int i=0;i<nhhinclude;i++) {

	    r2 += llqhh1[hhinclude[i]-1] + llbrand1[hhinclude[i]-1];
	    r3 += llqhh[hhinclude[i]-1] + llbrand[hhinclude[i]-1];
	    r4 += llqhh[hhinclude[i]-1];
	    r5 += llbrand[hhinclude[i]-1];

	  }
	  Rcpp::Rcout << "r1: " << r1 << std::endl;
	  Rcpp::Rcout << "new ll: " << r2 << " old ll: " << r3  << ", " << r4 << ", " << r5 << std::endl;
	  }*/

	r1 = exp(r1);

	double u1 = R::runif(0,1);

	/*Rcpp::Rcout << "Fixed draws:" << std::endl;
	for(int i=0;i<npsmall;i++) {
	  if(fixedpop[i]) {
	    Rcpp::Rcout << b(i,0) << " ";
	  }
	  }*/
	//Rcpp::Rcout << std::endl;
	//Rcpp::Rcout << b2 << std::endl;
	//Rcpp::Rcout << "r1: " << r1 << ", u: " << u1 << std::endl;
	//if(rep == 9) {throw std::range_error("stopping");}

	if(u1 < dmin(1,r1) && !std::isnan(r1)) {
	  //if(r1 > 0) {

	  i1=0;
	  for(int i=0;i<npsmall;i++) {
	    if(fixedpop[i]) {
	      b(i,0) = b2[i1];
	      i1++;
	    }
	  }

	  std::copy(pnew.begin(),pnew.end(),pold.begin());
	  std::copy(pnewtf.begin(),pnewtf.end(),poldtf.begin());

	  xaccept(rep%npropsave,0) = 1;

	  std::copy(ivbig1.begin(),ivbig1.end(),ivbig.begin());

	  std::copy(rcutmatout1.begin(),rcutmatout1.end(),rcutmatout.begin());
	  std::copy(rivcoefout1.begin(),rivcoefout1.begin(),rivcoefout.begin());
	  std::copy(rivvariout1.begin(),rivvariout1.begin(),rivvariout.begin());
	  std::copy(rdetivout1.begin(),rdetivout1.end(),rdetivout.begin());

	  std::copy(llqhh1.begin(),llqhh1.end(),llqhh.begin());
	  std::copy(llbrand1.begin(),llbrand1.end(),llbrand.begin());

	} else {

	  xaccept(rep%npropsave,0) = 0;

	}

	if(rep%npropsave == 0 && rep > 0) {

	  double facceptfixed = 0.0;

	  for(int k=0;k<npropsave;k++)
	    facceptfixed += xaccept(k,0)/dmin(rep+1,npropsave);

	  //double facceptfixed = std::accumulate(xaccept2.begin(),xaccept2.end(),0.0)/dmin(rep+1,npropsave);

	  if(facceptfixed < 0.3) {
	    rho1(0) *= 0.9;
	  } else {
	    rho1(0) *= 1.1;
	  }

	}

      }

    }



    // update the value function

    int indx1 = rep%nsave;

    if(!llinfo.myopic && !brandonly) {
      if(rep == niterstop)
        nitervf = 1;
      //Rcpp::Rcout << "Updating VF" << std::endl;
      vfupdateijc(hhinds,bstates,pgrid,ngrid,poldtf,vfsave,rgrid,cosave,rgridsave,
                  packsize,rcutmatout,rivcoefout,rivvariout,rdetivout,gvari,gmean,
                  indexes,brsize,impdist,bwmatrix,revarray,badmissable,bindex,ivdrop2,
		  llinfo,rep,nitervf,gridmethod,
                  vfinterpmethod,vf);


      // see if can speed this up w std::copy - will need to reindex
      for(int hh=0;hh<nhhs;hh++) {
        for(int i=0;i<vfblocksize2;i++) {
          indcheck(hh*vfblocksize1 + indx1 + nsave*i,0,vfsave.length(),"vfsave copying bug: hh: %d; rep: %d; i: %d; vfblocksize1: %d; nsave: %d",hh,rep,i,
                   vfblocksize1,nsave);
          indcheck(hh*vfblocksize2+i,0,vf.length(),"vf copying bug: hh: %d; rep: %d; i: %d",hh,rep,i);
          vfsave[hh*vfblocksize1 + indx1 + nsave*i] = vf[hh*vfblocksize2+i];
        }
      }

    }

    for(int s=0;s<rgridbd;s++) {
      for(int j=0;j<nrgrid;j++) {
        indcheck(s+rgridbd*(j+nrgrid*indx1),0,rgridsave.length(),"rgridsave copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        indcheck(s+rgridbd*j,0,rgridsave.length(),"rgrid copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        rgridsave[s+rgridbd*(j+nrgrid*indx1)] = rgrid[s+rgridbd*j];
      }
    }

    for(int hh=0;hh<nhhs;hh++) {
      int pind = hh*npbig*nsave;
      for(int j=0;j<npbig;j++) {
        indcheck(j+pind+indx1*npbig,0,cosave.length(),"cosave copying bug: hh: %d; rep: %d; j: %d",hh,rep,j);
        indcheck(j+npbig*hh,0,poldtf.length(),"poldtf copying bug: hh: %d; rep: %d; j: %d",hh,rep,j);
        cosave[j+pind+indx1*npbig] = poldtf[j+npbig*hh];
      }
    }

    if(gridmethod != 2) {
      for(int i=0;i<rgridbd;i++) {
	rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
      }

      tempm = rgrida.t()*gchol;
      rgrida = tempm.t();
      rgrida.each_col() += gmeana;

    }

    if(nkeep == 1 || rep%nkeep == 0) {

      for(int i=0;i<b.nrow()*b.ncol();i++)
        bdraw(keep,i) = b[i];
      for(int i=0;i<W.nrow()*W.ncol();i++) {
        Wdraw[keep+ndrawkeep*i] = W[i];
      }

      int j=0;
      for(int i=0;i<npbig;i++) {
        if(!fixed(i)) {
          for(int hh=0;hh<nhhs;hh++) {
            betadraw[keep+ndrawkeep*(j+npsmall*hh)] = poldtf(i,hh);
          }
          j++;
        }
      }

      rho1save(keep,Rcpp::_) = rho1;

      //llqsave(keep) = std::accumulate(llqhh.begin(),llqhh.end(),0.0);
      //llbrsave(keep) = std::accumulate(llbrand.begin(),llbrand.end(),0.0);

      double rr = 0.0;
      for(int i=0;i<nhhinclude;i++) {
        rr += llbrand(hhinclude(i)-1);
      }
      llbrsave(keep) = rr;

      rr = 0.0;
      for(int i=0;i<nhhinclude;i++) {
        rr += llqhh(hhinclude(i)-1);
      }
      llqsave(keep) = rr;
      
      // save results to file
      if(rep%nwrite == 0) {
	// do we need to do this?? are the els of the list just pointers?

	//Rcpp::List saveoutput;

	base["bdraw"] = bdraw;
	base["Wdraw"] = Wdraw;
	base["llbrsave"] = llbrsave;
	base["llqsave"] = llqsave;
	base["betadraw"] = betadraw;

	std::vector<std::string> lst;
	lst.push_back("bdraw");
	lst.push_back("Wdraw");
	lst.push_back("llbrsave");
	lst.push_back("llqsave");
	lst.push_back("betadraw");
	Rcpp::CharacterVector saveoutput = Rcpp::wrap(lst);

	save(Rcpp::Named("list",saveoutput), Rcpp::Named("envir", base) , Rcpp::Named("file","temp_mcmc.RData"));

      }

      keep++;

    }
    base["poldtf"] = poldtf;
    base["vf"] = vf;
    //Rcpp::Rcout << "abc" << std::endl;
    //Rcpp::Rcout << "avg cons rate(1 end): " << Rcpp::mean(poldtf(nbrand,Rcpp::_)) << std::endl;

    //Rcpp::Rcout << "avg cons rate(2 end): " << Rcpp::mean(poldtf(nbrand+1,Rcpp::_)) << std::endl;

    if(((rep+1)%nprint == 0 && nprint > 1) || nprint == 1) {
      Rcpp::Rcout << "Information for iteration " << rep+1 << ":" << std::endl;

      Rcpp::NumericMatrix ans(npbig,4+ndvars);

      //Rcpp::colnames(ans) = Rcpp::CharacterVector::create("Median","0.05 %","0.95 %","b","diag(W^0.5)");

      Rcpp::CharacterVector colnamesans(4+ndvars);
      colnamesans(0) = "Median";
      colnamesans(1) = "0.05 %";
      colnamesans(2) = "0.95 %";
      for(int i=0;i<ndvars;i++) {
        colnamesans(3+i) = dvarnames(i);
      }
      colnamesans(3+ndvars) = "diag(W^0.5)";

      Rcpp::colnames(ans) = colnamesans;
      Rcpp::rownames(ans) = pnames;

      Rcpp::NumericMatrix bbig(npbig,ndvars);
      Rcpp::NumericVector Wdiag(npbig);

      i1=0;
      for(int i=0;i<npbig;i++) {
        if(fixed[i]) {
          bbig(i,0) = xfull[i];
          Wdiag[i] = 0;
        } else {
          for(int j=0;j<ndvars;j++)
            bbig(i,j) = b(i1,j);
          Wdiag[i] = sqrt(W(i1,i1));
          i1++;
        }
      }
      if(crfix) {
        for(int j=0;j<ndvars;j++)
          bbig(nbrand+1,j) = bbig(nbrand,j);
        Wdiag[nbrand+1] = Wdiag[nbrand];
      }

      double quds[] = {0.5,0.05,0.95};

      for(int j=0;j<3;j++) {
        i1=0;
        for(int i=0;i<npbig;i++) {
	  if(datacrate && (i == nbrand || i == nbrand+1)) {
	    ans(i,j) = quantile(poldtf(i,Rcpp::_),quds[j]);
	  } else {
	    if(fixed[i] || fixedpop[i1]) {
	      if(crfix && i == nbrand + 1) {
		if(fixed[i-1]) {
		  ans(i,j) = poldtf(i-1,0);
		} else {
		  ans(i,j) = quantile(poldtf(i-1,Rcpp::_),quds[j]);
		}
	      } else {
		//ans(i,j) = poldtf(i,0);
		ans(i,j) = quantile(poldtf(i,Rcpp::_),quds[j]);
	      }
	      if(!fixed[i]) i1++;
	    } else {
	      ans(i,j) = quantile(poldtf(i,Rcpp::_),quds[j]);
	      i1++;
	    }
	  }
        }
      }
      //ans(Rcpp::_,3) = bbig;
      for(int i=0;i<npbig;i++) {
        for(int j=0;j<ndvars;j++) {
          ans(i,3+j) = bbig(i,j);
        }
      }
      ans(Rcpp::_,3+ndvars) = Wdiag;
      //Rcpp::Rcout << "def" << std::endl;
      Rcpp::print(ans);
      Rcpp::Rcout << std::endl;

      Rcpp::Rcout << "End of Step Log Likelihoods:" << std::endl;
      //Rcpp::Rcout << "Brand: " << std::accumulate(llbrand.begin(),llbrand.end(),0.0) << std::endl;
      //Rcpp::Rcout << "Quantity: " << std::accumulate(llqhh.begin(),llqhh.end(),0.0) << std::endl;

      // new code with hhinclude
      double rr = 0.0;
      for(int i=0;i<nhhinclude;i++) {
        rr += llbrand(hhinclude(i)-1);
      }
      Rcpp::Rcout << "Brand: " << rr << std::endl;

      double rr1 = 0.0;
      for(int i=0;i<nhhinclude;i++) {
        rr1 += llqhh(hhinclude(i)-1);
      }
      Rcpp::Rcout << "Quantity: " << rr1 << std::endl;
      Rcpp::Rcout << "Total: " << rr+rr1 << std::endl;

      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << "Percent of varying acceptances: " << ((double)naccept)/((double)nhhs) << std::endl;
      //Rcpp::Rcout << "Percent of fixed acceptances: " << std::accumulate(xaccept.begin(),xaccept.end(),0.0)/dmin(rep+1,npropsave) << std::endl;
      Rcpp::Rcout << "Percent of fixed acceptances: ";
      for(int iis=0;iis<2;iis++) {
        double ff = 0;
        for(int k=0;k<npropsave;k++)
          ff += xaccept(k,iis)/dmin(rep+1,npropsave);
        Rcpp::Rcout << " Step " << iis+1 << ": " << ff;
      }
      Rcpp::Rcout << std::endl;

      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << "Varying coefficient RW Std Err: " << sqrt(rho(0)) << std::endl;
      Rcpp::Rcout << "Fixed coefficient RW Std Err: ";
      for(int iis=0;iis<2;iis++) {
        Rcpp::Rcout << " Step " << iis+1 << ": " << sqrt(rho1(iis));
      }

      //Rcpp::Rcout << "Fixed coefficient RW Std Err: " << sqrt(rho1) << std::endl;

      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << "Average Value Function:" << std::accumulate(vf.begin(),vf.end(),0.0)/((double)vf.size()) << std::endl;


      timer.step("Iteration");

      Rcpp::NumericVector tt(timer);

      tt[0] /= 1e9;
      Rcpp::Rcout << "Elapsed Seconds: "  << tt << std::endl;
      Rcpp::Rcout << std::endl;
      //throw std::runtime_error("stopping");
    }

  }


  Rcpp::List ret;
  ret["poldtf"] = poldtf;
  ret["pnewtf"] = pnewtf;
  ret["llbrand"] = llbrand;
  ret["llbrand1"] = llbrand1;
  ret["riv"] = rivout;
  ret["rivcoef"] = rivcoefout;
  ret["rivvari"] = rivvariout;
  ret["rdetiv"] = rdetivout;
  ret["llqhh"] = llqhh;
  ret["llqhh1"] = llqhh1;
  ret["naccept"] = naccept;
  ret["r"] = r;
  ret["b"] = b;
  ret["W"] = W;
  ret["acceptflag"] = acceptflag;
  ret["llqhh1a"] = llqhh1a;
  ret["llbrand1a"] = llbrand1a;
  ret["xaccept"] = xaccept;
  ret["vf"] = vf;
  ret["bdraw"] = bdraw;
  ret["Wdraw"] = Wdraw;
  ret["betadraw"] = betadraw;
  ret["rho1"] = rho1save;
  ret["pstart"] = pstart;

  Rcpp::NumericVector vfout(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  for(int hh=0;hh<nhhs;hh++) {
    for(int j=0;j<nrgrid;j++) {
      for(int i=0;i<nistates;i++) {
        for(int b=0;b<nbsmall;b++) {
          vfout(hh + nhhs*(j+nrgrid*(i+nistates*b))) = vf(hh*vfblocksize2 + j+nrgrid*(i+nistates*b));
        }
      }
    }
  }

  Rcpp::IntegerVector dims(4);
  dims(0) = nhhs;
  dims(1) = nrgrid;
  dims(2) = nistates;
  dims(3) = nbsmall;

  vfout.attr("dim") = dims;

  ret["vfsmall"] = vfout;
  ret["llqsave"] = llqsave;
  ret["llbrsave"] = llbrsave;

  Rcpp::IntegerVector dim2(3);
  dim2[0] = llinfo.nobs;
  dim2[1] = 1+llinfo.capj;
  dim2[2] = nsize;

  vfll.attr("dim") = dim2;
  utll.attr("dim") = dim2;

  ret["vfll"] = vfll;
  ret["utll"] = utll;

  ret["risave2"] = risave2;

  return ret;

  END_RCPP

}


// calculate value function using IJC style updating

// [[Rcpp::export]]
Rcpp::List computevf(Rcpp::NumericMatrix param, Rcpp::NumericVector nreptot, Rcpp::List info,
                     Rcpp::List hhbig, Rcpp::List hhmerge, Rcpp::List hhmergebr, Rcpp::NumericVector vmethod, Rcpp::IntegerVector grmethod1, Rcpp::IntegerVector ivfmeth1,
		     Rcpp::IntegerVector nsave1, Rcpp::NumericMatrix bwmatrix) {

  BEGIN_RCPP

    //int maxrep = Rcpp::as< std::vector<int> >(Mcmc["R"])[0];
    int nrep = Rcpp::as< std::vector<int> >(nreptot)[0];
  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];
  int crfix = Rcpp::as< std::vector<int> >(info["crfix"])[0];
  int datacrate = Rcpp::as< std::vector<int> >(info["data.crate"])[0];
  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }

  double tol = Rcpp::as< std::vector<double> >(info["tol"])[0];

  int method = Rcpp::as< std::vector<int> >(vmethod)[0];

  Rcpp::NumericVector crate(nhhs);
  if(datacrate) {
    SEXP temp(info["crate"]);
    Rcpp::NumericVector temp1(temp);
    crate = clone(temp1);
  }

  SEXP xfull1 = info["xfull"];
  Rcpp::NumericVector xfull(xfull1);

  int npbig = xfull.size();

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdropb(ivdrop1);
    ivdrop = clone(ivdropb);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = clone(ivdrop);

  SEXP tf1 = info["tform"];
  Rcpp::NumericVector tf(tf1);
  SEXP fixed1 = info["fixed"];
  Rcpp::IntegerVector fixed(fixed1);
  SEXP pequal1 = info["paramequal"];
  Rcpp::IntegerVector paramequal(pequal1);
  SEXP lbounds1 = info["lbounds"];
  Rcpp::NumericVector lbounds(lbounds1);
  SEXP ubounds1 = info["ubounds"];
  Rcpp::NumericVector ubounds(ubounds1);
  int crhhlb = Rcpp::as< std::vector<int> >(info["crhhlb"])[0];

  //SEXP inputpstart = inputs["pstart"];
  //Rcpp::NumericMatrix pstart1(inputpstart);
  //Rcpp::NumericMatrix pstart = Rcpp::clone(inputpstart);
  int npsmall = 0;
  for(int i = 0;i < npbig;i++) {
    npsmall += !fixed[i];
  }

  int useparamstart = 0;
  Rcpp::NumericMatrix paramstart(npbig,nhhs);
  if(info.containsElementNamed("paramstart")) {
    SEXP paramstart1 = info["paramstart"];
    Rcpp::NumericMatrix paramstart1a(paramstart1);
    //pstart = Rcpp::clone(paramstart1a);
    paramstart = Rcpp::clone(paramstart1a);
    useparamstart = 1;
  }


  //if(datacrate) {
  //  pstart(nbrand,Rcpp::_) = crate;
  //}

  //Rcpp::NumericMatrix pold(npsmall,nhhs);
  //int j=-1;
  //for(int i=0;i<npbig;i++) {
  //  if(!fixed[i]) {
  //    j++;
  //    pold(j,Rcpp::_) = pstart(i,Rcpp::_);
  //  }
  //}

  Rcpp::NumericMatrix poldtf(npbig,nhhs);

  // this will not copy the object

  SEXP tunitstemp(hhmergebr["totunits"]);
  Rcpp::NumericVector tunitsbr(tunitstemp);

  SEXP panidtemp(hhmergebr["PANID"]);
  Rcpp::NumericVector panidbr(panidtemp);

  SEXP brindextemp(hhmergebr["brindex"]);
  Rcpp::IntegerVector brindexbr(brindextemp);

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  //int pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  //int retindll = Rcpp::as< std::vector<int> >(info["retindll"])[0];
  int varflag = 0;
  if(info.containsElementNamed("varflag"))
    varflag = Rcpp::as< std::vector<int> >(info["varflag"])[0];

  SEXP panidmtemp(hhmerge["PANID"]);
  Rcpp::NumericVector panidmerge(panidmtemp);

  Rcpp::IntegerVector hhinds(2);
  hhinds[0] = 1;
  hhinds[1] = nhhs;

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  int nobsbr = tunitsbr.size();
  int nobsmerge = panidmerge.size();

  int prindex = Rcpp::as< std::vector<int> >(info["prindex"])[0];
  Rcpp::NumericMatrix pricematbr(nobsbr,nbrand);
  SEXP ptemp = hhmergebr[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = hhmergebr[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }

  // this copies - should just work with numericmatrix directly
  //std::vector<double> pricevecbr = Rcpp::as< std::vector<double> >(pricematbr);

  Rcpp::NumericMatrix pricematmerge(nobsmerge,nbrand);
  SEXP ptemp2 = hhmerge[prindex-1];
  Rcpp::NumericVector ptemp3(ptemp2);
  for(int i=0;i<nbrand;i++) {
    ptemp2 = hhmerge[prindex-1+i];
    ptemp3 = ptemp2;
    pricematmerge(Rcpp::_,i) = ptemp3;
  }

  // define inputs for quantity log likelihood
  llinputs llinfo;

  SEXP hhbigcol1 = hhbig[0];
  Rcpp::NumericVector hhbigc1(hhbigcol1);
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];

  llinfo.nb = nb;
  llinfo.nbsmall = nbsmall;
  llinfo.nco = npbig;
  llinfo.nobs = hhbigc1.size();
  llinfo.nbrand = nbrand;
  llinfo.nsize = nsize;
  llinfo.nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  llinfo.sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  llinfo.ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  llinfo.inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];
  int nsave = nsave1(0);
  //if(Mcmc.containsElementNamed("nsave")) {
  //  llinfo.nsave = Rcpp::as< std::vector<int> >(Mcmc["nsave"])[0];
  //  nsave = llinfo.nsave;
  //} else {
  //  llinfo.nsave = 10;
  //}
  //if(Mcmc.containsElementNamed("ntilde")) {
  //  llinfo.ntilde = Rcpp::as< std::vector<int> >(Mcmc["ntilde"])[0];
  //} else {
  //  llinfo.ntilde = 3;
  //}
  llinfo.nsave = nsave;
  llinfo.ntilde = nsave;

  llinfo.nrgrid = Rcpp::as< std::vector<int> >(info["nrgrid"])[0];
  llinfo.necoef = Rcpp::as< std::vector<int> >(info["necoef"])[0];
  llinfo.capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];
  llinfo.nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];
  llinfo.ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];
  llinfo.retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];
  llinfo.myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];
  llinfo.cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];
  llinfo.debug = Rcpp::as< std::vector<int> >(info["debug"])[0];
  llinfo.usevf = 1;
  llinfo.idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];
  llinfo.first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];
  llinfo.hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];
  llinfo.genrnd = Rcpp::as< std::vector<int> >(info["genrnd"])[0];
  llinfo.nd = Rf_length(info["ngrid"]);
  llinfo.nq = INTEGER(Rf_getAttrib(info["qpts"], R_DimSymbol))[0];
  llinfo.ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];
  llinfo.pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  llinfo.hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  llinfo.repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0]-1;
  llinfo.ncpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[1];
  llinfo.rcpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[0];
  llinfo.ncrate = Rcpp::as< std::vector<int> >(info["ncost"])[0];
  llinfo.dg = Rcpp::as< std::vector<double> >(info["dg"])[0];
  llinfo.maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];
  llinfo.hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];
  llinfo.invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  llinfo.maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];

  if(llinfo.invmodel == 1) {
    llinfo.lomega = 2;
  } else {
    llinfo.lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }

  llinfo.initvf = 0;

  if(llinfo.debug && llinfo.retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  int gridmethod = grmethod1(0);

  // note that either method 1 or 2 should give the same result, but
  // it can be useful to compare the two for debugging purposes
  int vfinterpmethod = ivfmeth1(0);

  if(bwmatrix.ncol() != npbig || bwmatrix.nrow() != npbig)
    throw std::runtime_error("Incorrect dimensions for bwmatrix");


  int rgridbdbig = gridmethod == 2 ? nbrand*llinfo.nrgrid*llinfo.nsave : nsize*llinfo.nrgrid*llinfo.nsave;

  int rgridbd = gridmethod == 2 ? nbrand : nsize;


  Rcpp::NumericVector impdist(llinfo.nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }


  Rcpp::NumericVector cosave(npbig*nhhs*llinfo.nsave);
  Rcpp::NumericVector rgridsave(rgridbdbig);
  Rcpp::IntegerVector indexes(llinfo.nsave);

  for(int i=0;i<llinfo.nsave;i++) {
    indexes[i] = i+1;
  }

  SEXP tunitsbigtemp(hhbig["totunits"]);
  Rcpp::NumericVector tunitsbig(tunitsbigtemp);

  SEXP panidbigtemp(hhbig["PANID"]);
  Rcpp::NumericVector panidbig(panidbigtemp);

  SEXP brindexbigtemp(hhbig["brindex"]);
  Rcpp::IntegerVector brindexbig(brindexbigtemp);

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP expandbig1(info["expandbig"]);
  Rcpp::IntegerVector expandbig(expandbig1);

  Rcpp::NumericMatrix ivbig(llinfo.nobs,nsize);
  Rcpp::NumericMatrix ivbig1(llinfo.nobs,nsize);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);

  int maxgridlength = ngrid[ngrid.size()-1];

  int vfblocksize1 = llinfo.nsave*llinfo.nrgrid*nbsmall*maxgridlength;
  int vfblocksize2 = llinfo.nrgrid*nbsmall*maxgridlength;
  Rcpp::Rcout << "vfblocksize1 " << vfblocksize1 << std::endl;
  Rcpp::NumericVector vfsave(vfblocksize1*nhhs);
  

  if(info.containsElementNamed("vfstart")) {
    SEXP tempvfinfo(info["vfstart"]);
    Rcpp::NumericVector tempvf(tempvfinfo);
    vfsave = Rcpp::clone(tempvf);
  }

  Rcpp::NumericVector vf(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);
  Rcpp::NumericVector vfold(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);
  Rcpp::NumericVector vfout(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  SEXP cdraws1(info["cdraws"]);
  Rcpp::NumericVector cdraws(cdraws1);

  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);

  SEXP iiv1(info["initivdraws"]);
  Rcpp::NumericVector iiv(iiv1);

  SEXP ldraws1(info["logitdraws"]);
  Rcpp::NumericVector ldraws(ldraws1);

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);

  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  SEXP qpts1(info["qpts"]);
  Rcpp::NumericVector qpts(qpts1);

  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);

  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);

  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);

  //std::vector<double> pricevecmerge = Rcpp::as< std::vector<double> >(pricematmerge);

  int llbrandsize = nhhs;
  if(varflag)
    llbrandsize = nobsbr;

  // declarations of variables that get returned/work variables

  Rcpp::NumericVector llbrand(llbrandsize);
  Rcpp::NumericVector llbrand1(llbrandsize);
  Rcpp::NumericVector llbrand1a(llbrandsize);

  Rcpp::NumericMatrix rivout(nobsmerge,nsize);
  Rcpp::NumericMatrix rivout1(nobsmerge,nsize);

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;
  Rcpp::NumericVector rcutmatout(cutblocksize*nhhs);
  Rcpp::NumericVector rcutmatout1(cutblocksize*nhhs);

  int nrowiv=0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize + 1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivcoefblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  Rcpp::NumericVector rivcoefout(ivcoefblocksize*nhhs);
  Rcpp::NumericVector rivcoefout1(ivcoefblocksize*nhhs);

  int ivvariblocksize = nsize*nsize;
  Rcpp::NumericVector rivvariout(ivvariblocksize*nhhs);
  Rcpp::NumericVector rivvariout1(ivvariblocksize*nhhs);

  Rcpp::NumericVector rdetivout(nhhs);
  Rcpp::NumericVector rdetivout1(nhhs);

  // setup of grids, starting values for B and W,

  //Rcpp::IntegerVector fixedpop(npsmall);
  //if(Mcmc.containsElementNamed("fixedcoef")) {
  //  SEXP tempfixed(Mcmc["fixedcoef"]);
  //  Rcpp::IntegerVector tempfixed1(tempfixed);
  //  fixedpop = clone(tempfixed1);
  //}

  //int nfixedpop = 0;
  //for(int i=0;i<npsmall;i++) {
  //  nfixedpop += fixedpop[i];
  //}

  //arma::mat proposal(nfixedpop,nfixedpop);
  //proposal.eye();
  //if(Mcmc.containsElementNamed("proposal")) {
  //  SEXP tempprop1(Mcmc["fixedcoef"]);
  //  Rcpp::NumericVector tempprop(tempprop1);
  //  if(tempprop.length() != nfixedpop*nfixedpop)
  //    throw std::range_error("Incorrect dimensions for proposal.  Make sure the dimensions correspond to the number of population fixed coefficient.");
  //  std::copy(tempprop.begin(),tempprop.end(),proposal.begin());
  //}

  //arma::mat cholprop = chol(proposal);
  
  arma::mat gvari1a(gvari.begin(),nsize,nsize,false);
  arma::mat gvar = arma::inv_sympd(gvari1a);
  arma::mat gchol = arma::chol(gvar);
  int nrgrid = llinfo.nrgrid;
  Rcpp::NumericMatrix rgrid(nsize,nrgrid);
  if(gridmethod==1) {
    for(int i=0;i<nsize;i++) {
      rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
    }
  } else {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  arma::mat rgrida(rgrid.begin(),rgridbd,nrgrid,false);
  arma::colvec gmeana(gmean.begin(),gmean.size(),false);

  arma::mat tempm(nsize,nsize);
  if(gridmethod != 2) {
    tempm = rgrida.t()*gchol;
    rgrida = tempm.t();
    rgrida.each_col() += gmeana;
  }

  int nistates = ngrid[llinfo.nd-1];

  //int nkeep = Rcpp::as< std::vector<int> >(Mcmc["keep"])[0];

  //int ndrawkeep = maxrep/nkeep;

  //Rcpp::NumericMatrix bdraw(ndrawkeep,npsmall);

  //Rcpp::NumericMatrix Wdraw(ndrawkeep,npsmall*npsmall);

  //Rcpp::NumericVector betadraw(ndrawkeep*nhhs*npsmall);

  //Rcpp::NumericMatrix W(npsmall,npsmall);

  //for(int i=0;i<npsmall;i++) {
  //  if(!fixedpop[i])
  //    W(i,i) = 1;
  //}

  //Rcpp::NumericVector b(npsmall);

  //for(int i=0;i<npsmall;i++) {
  //  b[i] = Rcpp::mean(pold(i,Rcpp::_));
  //}

  //arma::mat W1(npsmall-nfixedpop,npsmall-nfixedpop);
  //W1.zeros();

  //arma::mat W1chol(npsmall-nfixedpop,npsmall-nfixedpop);
  //W1chol.zeros();

  //arma::mat W1inv(npsmall-nfixedpop,npsmall-nfixedpop);
  //W1inv.zeros();

  //arma::colvec bnotfixed(npsmall-nfixedpop);
  //bnotfixed.zeros();

  //arma::mat temppnew(npsmall-nfixedpop,nhhs);
  //temppnew.zeros();

  //arma::mat temppold(npsmall-nfixedpop,nhhs);
  //temppold.zeros();

  //Rcpp::NumericVector bdiff0(nhhs);
  //Rcpp::NumericVector bdiff1(nhhs);

  //arma::mat bdiff0a(bdiff0.begin(),bdiff0.size(),false);
  //arma::mat bdiff1a(bdiff1.begin(),bdiff1.size(),false);

  //arma::mat betamean(npsmall-nfixedpop,1);
  //betamean.zeros();

  //arma::mat b1(npsmall-nfixedpop,1);
  //b1.zeros();

  //arma::mat b2(nfixedpop,1);
  //b2.zeros();

  //arma::mat btempdiff(npsmall-nfixedpop,1);
  //btempdiff.zeros();

  //Rcpp::NumericVector pnewdrawvec((npsmall-nfixedpop)*nhhs);
  //arma::mat pnewdrawmat(pnewdrawvec.begin(),npsmall-nfixedpop,nhhs,false);

  Rcpp::NumericVector llqhh(nhhs);
  Rcpp::NumericVector llqhh1(nhhs);
  Rcpp::NumericVector llqhh1a(nhhs);
  Rcpp::NumericVector rllrun(llinfo.nsim*llinfo.nobs);
  Rcpp::NumericVector risave1(llinfo.nsim*nhhs*llinfo.ninitt);
  Rcpp::NumericVector risave2(llinfo.nsim*llinfo.nobs);

  Rcpp::NumericMatrix pnew(npsmall,nhhs);
  Rcpp::NumericMatrix pnewtf(npbig,nhhs);

  Rcpp::NumericVector r(nhhs);
  Rcpp::NumericVector u(nhhs);

  //arma::mat npsmallnf(npsmall-nfixedpop,npsmall-nfixedpop);
  //npsmallnf.eye();
  //npsmallnf *= npsmall-nfixedpop;

  //arma::mat S(npsmall-nfixedpop,npsmall-nfixedpop);
  //S.zeros();

  //arma::mat Schol(npsmall-nfixedpop,npsmall-nfixedpop);
  //Schol.zeros();

  //arma::mat XX(npsmall-nfixedpop,npsmall-nfixedpop);
  //XX.zeros();

  //arma::mat Stemp(npsmall-nfixedpop,npsmall-nfixedpop+nhhs);
  //Stemp.zeros();

  //arma::mat tempb1(1,npsmall - nfixedpop);
  //tempb1.zeros();

  //arma::mat tempb2(1,nfixedpop);
  //tempb1.zeros();


  //double rho = 0.1;
  //if(Mcmc.containsElementNamed("rho")) {
  //  rho = Rcpp::as< std::vector<double> >(Mcmc["rho"])[0];
  //}

  //double rho1 = 0.1;
  //if(Mcmc.containsElementNamed("rho1")) {
  //  rho1 = Rcpp::as< std::vector<double> >(Mcmc["rho1"])[0];
  //}

  //int npropsave = 100;
  //if(Mcmc.containsElementNamed("npropsave")) {
  //  npropsave = Rcpp::as< std::vector<int> >(Mcmc["npropsave"])[0];
  //}

  //int nprint = 10;
  //if(Mcmc.containsElementNamed("nprint")) {
  //  nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  //}

  //std::vector<std::string> pnamevec(npbig,"V");
  //for(int i=0;i<npbig;i++) {
  //  std::stringstream ss;
  //  ss << i+1;
  //  pnamevec[i] += ss.str();
  //}

  //Rcpp::CharacterVector pnames = Rcpp::wrap(pnamevec);

  //if(Mcmc.containsElementNamed("pnames")) {
  //  nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  //}

  //Rcpp::IntegerVector xaccept(npropsave);
  //Rcpp::IntegerVector acceptflag(nhhs);


  // vf setup stuff

  tform(param,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
	nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
        sizebrand.begin(),poldtf);

  for(int indx1=0;indx1<nsave;indx1++) {
    for(int hh=0;hh<nhhs;hh++) {
      int pind = hh*npbig*nsave;
      for(int j=0;j<npbig;j++) {
        indcheck(j+pind+indx1*npbig,0,cosave.length(),"cosave copying bug: hh: %d; rep: %d; j: %d",hh,indx1,j);
        indcheck(j+npbig*hh,0,poldtf.length(),"poldtf copying bug: hh: %d; rep: %d; j: %d",hh,indx1,j);
        cosave[j+pind+indx1*npbig] = poldtf[j+npbig*hh];
      }
    }
  }

  if(llinfo.hmodel) {
    fitiv(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
          obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
          ncutiv, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout);
  } else {
    fitivhh(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
            obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
            ncutiv, dropmat, ivdrop, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout,
	    ivdrop2);
  }


  //int naccept = 0;

  Rcpp::Rcout << "starting VF loop" <<std::endl;

  Rcpp::NumericMatrix convinfo(nrep,2);
  Rcpp::NumericMatrix hhconvmean(nrep,nhhs);
  Rcpp::NumericMatrix hhconvmax(nrep,nhhs);

  int convrep = -1;

  for(int rep=0;rep<nrep;rep++) {

    std::copy( vf.begin(), vf.end(), vfold.begin() ) ;

    Rcpp::Timer timer;

    int indx1 = rep%nsave;

    vfupdateijc(hhinds,bstates,pgrid,ngrid,poldtf,vfsave,rgrid,cosave,rgridsave,
                packsize,rcutmatout,rivcoefout,rivvariout,rdetivout,gvari,gmean,
                indexes,brsize,impdist,bwmatrix,revarray,badmissable,bindex,ivdrop2,
		llinfo,rep,1,gridmethod,vfinterpmethod,vf);


    // see if can speed this up w std::copy - will need to reindex
    for(int hh=0;hh<nhhs;hh++) {
      for(int i=0;i<vfblocksize2;i++) {
        indcheck(hh*vfblocksize1 + indx1 + nsave*i,0,vfsave.length(),"vfsave copying bug: hh: %d; rep: %d; i: %d; vfblocksize1: %d; nsave: %d",hh,rep,i,
                 vfblocksize1,nsave);
        indcheck(hh*vfblocksize2+i,0,vf.length(),"vf copying bug: hh: %d; rep: %d; i: %d",hh,rep,i);
        vfsave[hh*vfblocksize1 + indx1 + nsave*i] = vf[hh*vfblocksize2+i];
      }
      for(int j=0;j<nrgrid;j++) {
        for(int i=0;i<nistates;i++) {
          for(int b=0;b<nbsmall;b++) {
            vfout(hh + nhhs*(j+nrgrid*(i+nistates*b))) = vf(hh*vfblocksize2 + j+nrgrid*(i+nistates*b));
          }
        }
      }
    }

    for(int s=0;s<rgridbd;s++) {
      for(int j=0;j<nrgrid;j++) {
        indcheck(s+rgridbd*(j+nrgrid*indx1),0,rgridsave.length(),"rgridsave copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        indcheck(s+rgridbd*j,0,rgridsave.length(),"rgrid copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        rgridsave[s+rgridbd*(j+nrgrid*indx1)] = rgrid[s+rgridbd*j];
      }
    }

    if(method == 1 && gridmethod != 2) {

      for(int i=0;i<nsize;i++) {
        rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
      }

      tempm = rgrida.t()*gchol;
      rgrida = tempm.t();
      rgrida.each_col() += gmeana;

    }

    timer.step("Iteration");

    Rcpp::NumericVector tt(timer);

    tt[0] /= 1e9;
    //Rcpp::Rcout << "VF Iteration " << rep+1 << ". Average Value Function: " << Rcpp::mean(vf) << ". Elapsed Seconds: "  << tt << std::endl;

    convinfo(rep,0) = Rcpp::mean(vf);
    convinfo(rep,1) = Rcpp::max(Rcpp::abs(vf-vfold));
    for(int hh=0;hh<nhhs;hh++) {
      hhconvmean(rep,hh) = 0;
      for(int i=0;i<vfblocksize2;i++) {
	hhconvmean(rep,hh) += vf[hh*vfblocksize2+i];
      }
      hhconvmean(rep,hh) /= (double)vfblocksize2;
      hhconvmax(rep,hh) = vf[hh*vfblocksize2];
      for(int i=0;i<vfblocksize2;i++) {
	hhconvmax(rep,hh) = vf[hh*vfblocksize2+i] > hhconvmax(rep,hh) ? vf[hh*vfblocksize2+i] : hhconvmax(rep,hh);
      }
    }

    Rcpp::Rcout << "VF Iteration " << rep+1 << ". Average Value Function: " << Rcpp::mean(vf) << ". Max: "  <<  Rcpp::max(vf) << ". diff: " << convinfo(rep,1) << std::endl;

    //if(rep > 0 && (convinfo(rep,1)-convinfo(rep-1,1))/convinfo(rep-1,1) < tol) {
    if(rep > 0 && convinfo(rep,1) < tol) {
      convrep = rep;
      break;
    }
    
  }

  Rcpp::List ret;
  ret["VF"] = vf;
  Rcpp::List ivret(4);
  ivret[0] = rivcoefout;
  ivret[1] = rcutmatout;
  ivret[2] = rivvariout;
  ivret[3] = rdetivout;

  ret["ivinfo"] = ivret;

  ret["ivbig"] = rivout;

  ret["vfbig"] = vfsave;

  ret["oldparam"] = cosave;

  ret["rgridsave"] = rgridsave;

  Rcpp::IntegerVector dims(4);
  dims(0) = nhhs;
  dims(1) = nrgrid;
  dims(2) = nistates;
  dims(3) = nbsmall;

  vfout.attr("dim") = dims;

  ret["vfsmall"] = vfout;

  Rcpp::List conv;

  conv["allinfo"] = convinfo;
  conv["hhmean"] = hhconvmean;
  conv["hhmax"] = hhconvmax;
  conv["convrep"] = convrep;
  
  ret["convinfo"] = conv;

  return ret;

  END_RCPP
}



// simulate choices to make sure estimation is working

// [[Rcpp::export]]
Rcpp::List choicesimijc(Rcpp::NumericMatrix param, Rcpp::List data, Rcpp::List info,
			Rcpp::List ivvars,
                        Rcpp::NumericMatrix ivbig, Rcpp::NumericVector vfin,
			Rcpp::NumericVector oldparam,
                        Rcpp::NumericVector oldrgrid, int nrgrid, Rcpp::IntegerVector indexes,
                        Rcpp::IntegerVector hhinds, int prindex, int gridmethod,
                        Rcpp::IntegerVector vfimethod, Rcpp::NumericMatrix bwmatrix,
			Rcpp::List hhinitmerge ,Rcpp::IntegerVector stvisitinit) {

  BEGIN_RCPP

    //double * co = REAL(param);

  //int * hhinds = INTEGER(hhinds1);

  //int prindex = *INTEGER(prindex1);

  //double * cosave = REAL(oldparam);

  //double * rgridsave = REAL(oldrgrid);

  //int rep = *INTEGER(rep1);
  int rep = 1;

  //int * indexes = INTEGER(indexes1);

  int nobs = Rf_length(VECTOR_ELT(data,0));

  int initnobs = Rf_length(VECTOR_ELT(hhinitmerge,0));

  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];
  int invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  int maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }

  Rcpp::NumericMatrix pricematbr(nobs,nbrand);
  SEXP ptemp = data[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = data[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }


  Rcpp::NumericMatrix pricematinit(initnobs,nbrand);
  SEXP ptempinit = hhinitmerge[prindex-1];
  Rcpp::NumericVector ptempinit1(ptempinit);
  for(int i=0;i<nbrand;i++) {
    ptempinit = hhinitmerge[prindex-1+i];
    ptempinit1 = ptempinit;
    pricematinit(Rcpp::_,i) = ptempinit1;
  }

  SEXP sizechoiceinit1 = hhinitmerge["sizechoice"];
  Rcpp::IntegerVector sizechoiceinit(sizechoiceinit1);

  SEXP nbottlesinit1 = hhinitmerge["nbottles"];
  Rcpp::IntegerVector nbottlesinit(nbottlesinit1);

  //int pflag = *INTEGER(info["print"]);

  SEXP tunitstemp(data["totunits"]);
  Rcpp::NumericVector tunits(tunitstemp);
  
  //double * tunits = REAL(data["totunits"]);
  //Rcpp::Rcout << "aaa" << std::endl;
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];

  int inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];

  double * iv[nsize];

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);
  //int * obshhinds = INTEGER(info["obshhinds"]);

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  //double *vfsave = REAL(vfin);

  //int nsave = *INTEGER(nsave1);

  //int ntilde = *INTEGER(ntilde1);

  //int nrgrid = *INTEGER(nrgrid1); //for the simulation, we might want to make nrgrid bigger?

  //int gridmethod = gridmethod1(0);

  int hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];

  int maxpind = 0;

  if(hmodel) {
    maxpind = Rcpp::as< std::vector<int> >(info["hmaxpind"])[0];
  } else {
    maxpind = Rcpp::as< std::vector<int> >(info["maxpind"])[0];
  }
  //Rcpp::Rcout << "abc" << std::endl;
  SEXP pindsbig1(info["idpriceindbig"]);
  Rcpp::IntegerVector pindsbig(pindsbig1);

  SEXP hpindsbig1(info["hpriceindbig"]);
  Rcpp::IntegerVector hpindsbig(hpindsbig1);

  Rcpp::IntegerVector pinds(pindsbig.size());
  if(hmodel) {
    pinds = Rcpp::clone(hpindsbig);
  } else {
    pinds = Rcpp::clone(pindsbig);
  }

  int rgridbd = gridmethod == 2 ? nbrand : nsize;

  Rcpp::NumericMatrix rgrid(rgridbd,nrgrid);
  //Rcpp::Rcout << "abc" << std::endl;
  if (gridmethod == 2) {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  int vfinterpmethod = vfimethod(0);



  Rcpp::NumericVector impdist(nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }




  int nsave = 1;
  int ntilde = 1;

  int nrowivbig = 0;

  for(int i=0;i<nsize;i++) {
    Rcpp::NumericVector temp = ivbig(Rcpp::_,i);

    iv[i] = temp.begin();
    nrowivbig = temp.size();

  }

  //double * brindex = REAL(data["brindex"]);

  SEXP panidtemp(data["PANID"]);
  Rcpp::NumericVector panid(panidtemp);
  //double * panid = REAL(data["PANID"]);

  SEXP initsize1(hhinitmerge["sizechoice"]);  // this looks like it repeats code above - is it even used?
  Rcpp::IntegerVector initsize(initsize1);

  SEXP initbottles1(hhinitmerge["nbottles"]);
  Rcpp::IntegerVector initbottles(initbottles1);
  
  //double * cdraws = REAL(info["cdraws"]);
  SEXP cdrawtemp(info["cdraws"]);
  Rcpp::NumericVector cdraws(cdrawtemp);

  int usecdraws = Rcpp::as< std::vector<int> >(info["usecdraws"])[0];

  //int necoef = *INTEGER(info["necoef"]);

  int capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];

  int nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];

  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];

  int ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];

  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];

  int retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdropb(ivdrop1);
    ivdrop = clone(ivdropb);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = clone(ivdrop);


  //double * ic = REAL(info["initcdraws"]);
  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);
  
  //double * iiv = REAL(info["initivdraws"]);

  //double * ldraws = REAL(info["logitdraws"]);

  //double * packsize = REAL(info["packsize"]);
  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);
  
  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);

  
  double maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];

  int myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];

  //double * ivdep = REAL(ivdep1);

  int cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];

  int debug = Rcpp::as< std::vector<int> >(info["debug"])[0];

  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];

  int idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];

  //int * brsize = INTEGER(info["brsize"]);
  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  int first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];

  int hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];

  //double * initx;
  //int linitx=0;
  //if (idrawtype == 1) {
  //  initx = REAL(info["initx"]);
  //  linitx = Rf_length(info["initx"]);
  //}
  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);
  int linitx = initx.length();

  //int genrnd = *INTEGER(info["genrnd"]);

  if(debug && retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  double ll = 0;

  //double * pgrid, * ivcoef, * cutmat, * ivvari, *detiv;
  //double *gvari, *gmean, dg;
  //int * ngrid, nd, nq, ncut;

  //double * pgrid = REAL(info["pgrid"]);
  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);
  //double * qpts = REAL(info["qpts"]);

  int siminitchoices = Rcpp::as< std::vector<int> >(info["siminitchoices"])[0];

  //double * ivcoef = REAL(VECTOR_ELT(ivvars,0));
  //double * cutmat = REAL(VECTOR_ELT(ivvars,1));
  //double * ivvari = REAL(VECTOR_ELT(ivvars,2));
  //double * detiv = REAL(VECTOR_ELT(ivvars,3));
  SEXP ivcoef1(ivvars[0]);
  Rcpp::NumericVector ivcoef(ivcoef1);
  SEXP cutmat1(ivvars[1]);
  Rcpp::NumericVector cutmat(cutmat1);
  SEXP ivvari1(ivvars[2]);
  Rcpp::NumericVector ivvari(ivvari1);
  SEXP detiv1(ivvars[3]);
  Rcpp::NumericVector detiv(detiv1);


  //int * ngrid = INTEGER(info["ngrid"]);
  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);
  int nd = Rf_length(info["ngrid"]);
  //int nq = INTEGER(Rf_getAttrib(info["qpts"], R_DimSymbol))[0];
  int ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];

  
  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);
  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);
  //double * gvari = REAL(info["gvari"]);
  //double * gmean = REAL(info["gmean"]);
  double dg = Rcpp::as< std::vector<int> >(info["dg"])[0];

  int nistates = ngrid(nd-1);

  double rgridhh[(nd-1)*nrgrid];

  double ivtprobs[nrgrid*maxpind];
  int ivtpfilled[maxpind];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }

  //int nrowivdep = 0;
  //if(sptype == 1) {
  //  nrowivdep = 3*ncutiv+3;
  //} else {
  //  nrowivdep = 2*ncutiv+5;
  //}

  int ivblocksize;
  if(sptype == 1) {
    if(ncutiv == 0) {
      ivblocksize = (2*ncutiv+2)*nsize;
    } else {
      ivblocksize = (2*ncutiv+2)*nsize*nsize;
    }
  } else {
    if(ncutiv == 0) {
      ivblocksize = (nsize + 1)*nsize;
    } else {
      ivblocksize = (ncutiv+4)*nsize*nsize;
    }
  }

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;

  int nrep, nreptilde;

  if(rep <= nsave) {
    nrep = rep-1;
  } else {
    nrep = nsave;
  }

  if(nrep < ntilde) {
    nreptilde = rep-1;
  } else {
    nreptilde = ntilde;
  }

  double initinv[nhhs];
  int initb[nhhs];

  memset(initb,0,nhhs*sizeof(int));

  int nco = INTEGER(Rf_getAttrib(param, R_DimSymbol))[0];

  //Rcpp::NumericMatrix bwmatrix(nco,nco);
  //Rcpp::Rcout << "abc" << std::endl;
  Rcpp::NumericVector llhh(nhhs);
  Rcpp::NumericVector llhhbr(nhhs);

  double ut[(1+capj)*nsize];
  int drop[(1+capj)*nsize];
  double vfhat[(1+capj)*nsize];
  double cerror[(1+capj)*nsize];
  int bchoice[(1+capj)*nsize];
  double ichoice[(1+capj)*nsize];

  int b=0;
  double ivout[nsize];
  double ivgrid[(1+capj)*nsize];
  double vftemp[nbsmall*nistates*maxpind];
  int vffilled[nbsmall*nistates*maxpind];
  //for(int i=0;i<nistates;i++) {
  //  for(int j=0;j<nb;j++) {
  //    vffilled[i+nistates*j] = -1;
  //  }
  //}
  //Rcpp::Rcout << "nbsmall: " << nbsmall << "; nistates: " << nistates << "; maxpind: " << maxpind << std::endl;
  // valgrind didn't like these lines - don't know why
  //memset(vffilled,-1,nbsmall*nistates*maxpind*sizeof(int));
  //memset(ivtpfilled,-1,maxpind*sizeof(int));
  for(int i=0;i<nbsmall*nistates*maxpind;i++) {
    vffilled[i] = -1;
  }
  for(int i=0;i<maxpind;i++) {
    ivtpfilled[i] = -1;
  }
  double invprime = 0;
  int coorder[nsave];

  //int lenic = Rf_length(info["initcdraws"]);
  //int leniiv = Rf_length(info["initivdraws"]);
  int lencdraws = Rf_length(info["cdraws"]);
  int leniv = Rf_length(data["PANID"]);

  if(lencdraws == 0 && usecdraws) {
    throw std::range_error("cdraws is length zero but usecdraws is true");
  }

  Rcpp::NumericVector rllrun(nsim*nobs);
  Rcpp::NumericVector risave1(nsim*nhhs*ninitt);
  Rcpp::List invsave(2);

  int hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  int repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0];

  int bprime = 0;

  // compute initial inventory
  // we might want to approach this a bit differently - probably should have an AR step??


  double clb, cub, alpha, gamma, nu, beta, ccost;
  double xi[nbrand];

  int lomega = 2;

  if(invmodel==2) {
    lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }

  double omega[lomega];
  double invub = maxinv;

  double inv = 0;
  
  // simulate initial choices

  int iobs = -1;
  int iobsbig = 0;

  double kweights[nsave];

  int vfind=0;

  Rcpp::NumericVector ichoiceout(initx.length());
  Rcpp::IntegerVector ichoicej(initx.length());
  Rcpp::IntegerVector ichoices(initx.length());
  Rcpp::IntegerVector ibsave(initx.length());
  Rcpp::IntegerVector inbr(initx.length());

  for(int n=hhinds(0)-1;n<hhinds(1);n++) {
    inv = 0;
    b = 0;
    for(int j=0;j<lomega;j++) {
      omega[j] = param(nbrand+9+j,n);
    }
    clb = param(nbrand,n);
    cub = param(nbrand+1,n);

    alpha = param(nbrand+2,n);
    gamma = param(nbrand+3,n);
    nu = param(nbrand+4,n);
    beta = param(nbrand+5,n);
    ccost = param(nbrand+7,n);

    vfind = n*nrgrid*nbsmall*nistates*nsave;

    if(gridmethod == 2) {
      for(int p=0;p<nrgrid;p++) {
        indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);

        computeiv(param.begin()+nco*n,rgrid.begin()+nbrand*p,brsize.begin(),
                  nsize,nbrand,rgridhh+(nd-1)*p);

      }
    }

    for(int t=0;t<ninitt;t++) {

      double chosenq = 0;
      double rcnew;
      if(usecdraws) {
	rcnew = ic(n+nhhs*t);
      } else {
	rcnew = Rcpp::runif(1)[0];
      }
      if(rcnew < 0.5) {
	rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
      } else {
	rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
      }
      indcheck(t+ninitt*n,0,linitx,"initx: n: %d; t %d",n,t);
      chosenq = initx(t+ninitt*n);  // since this vector is ordered by id, then week, this should be the right indexing
      // chosenq will get overwritten if we're simulating initial choices
    
      invub = hinvbound ? omega[0] : maxinv;
      
      double ivinit[nsize];

      if(stvisitinit(iobsbig)) {
	iobs++;
      }

      if(stvisitinit(iobsbig) && siminitchoices) {

	// If we are resimulating choices we should create new inclusive values

	int nprint = -1; //282;
	int tprint = -1; //1;
	
	int nbravail = 0;
	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << "****** n: " << n << ", t: " << t << std::endl;
	}
	int drop1[nbrand];
	
	for(int size=0;size<nsize;size++) {
	  double umax = -10000000;
	  for(int j=0;j<nbrand;j++) {
	    double price = pricematinit(iobs,j);
	    if(price >= 999 || brsize[j+nbrand] != size+1) {
	      drop1[j] = 1;
	    } else {
	      nbravail++;
	      ut[j] = param(j,n) + alpha*price;
	      drop1[j] = 0;
	      if(ut[j] > umax)
		umax = ut[j];
	      if(n == nprint && t <= tprint) {
		Rcpp::Rcout << " size: " << size << " ut[" << j << "]: " << ut[j] << " price: " << price << ". ";
	      }
	    }
	  }
	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << std::endl;
	  }
	  double s = 0;
	  for(int j=0;j<nbrand;j++) {
	    if(!drop1[j]) {
	      s += exp(ut[j]-umax);
	      if(n == nprint && t <= tprint) {
		Rcpp::Rcout << "j: " << j << "; ut[j]: " << ut[j] << "; udiff" << ut[j] - umax << "; s: " << s << std::endl;
	      }
	    }
	  }
	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << "s: " << s << "; umax: " << umax << std::endl;
	  }
	  s=log(s)+umax;
	  ivinit[size] = s;
	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << "ivinit[" << size << "]: " << ivinit[size] << std::endl;
	  }

	  for(int j=0;j<1+capj;j++) {
            ivgrid[j+size*(1+capj)] = ivinit[size]*((double)j)/((double)capj);
          }
	}

	//Note that the code below doesn't make use of price states.  It's not worth the trouble
	//since this part only gets run once

	if(!myopic) {
	  if(vfinterpmethod == 1) {
	    getclosestind(param.begin() + nco*n, oldparam.begin() + nsave*nco*n, indexes.begin(), nco, nrep, nsave, rep, coorder);
	  } else {
	    if(!(hmodel && n > hhinds[0]-1)) {
	      getkweights(param.begin() + nco*n, oldparam.begin() + nsave*nco*n,
			  indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
			  rep, coorder, kweights,0);
	    }
	  }
	}

	if(!myopic && beta > 0) {
          ivpred(ivinit, ivcoef.begin() + ivblocksize*n, cutmat.begin() + cutblocksize*n, nd, ncut*(ivdrop(n)==0), sptype, inttype, capj, ivout);
          if(gridmethod == 2) {
	    for(int p=0;p<nrgrid;p++) {
	      ivtprobs[p] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,ivvari.begin()+(nd-1)*(nd-1)*n,
						       detiv(n),impdist(p));
	      if(n==nprint && t <= tprint && p < 10) {
		  
		  Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p] << ", " << impdist(p) <<
		    ", " << detiv(n) << std::endl;
		  for(int ii=0;ii<nd-1;ii++) {
		    Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
		  }
		  Rcpp::Rcout << std::endl;
		}
            }
          }

	  int printflag = 0;

	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << "vfin: " << vfin[0] << ", " << vfin[1] << ", " << vfin[2] << std::endl;
	  }

          int nextbflag = 0;
	  //Rcpp::Rcout << "n, t: " << n << ", " << t << std::endl;
	  for(int size=0;size<nsize;size++) {
            for(int j=0;j<capj+1;j++) {
              if( (j == 0 && size == 0) || j > 0 ) {
		int bprime1 = 0;
		if(invmodel == 1) {
		  invprime = dmax(inv + packsize(size)*((double)j) - rcnew,0);
		  invprime = invprime > invub ? invub : invprime;
		} else {
		  nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		  bprime1 = bprime;
		  bprime = bindex[bprime]-1;
		  indcheck(bprime,0,nbsmall,"bprime (init sim) error");
		}

		int iind = itoiind(invprime,istates,nistates,first20);
                //int vfindx1 = iind + nistates*(bprime+nb*(pinds[n]-1));
                

                // uncomment if not using price states
                int vfindx1 = iind + nistates*bprime;
		int ok = iind < nistates; //check to make sure the state above this one is admissable (I think the second part should always be true)
                indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug (choice sim - 1). iind-1 %d, bprime %d",iind-1,bprime);

		if(invmodel == 2) {
		  ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
		}
		if(ok) {
		  if(vffilled[vfindx1] < 0) {
		    vffilled[vfindx1] = 1;
		    if(gridmethod==2) {
		      vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
						    ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
						    gvari.begin(), detiv(n), dg, gmean.begin(),
						    nd-1, istates, nistates, first20, bprime, nbsmall, 0,0,vfinterpmethod);

		    } else {
		      vftemp[vfindx1] = vfavgijchh(ivout, iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
						   ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
						   gvari.begin(), detiv(n), dg, gmean.begin(),
						   nd-1, istates, nistates, first20, bprime, nbsmall, printflag);
		    }

		    if(n == nprint && t <= tprint) {
		      Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		    }


		    //vfout(n+obsend*iind) = vftemp[vfindx1];

		  }
		}

                //vfindx1 = iind-1 + nistates*(bprime+nb*(pinds[n]-1));
                vfindx1 = iind-1 + nistates*bprime;
                indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug (choice sim - 2). iind-1 %d, bprime %d",iind-1,bprime);
                if(vffilled[vfindx1] < 0) {
                  vffilled[vfindx1] = 1;
                  if(gridmethod == 2) {
                    vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
                                                  ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
                                                  gvari.begin(), detiv(n), dg, gmean.begin(),
                                                  nd-1, istates, nistates, first20, bprime, nbsmall, 0, 0,vfinterpmethod);
                  } else {
                    vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
                                                 ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
                                                 gvari.begin(), detiv(n), dg, gmean.begin(),
                                                 nd-1, istates, nistates, first20, bprime, nbsmall, 0);
                  }

		  if(n == nprint && t <= tprint) {
		    Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		  }

                  //vfout(n+obsend*(iind-1)) = vftemp[vfindx1];

                }

              }
            }
          }
        }

        int first = 1;

        int vfoffset = 0;
        double maxutil = -1000000;
        int nextbflag = 0;
        
        for(int size=0;size<nsize;size++) {
          for(int j=0;j<capj+1;j++) {
            int indx = j+(1+capj)*size;
            drop[indx] = 1;
            if( (j == 0 && size == 0) || j > 0 ) {
              drop[indx] = 0;
	      if(invmodel == 1) {
		invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		invprime = invprime > invub ? invub : invprime;
		ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize(size), rcnew, gamma, nu, ccost, omega, cmodel, 0);
	      } else {
		nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		bchoice[indx] = bprime;
		indcheck(bprime,0,nb,"bprime (init sim) error");
		ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				   nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
		bprime = bindex[bprime]-1;		
		ichoice[indx] = invprime;
	      }
              
              if(!myopic && beta > 0) {
                vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20,vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
              } else {
                vfhat[indx] = 0.0;
              }
	      if(n == nprint && t <= tprint) {
		Rcpp::Rcout << "ut[" << indx << "]: " << ut[indx] << ", vf: " << vfhat[indx] << ", ivgrid: " << ivgrid[indx] <<std::endl;
	      }
              maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);

            }
          }
        }

        double s1 = 0.0;

	int schoice=0;
	int jchoice=0;

        double choiceumax = -1000000;

        for(int size=0;size<nsize;size++) {
          for(int j=0;j<capj+1;j++) {
            int indx = j+(1+capj)*size;
            if(!drop[indx]) {
              ut[indx] -= maxutil;
              cerror[indx] = -log(-log(R::runif(0,1)));
              s1 += exp(ut[indx]+vfhat[indx]);
              if(ut[indx]+vfhat[indx]+cerror[indx] > choiceumax) {
                choiceumax = ut[indx]+vfhat[indx]+cerror[indx];
		if(invmodel==2) {
		  bprime = bchoice[indx];
		  invprime = ichoice[indx];
		}
                if(j > 0) {
                  schoice = size+1;
		  jchoice = j;
                }
              }
            }
          }
        }

	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << "DEBUG: b: " << b << ", inv : " << inv <<
	    ", bprime: " << bprime << ", invprime; " << invprime <<
	    ", jchoice: " << jchoice << ", schoice: " << schoice << std::endl;

	  int tempb; double tempi;
	  if(invmodel==2) {	  
	    nextinv(jchoice,schoice-1,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		    rcnew,revarray.begin(),nsize,badmissable.begin(),&tempb,&tempi,&nextbflag,0);
	    Rcpp::Rcout << "TEST: " << tempb << ", " << tempi << std::endl;
	  }
	  for(int size=0;size<nsize;size++) {
	    for(int j=0;j<capj+1;j++) {
	      int indx = j+(1+capj)*size;
	      if(!drop[indx]) {
		Rcpp::Rcout << "choice: " << indx << " (" << j << ", " << size << "): utility: " << ut[indx]
			    << " vf: " << vfhat[indx] << " error: " << cerror[indx]
			    << " bchoice: " << bchoice[indx] << " ichoice: " << ichoice[indx]
			    << "choice prob: " << exp(ut[indx]+vfhat[indx])/s1 
			    << "tutil: " << ut[indx]+vfhat[indx]+cerror[indx] << std::endl;		
	      }
	    }
	  }
	  //if(t==tprint) {throw std::range_error("stopping");}
	}

        double obspack = 0;

        if(nbravail > 0 && schoice > 0) {
          chosenq = (double)(packsize(schoice-1)*jchoice);
	  ichoicej(t+ninitt*n) = jchoice;
	  ichoices(t+ninitt*n) = schoice;
        } else {
	  chosenq = 0;
	  ichoicej(t+ninitt*n) = 0;
	  ichoices(t+ninitt*n) = 0;
	  if(invmodel==2) {
	    bprime = bchoice[0];
	    invprime = ichoice[0];
	  }
	}

	// reset vffilled
	memset(vffilled,-1,nb*nistates*maxpind*sizeof(int));

      } else if(invmodel == 2) {
        int nextbflag = 0;
	if(siminitchoices || !stvisitinit(iobsbig)) {
	  nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
	} else {
	  nextinv(nbottlesinit(iobs),sizechoiceinit(iobs)-1,b,inv,bstates.begin(),
		  nb,maxbottles,packsize.begin(),
		  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
	}
      }

      ichoiceout(t+ninitt*n) = chosenq;

      //if(retinv)
      risave1[n+nhhs*t] = inv;
      ibsave(n+nhhs*t) = b;
      if(invmodel == 1) {
	inv = dmax(inv + chosenq - rcnew,0);
	inv = inv > invub ? invub : inv;
      } else {
	inv = invprime;
	b = bprime;
      }
      //if(n == 127 && invmodel == 2) {
      //Rcpp::Rcout << "init for " << n << ", " << t << ": inv: " << inv << ", b " << b <<
      //  ", j: " << ichoicej(t+ninitt*n) << ", s: " << ichoices(t+ninitt*n) <<
      //  ", stvisit: " << stvisitinit(iobsbig) << std::endl;
      //}
      
      iobsbig++;

      
    }
    initinv[n] = inv;
    if(invmodel==2) {
      initb[n] = b;
    }
  }

  // compute likelihood at data

  //int lomega = 2;  // update this later
  //double omega[lomega];

  vfind = 0;

  int hh = hhinds[0]-2;
  double llrun = 0;
  int obsend;
  if(hhinds(1)==nhhs) {
    obsend = nobs;
  } else {
    obsend = obshhinds(hhinds(1));
  }

  //double invub = maxinv;

  int nobsdata = obsend - (obshhinds(hhinds(0)-1)-1);
  Rcpp::IntegerVector simchoice(nobsdata*3);

  Rcpp::NumericVector risave2(nobsdata);
  Rcpp::IntegerVector bsave2(nobsdata);

  // value function, utilities we save to look at
  Rcpp::NumericVector vfout(obsend*nsize*(1+capj));
  Rcpp::NumericVector utout(obsend*nsize*(1+capj));

  // importance probs - for checking
  Rcpp::NumericVector tprobout(obsend*nrgrid);

  int n2 = obshhindsmerge(hhinds(0)-1)-1;

  //Rcpp::Rcout << "starting new loop" << std::endl;
  for(int hh=hhinds(0)-1;hh<hhinds(1);hh++) {
    //Rcpp::Rcout << "hh: " << hh << std::endl;

    if(!hmodel) {
      memset(vffilled,-1,nbsmall*nistates*maxpind*sizeof(int));
      memset(ivtpfilled,-1,maxpind*sizeof(int));
    }

    //indcheck(hh,0,nhhs,"initinv. n: %d; hh %d",n,hh);
    inv = initinv[hh];
    b = initb[hh];
    llrun = 0;
    llhh(hh) = 0;
    vfind = hh*nrgrid*nbsmall*nistates*nsave;

    for(int j=0;j<nbrand;j++) {
      xi[j] = param(j,hh);
    }

    clb = param(nbrand,hh);
    cub = param(nbrand+1,hh);

    alpha = param(nbrand+2,hh);
    gamma = param(nbrand+3,hh);
    nu = param(nbrand+4,hh);
    beta = param(nbrand+5,hh);
    ccost = param(nbrand+7,hh);

    for(int j=0;j<lomega;j++) {
      omega[j] = param(nbrand+9+j,hh);
    }

    invub = hinvbound ? omega[0] : maxinv;

    if(!myopic) {
      if(vfinterpmethod == 1) {
        getclosestind(param.begin() + nco*hh, oldparam.begin() + nsave*nco*hh, indexes.begin(), nco, nrep, nsave, rep, coorder);
      } else {
        if(!(hmodel && hh > hhinds[0]-1)) {
          getkweights(param.begin() + nco*hh, oldparam.begin() + nsave*nco*hh,
                      indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
                      rep, coorder, kweights,0);
        }
      }
    }

    int obsend;
    if(hh == nhhs-1) {
      obsend = nobs;
    } else {
      obsend = obshhinds(hh+1)-1;
    }

    if(gridmethod == 2) {
      for(int p=0;p<nrgrid;p++) {
        indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);

        computeiv(param.begin()+nco*hh,rgrid.begin()+nbrand*p,brsize.begin(),
                  nsize,nbrand,rgridhh+(nd-1)*p);

      }
    }
    //Rcpp::Rcout << "abc" << std::endl;
    for(int n=obshhinds(hh)-1;n<obsend;n++) {

      int units = (int)tunits(n);

      double rcnew;
      //indcheck(n,0,lencdraws,"cdraws. n: %d",n);
      if(usecdraws) {
        rcnew = cdraws(n);
      } else {
        rcnew = Rcpp::runif(1)[0];
      }

      if(rcnew < 0.5) {
        rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
      } else {
        rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
      }

      int nextbflag = 0;

      if(units < 0) {
        // household didn't go to the store - there's nothing to model here.
        risave2[n] = inv;
	if(invmodel==1) {
	  inv = dmax(inv - rcnew,0);
	} else {
	  nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		  rcnew,revarray.begin(),nsize,badmissable.begin(),&b,&inv,
                  &nextbflag,0);
	}
        if(debug) {
          rllrun(n) = 1;
        }
        simchoice(n) = 0;
        simchoice(n+nobsdata) = -1;
        simchoice(n+nobsdata*2) = -1;
      } else {
        indcheck(n,0,leniv,"iv. n: %d",n);
        //Rcpp::Rcout << "1" << std::endl;
        double maxutil = -1000000;

        // compute avg predicted inclusive value, needed for importance weighting
        double rgrid1[nsize];
        //indcheck(n2,0,nrowivbig,"ivbig: n: %d",n);
        for(int size=0;size<nsize;size++) {
          rgrid1[size] = ivbig(n,size);
        }

	//int printobs = n <= 51141 ? 51141 : 52607;
	int printobs = -1; //38943; //n <= 158 ? 158 : 17378;;

        if(!myopic && beta > 0) {
          ivpred(rgrid1, ivcoef.begin() + ivblocksize*hh, cutmat.begin() + cutblocksize*hh, nd, ncut*(ivdrop(hh)==0), sptype, inttype, capj, ivout);
          if(gridmethod == 2) {
	    if(n==printobs) {
	      Rcpp::Rcout << "pinds[" << n << "]: " << pinds[n] << "; ivptfilled[" << n << "]: " << ivtpfilled[pinds[n]-1] << std::endl;
	      for(int ii=0;ii<nsize;ii++) {
		Rcpp::Rcout << "rgrid1[" << ii << "]: " << rgrid1[ii] << ", ";
	      }
	      Rcpp::Rcout << std::endl;
	      for(int ii=0;ii<nd-1;ii++) {
		Rcpp::Rcout << "ivout[" << ii << "]: " << ivout[ii] << ", ";
	      }
	      Rcpp::Rcout << std::endl;
	      Rcpp::Rcout << "ivcoef:" << std::endl;
	      for(int ii =0;ii<ivblocksize;ii++) {
		Rcpp::Rcout << "[" << ii << "]: " << *(ivcoef.begin() + ivblocksize*hh + ii) << "; ";
	      }
	      Rcpp::Rcout << std::endl;
	    }
            if(ivtpfilled[pinds[n]-1] < 0) {
              ivtpfilled[pinds[n]-1] = n;
              for(int p=0;p<nrgrid;p++) {
                ivtprobs[p+nrgrid*(pinds[n]-1)] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,ivvari.begin()+(nd-1)*(nd-1)*hh,
							 detiv(hh),impdist(p));
                tprobout(n+obsend*p) = ivtprobs[p];
		if(n==printobs && p < 10) {
		  
		  Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p+nrgrid*(pinds[n]-1)] << ", " << impdist(p) <<
		    ", " << detiv(hh) << std::endl;
		  for(int ii=0;ii<nd-1;ii++) {
		    Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
		  }
		  Rcpp::Rcout << std::endl;
		}
              }
            }
          }
        }
        //Rcpp::Rcout << "2" << std::endl;
        
        double iv0[nsize];
        for(int size=0;size<nsize;size++) {
          iv0[size] = ivbig(n,size);
          for(int j=0;j<1+capj;j++) {
            ivgrid[j+size*(1+capj)] = iv0[size]*((double)j)/((double)capj);
            if(n==printobs) {
              Rcpp::Rcout << "iv(" << n << "," << size << "): " << ivbig(n,size) <<
                " iv0[" << size << "]: " << iv0[size] << " ivgrid[" <<
                j+size*(1+capj) << "]: " << ivgrid[j+size*(1+capj)] << std::endl;
            }
          }
        }
        //Rcpp::Rcout << "3" << std::endl;
        // compute value function at points that get visited - could fold this into next loop

        int printflag = hh == hhprint && rep == repprint;

        //int pobs = 0; //n == 1 || n == obshhinds[hh+1]-2;

	if(n==printobs) {
	  Rcpp::Rcout << "inv: " << inv << ", b: " << b << ", init inv: " << initinv[hh] << ", init b: " << initb[hh] << std::endl;
	}
	
        if(!myopic && beta > 0) {
          for(int size=0;size<nsize;size++) {
            for(int j=0;j<capj+1;j++) {
              if( (j == 0 && size == 0) || j > 0 ) {
		int bprime1 = 0;
		if(invmodel == 1) {
		  invprime = dmax(inv + packsize(size)*((double)j) - rcnew,0);
		  //drop[j+(1+capj)*size] = invprime > maxinv;
		  invprime = invprime > invub ? invub : invprime;
		} else {
		  nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                          &nextbflag,0);
		  bprime1 = bprime;
		  bprime = bindex[bprime]-1;
		  indcheck(bprime,0,nb,"bprime (sim) error");
		}

		int iind = itoiind(invprime,istates,nistates,first20);
                /*int iind = 0;
                double dif = 0;
                if(invprime > istates[first20-1]) {
                  dif = istates[first20] - istates[first20-1];
                  iind = first20 + (int)((invprime-istates[first20-1])/dif);
                } else {
                  dif = istates[1] - istates[0];
                  iind = 1 + (int)((invprime-istates[0])/dif);
		  }*/
                int vfindx1 = iind + nistates*(bprime+nbsmall*(pinds[n]-1));
                

                // uncomment if not using price states
                //vfindx1 = iind + nistates*bprime;
                //if(vffilled[vfindx1] != n) {
		int ok = iind < nistates;
		if(invmodel == 2) {
		  ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
		}
		if(ok) {
		  indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
		  if(vffilled[vfindx1] < 0) {
		    vffilled[vfindx1] = n;
		    if(gridmethod==2) {
		      vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
						    ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						    gvari.begin(), detiv(hh), dg, gmean.begin(),
						    nd-1, istates, nistates, first20, bprime, nbsmall, printflag,0,vfinterpmethod);

		      if(n==printobs) {
		      //  Rcpp::Rcout << "nreptilde " << nreptilde << std::endl;
			Rcpp::Rcout << "ivtprobs: " << *(ivtprobs+nrgrid*(pinds[n]-1)) << std::endl;
			Rcpp::Rcout << "kweights: " << kweights[0] << std::endl;
			Rcpp::Rcout << "vfin: " << vfin[0] << std::endl;
		        Rcpp::Rcout << "vf output (" << n << ", " << iind << "): " << vftemp[iind+nistates*bprime] << std::endl;
		      }

		    } else {
		      vftemp[vfindx1] = vfavgijchh(ivout, iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
						   ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
						   gvari.begin(), detiv(hh), dg, gmean.begin(),
						   nd-1, istates, nistates, first20, bprime, nbsmall, printflag);
		    }
		    if(n == printobs) {
		      Rcpp::Rcout << " hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
			" i': " << invprime << "b': " << bprime << " vf: " << vftemp[vfindx1] << std::endl;
		    }

		    //vfout(n+obsend*iind) = vftemp[vfindx1];

		  }
		}

                vfindx1 = iind-1 + nistates*(bprime+nbsmall*(pinds[n]-1));
                indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
                //vfindx1 = iind-1 + nistates*bprime;
                //if(vffilled[vfindx1] != n) {
                if(vffilled[vfindx1] < 0) {
                  vffilled[vfindx1] = n;
                  if(gridmethod == 2) {
                    vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
                                                  ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
                                                  gvari.begin(), detiv(hh), dg, gmean.begin(),
                                                  nd-1, istates, nistates, first20, bprime, nbsmall, 0, 0,vfinterpmethod);
                  } else {
                    vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
                                                 ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*hh,
                                                 gvari.begin(), detiv(hh), dg, gmean.begin(),
                                                 nd-1, istates, nistates, first20, bprime, nbsmall, 0);
                  }

                  if(n == printobs) {
                    Rcpp::Rcout << " hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
                      " i': " << invprime << "b': " << bprime << " vf: " << vftemp[vfindx1] << std::endl;
                  }

                  //vfout(n+obsend*(iind-1)) = vftemp[vfindx1];

                }

              }
            }
          }
        }
        //Rcpp::Rcout << "4" << std::endl;
        int first = 1;

        int vfoffset = nistates*nbsmall*(pinds[n]-1);

        for(int size=0;size<nsize;size++) {
          for(int j=0;j<capj+1;j++) {
            int indx = j+(1+capj)*size;
            drop[indx] = 1;
            if( (j == 0 && size == 0) || j > 0 ) {
              drop[indx] = 0;
	      if(invmodel == 1) {
		invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		invprime = invprime > invub ? invub : invprime;
		ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize(size), rcnew, gamma, nu, ccost, omega, cmodel, 0);
	      } else {
		nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                        &nextbflag,0);
		indcheck(bprime,0,nb,"bprime (sim) error");
		ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				   nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
		bchoice[indx] = bprime;
		bprime = bindex[bprime]-1;
		ichoice[indx] = invprime;
	      }
              if(!myopic && beta > 0) {
                vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20,vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
              } else {
                vfhat[indx] = 0.0;
              }
              vfout(n+obsend*(indx)) = vfhat[indx];
              utout(n+obsend*(indx)) = ut[indx];
              maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);

              if(n==printobs) {
                Rcpp::Rcout << " ivgrid[" << indx << "]: " << ivgrid[indx] << " ut[" << indx << "]: " << ut[indx] << " vfhat[" << indx << "]: " << vfhat[indx] << std::endl;
              }

            }
          }
        }
        //Rcpp::Rcout << "4.5" << std::endl;
        //if(printflag)
        //  Rprintf("\n");

        //ut[0] -= maxutil;
        //cerror[0] = -log(-log(R::runif(0,1)));
        double s1 = 0.0; //exp(ut[0]+vfhat[0]);

        double choiceumax = -1000000;  //ut[0]+vfhat[0] + cerror[0];
        for(int j=0;j<3;j++) {
          simchoice(n+nobsdata*j) = 0;
        }

        for(int size=0;size<nsize;size++) {
          for(int j=0;j<capj+1;j++) {
            int indx = j+(1+capj)*size;
            if(!drop[indx]) {
              ut[indx] -= maxutil;
              cerror[indx] = -log(-log(R::runif(0,1)));
              s1 += exp(ut[indx]+vfhat[indx]);
              if(ut[indx]+vfhat[indx]+cerror[indx] > choiceumax) {
                choiceumax = ut[indx]+vfhat[indx]+cerror[indx];
		if(invmodel==2) {
		  bprime = bchoice[indx];
		  invprime = ichoice[indx];
		}
                if(j > 0) {
                  simchoice[n] = size+1;
                  simchoice[n+nobsdata] = j;
                }
              }
            }
          }
        }
        //Rcpp::Rcout << "4.6" << std::endl;
        // simulate brand choice if a size choice > 0 is chosen

        if(simchoice(n) > 0) {
          choiceumax = -1000000;
          double ut1;
	  double s=0;
	  double utnoerr = 0;
          for(int j=0;j<nbrand;j++) {
            //Rcpp::Rcout << "column index " << prindex << ", " << j << ", " << n << std::endl;
            //double price = REAL(VECTOR_ELT(data,prindex-1+j))[n];
	    double price = pricematbr(n,j);
	    //Rcpp::Rcout << "y" << std::endl;
            if(!(price >= 999 || simchoice(n) != brsize(j+nbrand))) {
              //Rcpp::Rcout << "z" << std::endl;
              ut1 = xi[j] + alpha*price - log(-log(R::runif(0,1)));
              if(ut1 > choiceumax) {
                choiceumax = ut1;
		utnoerr = xi[j] + alpha*price;
                simchoice(n+2*nobsdata) = j+1;
              }
	      s+=exp(xi[j] + alpha*price);
            }
          }
	  llhhbr(hh) += log(exp(utnoerr)/s);
	  
        }
        //Rcpp::Rcout << "4.7" << std::endl;
        double obspack = 0;
        int obschoice = 0;

        if(simchoice[n] > 0) {
          obspack = packsize(simchoice(n)-1);
          obschoice = simchoice(n+nobsdata) + (1+capj)*(simchoice(n)-1);
        }

        indcheck(obschoice,0,(1+capj)*nsize,"ut. tunits: %d",((int)simchoice[n+nobsdata]));
        //if( dmax(inv + obspack*((double)simchoice[n+nobsdata]) - rcnew,0) > invub ) {
        //  llrun -= 1000000;
        //} else {
        //llrun += ut[obschoice] + vfhat[obschoice] - log(s1);
        //}

	llrun += log(exp(ut[obschoice] + vfhat[obschoice])/s1);

        if(debug) {
          //if( dmax(inv + obspack*((double)simchoice[n+nobsdata]) - rcnew,0) > invub ) {
          //  rllrun[n] = -1000000;
          //} else {
          rllrun(n) = exp(ut[obschoice]+vfhat[obschoice])/s1;
          //}
        }
        risave2(n) = inv;
	bsave2(n) = b;
	if(invmodel == 1) {
	  inv = dmax(inv + obspack*((double)simchoice[n+nobsdata]) - rcnew,0);
	  inv = inv > invub ? invub : inv;
	} else {
	  inv = invprime;
	  b = bprime;
	}
        n2++;
      }
      if(n==(nobs-1) || panid[n] != panid[n+1]) {
        llhh(hh) = llrun;  // I changed to return log likelihood rather than likelihood
        ll += log(llhh(hh));
      }
      //Rcpp::Rcout << "5" << std::endl;
    }
  }

  //Rcpp::Rcout << "done " << std::endl;
  if(retinv) {
    Rcpp::IntegerVector dim1(3), dim2(2);
    dim1[0] = nsim;
    dim1[1] = nhhs;
    dim1[2] = ninitt;
    dim2[0] = nsim;
    dim2[1] = nobs;
    risave1.attr("dim") = dim1;
    risave2.attr("dim") = dim2;
    invsave[0] = risave1;
    invsave[1] = risave2;

    return invsave;

  }

  if(debug) {
    Rcpp::List res(3);
    Rcpp::NumericVector lltot(1);
    lltot[0] = ll;

    Rcpp::IntegerVector dim1(2);
    dim1[0] = nobs;
    dim1[1] = nsim;

    rllrun.attr("dim") = dim1;

    res[0] = lltot;
    res[1] = llhh;
    res[2] = rllrun;

    return res;

  } else {

    int nres = 8;
    if(invmodel == 2) {
      nres = 13;
    }
    
    Rcpp::List res(nres);
    Rcpp::IntegerVector dim1(2);
    dim1[0] = nobsdata;
    dim1[1] = 3;

    simchoice.attr("dim") = dim1;
    res[0] = simchoice;
    res[1] = risave2;

    Rcpp::IntegerVector dim2(3);
    dim2[0] = obsend;
    dim2[1] = 1+capj;
    dim2[2] = nsize;
    vfout.attr("dim") = dim2;
    res[2] = vfout;

    utout.attr("dim") = dim2;
    res[3] = utout;

    Rcpp::IntegerVector dim3(2);
    dim3[0] = obsend;
    dim3[1] = nrgrid;
    tprobout.attr("dim") = dim3;
    res[4] = tprobout;

    res[5] = ichoiceout;

    res[6] = llhh;
    res[7] = llhhbr;
    if(invmodel==2) {
      res[8] = bsave2;
      res[9] = ichoicej;
      res[10] = ichoices;
      res[11] = risave1;
      res[12] = ibsave;
    }

    return res;
  }

  END_RCPP

}


// brand likelihood

// [[Rcpp::export]]
Rcpp::List brandlik(Rcpp::List inputs, Rcpp::List Mcmc, Rcpp::List info,
                    Rcpp::List hhmergebr) {
  BEGIN_RCPP

  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int crfix = Rcpp::as< std::vector<int> >(info["crfix"])[0];
  int datacrate = Rcpp::as< std::vector<int> >(info["data.crate"])[0];
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }
  
  Rcpp::NumericVector crate(nhhs);
  if(datacrate) {
    SEXP temp(info["crate"]);
    Rcpp::NumericVector temp1(temp);
    crate = clone(temp1);
  }


  SEXP xfull1 = info["xfull"];
  Rcpp::NumericVector xfull(xfull1);

  int npbig = xfull.size();

  SEXP tf1 = info["tform"];
  Rcpp::NumericVector tf(tf1);
  SEXP fixed1 = info["fixed"];
  Rcpp::IntegerVector fixed(fixed1);
  SEXP pequal1 = info["paramequal"];
  Rcpp::IntegerVector paramequal(pequal1);
  SEXP lbounds1 = info["lbounds"];
  Rcpp::NumericVector lbounds(lbounds1);
  SEXP ubounds1 = info["ubounds"];
  Rcpp::NumericVector ubounds(ubounds1);
  int crhhlb = Rcpp::as< std::vector<int> >(info["crhhlb"])[0];

  if(!inputs.containsElementNamed("pstart"))
    throw std::range_error("pstart is NULL in inputs");

  SEXP inputpstart = inputs["pstart"];
  Rcpp::NumericMatrix pstart1(inputpstart);
  Rcpp::NumericMatrix pstart = Rcpp::clone(inputpstart);
  int npsmall = 0;
  for(int i = 0;i < npbig;i++) {
    npsmall += !fixed(i);
  }

  int useparamstart = 0;
  Rcpp::NumericMatrix paramstart(npbig,nhhs);
  if(info.containsElementNamed("paramstart")) {
    SEXP paramstart1 = info["paramstart"];
    Rcpp::NumericMatrix paramstart1a(paramstart1);
    pstart = Rcpp::clone(paramstart1a);
    paramstart = Rcpp::clone(paramstart1a);
    useparamstart = 1;
  }

  if(datacrate) {
    pstart(nbrand,Rcpp::_) = crate;
  }

  Rcpp::NumericMatrix pold(npsmall,nhhs);
  int j=-1;
  for(int i=0;i<npbig;i++) {
    if(!fixed(i)) {
      j++;
      pold(j,Rcpp::_) = pstart(i,Rcpp::_);
    }
  }

  Rcpp::NumericMatrix poldtf(npbig,nhhs);

  SEXP tunitstemp(hhmergebr["totunits"]);
  Rcpp::NumericVector tunitsbr(tunitstemp);

  SEXP panidtemp(hhmergebr["PANID"]);
  Rcpp::NumericVector panidbr(panidtemp);

  SEXP brindextemp(hhmergebr["brindex"]);
  Rcpp::IntegerVector brindexbr(brindextemp);

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);


  int nobsbr = tunitsbr.size();

  int varflag = 0;
  if(info.containsElementNamed("varflag"))
    varflag = Rcpp::as< std::vector<int> >(info["varflag"])[0];

  int pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  int retindll = Rcpp::as< std::vector<int> >(info["retindll"])[0];

  int prindex = Rcpp::as< std::vector<int> >(info["prindex"])[0];
  Rcpp::NumericMatrix pricematbr(nobsbr,nbrand);
  SEXP ptemp = hhmergebr[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = hhmergebr[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }

  int llbrandsize = nhhs;
  if(varflag)
    llbrandsize = nobsbr;

  // declarations of variables that get returned/work variables

  Rcpp::NumericVector llbrand(llbrandsize);
  //Rcpp::Rcout << "bbb" << std::endl;
  tform(pold,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
        nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
        sizebrand.begin(),poldtf);

  /*Rcpp::Rcout << "pold: ";
  for(int i=0;i<npsmall;i++) {
    Rcpp::Rcout << " " << pold(i,0);
  }
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << "poldtf: ";
  for(int i=0;i<npbig;i++) {
    Rcpp::Rcout << " " << poldtf(i,0);
  }
  Rcpp::Rcout << std::endl;*/
  //Rcpp::Rcout << "ccc" << std::endl;
  brchoicellhh(poldtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
               nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand);

  Rcpp::List ret;
  ret["llbrand"] = llbrand;

  return ret;

  END_RCPP

}



// estimation using MLE - this is just for the no heterogeneity version of the model
// make sure to set nsave, ntilde etc to 1

// [[Rcpp::export]]
Rcpp::List mlelik(Rcpp::List inputs, Rcpp::List Mcmc, Rcpp::List info,
		  Rcpp::List hhbig, Rcpp::List hhmerge, Rcpp::List hhmergebr) {

  BEGIN_RCPP

#ifdef defined __unix__

  const rlim_t kStackSize = 48 * 1024 * 1024;   // min stack size = 16 MB
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
    {
      if (rl.rlim_cur < kStackSize)
        {
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0)
            {
              Rcpp::Rcout << "setrlimit returned result = " << result << std::endl;
              //fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }

#endif

  int maxrep = Rcpp::as< std::vector<int> >(Mcmc["R"])[0];
  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];
  int crfix = Rcpp::as< std::vector<int> >(info["crfix"])[0];
  int datacrate = Rcpp::as< std::vector<int> >(info["data.crate"])[0];
  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  int usecdraws = Rcpp::as< std::vector<int> >(info["usecdraws"])[0];
  double tol = Rcpp::as< std::vector<double> >(info["tol"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }

  Rcpp::NumericVector crate(nhhs);
  if(datacrate) {
    SEXP temp(info["crate"]);
    Rcpp::NumericVector temp1(temp);
    crate = clone(temp1);
  }

  SEXP xfull1 = info["xfull"];
  Rcpp::NumericVector xfull(xfull1);

  //Rcpp::Rcout << "maxrep: " << maxrep << std::endl;
  //Rcpp::Rcout << "nhhs: " << nhhs << std::endl;
  //Rcpp::Rcout << "xfull: ";

  int npbig = xfull.size();
  //for(int i=0;i<npbig;i++) {
  //  Rcpp::Rcout << xfull[i] << " ";
  //}
  //Rcpp::Rcout << std::endl;

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdrop2(ivdrop1);
    ivdrop = clone(ivdrop2);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = clone(ivdrop);


  SEXP tf1 = info["tform"];
  Rcpp::NumericVector tf(tf1);
  SEXP fixed1 = info["fixed"];
  Rcpp::IntegerVector fixed(fixed1);
  SEXP pequal1 = info["paramequal"];
  Rcpp::IntegerVector paramequal(pequal1);
  SEXP lbounds1 = info["lbounds"];
  Rcpp::NumericVector lbounds(lbounds1);
  SEXP ubounds1 = info["ubounds"];
  Rcpp::NumericVector ubounds(ubounds1);
  int crhhlb = Rcpp::as< std::vector<int> >(info["crhhlb"])[0];

  if(!inputs.containsElementNamed("pstart"))
    throw std::range_error("pstart is NULL in inputs");

  SEXP inputpstart = inputs["pstart"];
  Rcpp::NumericMatrix pstart1(inputpstart);
  Rcpp::NumericMatrix pstart = Rcpp::clone(inputpstart);
  int npsmall = 0;
  for(int i = 0;i < npbig;i++) {
    npsmall += !fixed(i);
  }

  int useparamstart = 0;
  Rcpp::NumericMatrix paramstart(npbig,nhhs);
  if(info.containsElementNamed("paramstart")) {
    SEXP paramstart1 = info["paramstart"];
    Rcpp::NumericMatrix paramstart1a(paramstart1);
    pstart = Rcpp::clone(paramstart1a);
    paramstart = Rcpp::clone(paramstart1a);
    useparamstart = 1;
  }

  if(datacrate) {
    pstart(nbrand,Rcpp::_) = crate;
  }

  Rcpp::NumericMatrix pold(npsmall,nhhs);
  int j=-1;
  for(int i=0;i<npbig;i++) {
    if(!fixed(i)) {
      j++;
      pold(j,Rcpp::_) = pstart(i,Rcpp::_);
    }
  }

  Rcpp::NumericMatrix poldtf(npbig,nhhs);

  // this will not copy the object

  SEXP tunitstemp(hhmergebr["totunits"]);
  Rcpp::NumericVector tunitsbr(tunitstemp);

  SEXP panidtemp(hhmergebr["PANID"]);
  Rcpp::NumericVector panidbr(panidtemp);

  SEXP brindextemp(hhmergebr["brindex"]);
  Rcpp::IntegerVector brindexbr(brindextemp);

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  int pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  int retindll = Rcpp::as< std::vector<int> >(info["retindll"])[0];
  int varflag = 0;
  if(info.containsElementNamed("varflag"))
    varflag = Rcpp::as< std::vector<int> >(info["varflag"])[0];

  SEXP panidmtemp(hhmerge["PANID"]);
  Rcpp::NumericVector panidmerge(panidmtemp);

  Rcpp::IntegerVector hhinds(2);
  hhinds(0) = 1;
  hhinds(1) = nhhs;

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  SEXP pindsmerge1(info["idpriceindmerge"]);
  Rcpp::IntegerVector pindsmerge(pindsmerge1);

  SEXP pindsbig1(info["idpriceindbig"]);
  Rcpp::IntegerVector pindsbig(pindsbig1);

  int hmaxpind = Rcpp::as< std::vector<int> >(info["hmaxpind"])[0];
  SEXP hpindsbig1(info["hpriceindbig"]);
  Rcpp::IntegerVector hpindsbig(hpindsbig1);

  int maxpind = Rcpp::as< std::vector<int> >(info["maxpind"])[0];
  //Rcpp::Rcout << "maxpind: " << maxpind << std::endl;
  int nobsbr = tunitsbr.size();
  int nobsmerge = panidmerge.size();

  int prindex = Rcpp::as< std::vector<int> >(info["prindex"])[0];
  Rcpp::NumericMatrix pricematbr(nobsbr,nbrand);
  SEXP ptemp = hhmergebr[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = hhmergebr[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }

  // this copies - should just work with numericmatrix directly
  //std::vector<double> pricevecbr = Rcpp::as< std::vector<double> >(pricematbr);

  Rcpp::NumericMatrix pricematmerge(nobsmerge,nbrand);
  SEXP ptemp2 = hhmerge[prindex-1];
  Rcpp::NumericVector ptemp3(ptemp2);
  for(int i=0;i<nbrand;i++) {
    ptemp2 = hhmerge[prindex-1+i];
    ptemp3 = ptemp2;
    pricematmerge(Rcpp::_,i) = ptemp3;
  }

  // define inputs for quantity log likelihood
  llinputs llinfo;

  SEXP hhbigcol1 = hhbig[0];
  Rcpp::NumericVector hhbigc1(hhbigcol1);
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];

  llinfo.nb = nb;
  llinfo.nbsmall = nbsmall;
  llinfo.nco = npbig;
  llinfo.nobs = hhbigc1.size();
  llinfo.nbrand = nbrand;
  llinfo.nsize = nsize;
  llinfo.nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  llinfo.sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  llinfo.ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  llinfo.inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];
  int nsave = 10;
  if(Mcmc.containsElementNamed("nsave")) {
    llinfo.nsave = Rcpp::as< std::vector<int> >(Mcmc["nsave"])[0];
    nsave = llinfo.nsave;
  } else {
    llinfo.nsave = 10;
  }
  if(Mcmc.containsElementNamed("ntilde")) {
    llinfo.ntilde = Rcpp::as< std::vector<int> >(Mcmc["ntilde"])[0];
  } else {
    llinfo.ntilde = 3;
  }

  int nitervf = 50; // number of iterations for value function update
  int niterstop = 500; // rep at which we switch to only a single update
  if(Mcmc.containsElementNamed("nitervf")) {
    nitervf = Rcpp::as< std::vector<int> >(Mcmc["nitervf"])[0];
  }

  if(Mcmc.containsElementNamed("niterstop")) {
    niterstop = Rcpp::as< std::vector<int> >(Mcmc["niterstop"])[0];
  }

  llinfo.nrgrid = Rcpp::as< std::vector<int> >(info["nrgrid"])[0];
  llinfo.necoef = Rcpp::as< std::vector<int> >(info["necoef"])[0];
  llinfo.capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];
  llinfo.nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];
  llinfo.ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];
  llinfo.retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];
  llinfo.myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];
  llinfo.cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];
  llinfo.debug = Rcpp::as< std::vector<int> >(info["debug"])[0];
  llinfo.usevf = 1;
  llinfo.idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];
  llinfo.first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];
  llinfo.hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];
  llinfo.genrnd = Rcpp::as< std::vector<int> >(info["genrnd"])[0];
  llinfo.nd = Rf_length(info["ngrid"]);
  llinfo.nq = INTEGER(Rf_getAttrib(info["qpts"], R_DimSymbol))[0];
  llinfo.ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];
  llinfo.pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  llinfo.hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  llinfo.repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0]-1;
  llinfo.ncpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[1];
  llinfo.rcpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[0];
  llinfo.ncrate = Rcpp::as< std::vector<int> >(info["ncost"])[0];
  llinfo.dg = Rcpp::as< std::vector<double> >(info["dg"])[0];
  llinfo.maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];
  llinfo.hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];
  llinfo.usecdraws = usecdraws;
  llinfo.initvf = 0;
  
  llinfo.invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  llinfo.maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];

  if(llinfo.invmodel == 1) {
    llinfo.lomega = 2;
  } else {
    llinfo.lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }
  
  int hmodel = llinfo.hmodel;

  Rcpp::NumericVector vfll(llinfo.nobs*nsize*(1+llinfo.capj));
  Rcpp::NumericVector utll(llinfo.nobs*nsize*(1+llinfo.capj));


  if(llinfo.debug && llinfo.retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  Rcpp::NumericVector cosave(npbig*nhhs*llinfo.nsave);
  Rcpp::IntegerVector indexes(llinfo.nsave);

  for(int i=0;i<llinfo.nsave;i++) {
    indexes(i) = i+1;
  }

  SEXP tunitsbigtemp(hhbig["totunits"]);
  Rcpp::NumericVector tunitsbig(tunitsbigtemp);

  SEXP panidbigtemp(hhbig["PANID"]);
  Rcpp::NumericVector panidbig(panidbigtemp);

  SEXP brindexbigtemp(hhbig["brindex"]);
  Rcpp::IntegerVector brindexbig(brindexbigtemp);

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP expandbig1(info["expandbig"]);
  Rcpp::IntegerVector expandbig(expandbig1);

  Rcpp::NumericMatrix ivbig(llinfo.nobs,nsize);
  Rcpp::NumericMatrix ivbig1(llinfo.nobs,nsize);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);

  int maxgridlength = ngrid[ngrid.size()-1];

  int vfblocksize1 = llinfo.nsave*llinfo.nrgrid*nbsmall*maxgridlength;
  int vfblocksize2 = llinfo.nrgrid*nbsmall*maxgridlength;

  Rcpp::NumericVector vfsave(vfblocksize1*nhhs);

  if(info.containsElementNamed("vfstart")) {
    SEXP tempvfinfo(info["vfstart"]);
    Rcpp::NumericVector tempvf(tempvfinfo);
    vfsave = Rcpp::clone(tempvf);
    llinfo.initvf = 1;
  }

  // if initvf is true, copy in supplied cosave
  if(llinfo.initvf) {
    SEXP tempcosave1(info["cosave"]);
    Rcpp::NumericVector tempcosave(tempcosave1);
    cosave = Rcpp::clone(tempcosave);
  }

  //Rcpp::Rcout << "aaa" << std::endl;

  Rcpp::NumericVector vf(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  //SEXP cdraws1(info["cdraws"]);
  //Rcpp::NumericVector cdraws(cdraws1);

  // cdraws needs to be reddrawn so we redefine it here
  Rcpp::NumericVector cdraws(llinfo.nobs);

  // icdraws won't be used anymore except for debugging
  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);

  SEXP cdrawtemp(info["cdraws"]);
  Rcpp::NumericVector cdrawinfo(cdrawtemp);

  SEXP iiv1(info["initivdraws"]);  // this and logitdraws can probably be removed
  Rcpp::NumericVector iiv(iiv1);

  SEXP ldraws1(info["logitdraws"]);
  Rcpp::NumericVector ldraws(ldraws1);

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);

  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);
  
  SEXP initb1(info["initb"]);
  Rcpp::IntegerVector initb(initb1);

  SEXP inits1(info["inits"]);
  Rcpp::IntegerVector inits(inits1);

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  SEXP qpts1(info["qpts"]);
  Rcpp::NumericVector qpts(qpts1);

  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);

  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);

  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);
  
  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);

  //std::vector<double> pricevecmerge = Rcpp::as< std::vector<double> >(pricematmerge);

  int llbrandsize = nhhs;
  if(varflag)
    llbrandsize = nobsbr;

  // declarations of variables that get returned/work variables

  Rcpp::NumericVector llbrand(llbrandsize);
  Rcpp::NumericVector llbrand1(llbrandsize);
  Rcpp::NumericVector llbrand1a(llbrandsize);

  Rcpp::NumericMatrix rivout(nobsmerge,nsize);
  Rcpp::NumericMatrix rivout1(nobsmerge,nsize);

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;
  Rcpp::NumericVector rcutmatout(cutblocksize*nhhs);
  Rcpp::NumericVector rcutmatout1(cutblocksize*nhhs);

  int nrowiv=0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivcoefblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  Rcpp::NumericVector rivcoefout(ivcoefblocksize*nhhs);
  Rcpp::NumericVector rivcoefout1(ivcoefblocksize*nhhs);

  int ivvariblocksize = nsize*nsize;
  Rcpp::NumericVector rivvariout(ivvariblocksize*nhhs);
  Rcpp::NumericVector rivvariout1(ivvariblocksize*nhhs);

  Rcpp::NumericVector rdetivout(nhhs);
  Rcpp::NumericVector rdetivout1(nhhs);

  // setup of grids, starting values for B and W,

  Rcpp::IntegerVector fixedpop(npsmall);
  if(Mcmc.containsElementNamed("fixedcoef")) {
    SEXP tempfixed(Mcmc["fixedcoef"]);
    Rcpp::IntegerVector tempfixed1(tempfixed);
    fixedpop = clone(tempfixed1);
  }

  int nfixedpop = 0;
  for(int i=0;i<npsmall;i++) {
    nfixedpop += fixedpop(i);
  }

  int allfixed = nfixedpop == npsmall;
  int allvarying = nfixedpop == 0;

  if(llinfo.hmodel && !allfixed)
    throw std::runtime_error("if hmodel is true all parameters must be fixed");


  int fixedmbound = nfixedpop;
  if(allvarying)
    fixedmbound = 1;

  arma::mat proposal(fixedmbound,fixedmbound);
  proposal.eye();
  if(Mcmc.containsElementNamed("proposal")) {
    SEXP tempprop1(Mcmc["fixedcoef"]);
    Rcpp::NumericVector tempprop(tempprop1);
    if(tempprop.length() != nfixedpop*nfixedpop)
      throw std::range_error("Incorrect dimensions for proposal.  Make sure the dimensions correspond to the number of population fixed coefficient.");
    std::copy(tempprop.begin(),tempprop.end(),proposal.begin());
  }

  arma::mat cholprop = chol(proposal);

  // split up the proposal into 2 pieces

  int nbrfixed = 0;
  int ii1 = 0;
  for(int i=0;i<nbrand;i++) {
    if(!fixed(i)) {
      if(fixedpop(ii1)) {
        nbrfixed++;
      }
      ii1++;
    }
  }

  //Rcpp::Rcout << "nbrfixed: " << nbrfixed << " nfixedpop: " << nfixedpop << std::endl;

  int bend1 = nbrfixed;
  int bstart2 = nbrfixed;
  int bend2 = nfixedpop;

  if(nbrfixed == 0 || allvarying)
    bend1 = 1;

  arma::mat cpropbrand = cholprop.submat(0,0,bend1-1,bend1-1);

  if( nfixedpop == nbrfixed ) {
    bstart2 = nfixedpop-1;
    bend2 = nfixedpop;
  }

  if(allvarying) {
    bstart2 = 0;
    bend2 = 1;
  }

  arma::mat cpropother = cholprop.submat(bstart2,bstart2,bend2-1,bend2-1);

  arma::mat gvari1a(gvari.begin(),nsize,nsize,false);
  arma::mat gvar = arma::inv_sympd(gvari1a);
  arma::mat gchol = arma::chol(gvar);

  int nrgrid = llinfo.nrgrid;

  // gridmethod == 1: old random grid of inclusive values
  // gridmethod == 2: fixed grid of price vectors (randomly drawn in advance)
  int gridmethod = 1;
  if(Mcmc.containsElementNamed("gridmethod")) {
    gridmethod = Rcpp::as< std::vector<int> >(Mcmc["gridmethod"])[0];
  }

  int vfinterpmethod = 1;
  if(Mcmc.containsElementNamed("vfinterpmethod")) {
    vfinterpmethod = Rcpp::as< std::vector<int> >(Mcmc["vfinterpmethod"])[0];
  }

  Rcpp::NumericMatrix bwmatrix(npbig,npbig);

  double hthumb = pow(4/(3*( (double)nsave )),0.2);
  for(int i=0;i<npbig;i++) {
    bwmatrix(i,i) = 1/hthumb;
  }

  if(Mcmc.containsElementNamed("bwmatrix")) {
    SEXP bwtemp1(Mcmc["bwmatrix"]);
    Rcpp::NumericMatrix bwtemp(bwtemp1);
    bwmatrix = Rcpp::clone(bwtemp);
  }

  int rgridbdbig = gridmethod == 2 ? nbrand*llinfo.nrgrid*llinfo.nsave : nsize*llinfo.nrgrid*llinfo.nsave;

  int rgridbd = gridmethod == 2 ? nbrand : nsize;

  Rcpp::NumericMatrix rgrid(rgridbd,nrgrid);
  Rcpp::NumericVector rgridsave(rgridbdbig);


  if(gridmethod == 1) {
    for(int i=0;i<nsize;i++) {
      rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
    }
  } else if (gridmethod == 2) {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  Rcpp::NumericVector impdist(nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }

  arma::mat rgrida(rgrid.begin(),rgridbd,nrgrid,false);
  arma::colvec gmeana(gmean.begin(),gmean.size(),false);

  arma::mat tempm(nsize,nsize);
  if(gridmethod != 2) {
    tempm = rgrida.t()*gchol;
    rgrida = tempm.t();
    rgrida.each_col() += gmeana;
  }

  int nistates = ngrid[llinfo.nd-1];

  int nkeep = Rcpp::as< std::vector<int> >(Mcmc["keep"])[0];

  int ndrawkeep = maxrep/nkeep;

  Rcpp::NumericMatrix bdraw(ndrawkeep,npsmall);

  Rcpp::NumericMatrix Wdraw(ndrawkeep,npsmall*npsmall);

  Rcpp::NumericVector llbrsave(ndrawkeep);
  Rcpp::NumericVector llqsave(ndrawkeep);

  // make these TRANSFORMED beta draws
  Rcpp::NumericVector betadraw(ndrawkeep*nhhs*npsmall);

  Rcpp::NumericMatrix W(npsmall,npsmall);


  for(int i=0;i<npsmall;i++) {
    if(!fixedpop[i])
      W(i,i) = 1;
  }

  Rcpp::NumericVector b(npsmall);

  for(int i=0;i<npsmall;i++) {
    b[i] = Rcpp::mean(pold(i,Rcpp::_));
  }

  arma::mat W1(npsmall-nfixedpop,npsmall-nfixedpop);
  W1.zeros();

  arma::mat W1chol(npsmall-nfixedpop,npsmall-nfixedpop);
  W1chol.zeros();

  arma::mat W1inv(npsmall-nfixedpop,npsmall-nfixedpop);
  W1inv.zeros();

  arma::colvec bnotfixed(npsmall-nfixedpop);
  bnotfixed.zeros();

  arma::mat temppnew(npsmall-nfixedpop,nhhs);
  temppnew.zeros();

  arma::mat temppold(npsmall-nfixedpop,nhhs);
  temppold.zeros();

  Rcpp::NumericVector bdiff0(nhhs);
  Rcpp::NumericVector bdiff1(nhhs);

  arma::mat bdiff0a(bdiff0.begin(),bdiff0.size(),false);
  arma::mat bdiff1a(bdiff1.begin(),bdiff1.size(),false);

  arma::mat betamean(npsmall-nfixedpop,1);
  betamean.zeros();

  arma::mat b1(npsmall-nfixedpop,1);
  b1.zeros();

  arma::mat b2(nfixedpop,1);
  b2.zeros();

  arma::mat b2a(nbrfixed,1);
  b2a.zeros();

  arma::mat b2b(nfixedpop-nbrfixed,1);
  b2b.zeros();
  //Rcpp::Rcout << "b2b size " << nfixedpop-nbrfixed << std::endl;

  arma::mat btempdiff(npsmall-nfixedpop,1);
  btempdiff.zeros();

  Rcpp::NumericVector pnewdrawvec((npsmall-nfixedpop)*nhhs);
  arma::mat pnewdrawmat(pnewdrawvec.begin(),npsmall-nfixedpop,nhhs,false);

  Rcpp::NumericVector llqhh(nhhs);
  Rcpp::NumericVector llqhh1(nhhs);
  Rcpp::NumericVector llqhh1a(nhhs);

  int nllrun = llinfo.debug == 2 ? llinfo.nobs*(1+llinfo.capj)*nsize : llinfo.nsim*llinfo.nobs;
  
  Rcpp::NumericVector rllrun(nllrun);
  Rcpp::NumericVector risave1(llinfo.nsim*nhhs*llinfo.ninitt);
  Rcpp::NumericVector risave2(llinfo.nsim*llinfo.nobs);

  Rcpp::NumericMatrix pnew(npsmall,nhhs);
  Rcpp::NumericMatrix pnewtf(npbig,nhhs);

  Rcpp::NumericVector r(nhhs);
  Rcpp::NumericVector u(nhhs);

  arma::mat npsmallnf(npsmall-nfixedpop,npsmall-nfixedpop);
  npsmallnf.eye();
  npsmallnf *= npsmall-nfixedpop;

  arma::mat S(npsmall-nfixedpop,npsmall-nfixedpop);
  S.zeros();

  arma::mat Schol(npsmall-nfixedpop,npsmall-nfixedpop);
  Schol.zeros();

  arma::mat XX(npsmall-nfixedpop,npsmall-nfixedpop);
  XX.zeros();

  arma::mat Stemp(npsmall-nfixedpop,npsmall-nfixedpop+nhhs);
  //arma::mat Stemp(npsmall-nfixedpop,nu+nhhinclude);
  Stemp.zeros();

  arma::mat tempb1(1,npsmall - nfixedpop);
  tempb1.zeros();

  arma::mat tempb2(1,nfixedpop);
  tempb2.zeros();

  arma::mat tempb2a(1,nbrfixed);
  tempb2a.zeros();

  arma::mat tempb2b(1,nfixedpop-nbrfixed);
  tempb2b.zeros();

  // varying coefficient rho (1st el is brands, and second is other params)
  Rcpp::NumericVector rho(2);
  if(Mcmc.containsElementNamed("rho")) {
    SEXP temprho(Mcmc["rho"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho in Mcmc must be 2");
    rho = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho(i) = 0.1;
  }

  // fixed coefficients rho1
  Rcpp::NumericVector rho1(2);
  if(Mcmc.containsElementNamed("rho1")) {
    SEXP temprho(Mcmc["rho1"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho1 in Mcmc must be 2");
    rho1 = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho1(i) = 0.1;
  }

  // saves of rho1 - so we can see when it stabilizes
  Rcpp::NumericMatrix rho1save(ndrawkeep,2);

  int splitfixed = 1;
  if(Mcmc.containsElementNamed("splitfixed")) {
    splitfixed = Rcpp::as< std::vector<int> >(Mcmc["splitfixed"])[0];
  }


  int npropsave = 100;
  if(Mcmc.containsElementNamed("npropsave")) {
    npropsave = Rcpp::as< std::vector<int> >(Mcmc["npropsave"])[0];
  }

  int nprint = 10;
  if(Mcmc.containsElementNamed("nprint")) {
    nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  }

  std::vector<std::string> pnamevec(npbig,"V");
  for(int i=0;i<npbig;i++) {
    std::stringstream ss;
    ss << i+1;
    pnamevec[i] += ss.str();
  }

  Rcpp::CharacterVector pnames = Rcpp::wrap(pnamevec);

  if(Mcmc.containsElementNamed("pnames")) {
    nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  }

  Rcpp::IntegerMatrix xaccept(npropsave,2);
  Rcpp::IntegerVector acceptflag(nhhs);

  int wmethod = 2;  //inverse wishart drawing method - all seem about the same
  // 1: Rossi code (with demographics)
  // 2: my code (derived from ken train's R code)
  // 3: inverse gamma

  // priors

  arma::mat S0 = arma::eye(npsmall-nfixedpop,npsmall-nfixedpop);
  int nu = npsmall-nfixedpop+3;
  S0 = S0*( (double)nu );
  if(Mcmc.containsElementNamed("S0")) {
    SEXP temp(Mcmc["S0"]);
    Rcpp::NumericMatrix tempS(temp);
    for(int i=0;i<npsmall-nfixedpop;i++) {
      for(int j=0;j<npsmall-nfixedpop;j++) {
        S0(i,j) = tempS(i,j);
      }
    }
  }
  if(Mcmc.containsElementNamed("nu")) {
    nu = Rcpp::as< std::vector<int> >(Mcmc["nu"])[0];
  }

  nu += nhhs;


  int naccept = 0;

  //Rcpp::Rcout << "done setting stuff" <<std::endl;


  int keep=0;

  int rep = 0;

  if(usecdraws) {
    cdraws = Rcpp::clone(cdrawinfo);
  } else {
    cdraws = Rcpp::runif(llinfo.nobs);
  }

  tform(pold,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
        nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
        sizebrand.begin(),poldtf);

  brchoicellhh(poldtf,tunitsbr,panidbr,pricematbr,brindexbr,brsize,nobsbr,
               nbrand,pflag,retindll,nhhs,varflag,npbig,llbrand);

  if(hmodel) {
    fitiv(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
          obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
          ncutiv, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout);
  } else {
    fitivhh(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
            obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
            ncutiv, dropmat, ivdrop, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout,
	    ivdrop2);
  }

  //Rcpp::Rcout << "abc" << std::endl;

  int stoprep = 10000;
  
  if(!llinfo.myopic) {

    double dif = tol + 1;

    while(dif > tol) {

      Rcpp::NumericVector vfold = Rcpp::clone(vf);

      nitervf = 1;

      int indx1 = rep%nsave;

      vfupdateijc(hhinds,bstates,pgrid,ngrid,poldtf,vfsave,rgrid,cosave,rgridsave,
                  packsize,rcutmatout,rivcoefout,rivvariout,rdetivout,gvari,gmean,
                  indexes,brsize,impdist,bwmatrix,revarray,badmissable,bindex,ivdrop2,
		  llinfo,rep,nitervf,gridmethod,vfinterpmethod,vf);

      // see if can speed this up w std::copy - will need to reindex
      for(int hh=0;hh<nhhs;hh++) {
        for(int i=0;i<vfblocksize2;i++) {
          indcheck(hh*vfblocksize1 + indx1 + nsave*i,0,vfsave.length(),"vfsave copying bug: hh: %d; rep: %d; i: %d; vfblocksize1: %d; nsave: %d",hh,rep,i,
                   vfblocksize1,nsave);
          indcheck(hh*vfblocksize2+i,0,vf.length(),"vf copying bug: hh: %d; rep: %d; i: %d",hh,rep,i);
          vfsave[hh*vfblocksize1 + indx1 + nsave*i] = vf[hh*vfblocksize2+i];
        }
      }

      Rcpp::NumericVector vfdiff = Rcpp::abs(vfold-vf);

      int mvfi;
      maxvec(vfdiff.begin(),vfdiff.size(),&dif,&mvfi,0);
      //Rcpp::Rcout << "Rep: " << rep << "; Mean vf: " << Rcpp::mean(vf) << "dif: "
      //            << dif << std::endl;
      rep++;

      if(rep >= stoprep) {
	break;
      }
      
    }

  }

  for(int i=0;i<llinfo.nobs;i++) {
    for(int j=0;j<nsize;j++) {
      ivbig(i,j) = rivout(expandbig[i]-1,j);
    }
  }
  Rcpp::Rcout << "computing the likelihood" << std::endl;
  if(hmodel) {
    llijc(poldtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
          hhinds,brsize,obshhinds,ivbig,vfsave,cdraws,ic,iiv,ldraws,
          packsize,initx,initb,inits,ngrid,rcutmatout,rivcoefout,rivvariout,rdetivout,
          pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,hpindsbig,
	  bstates,revarray,badmissable,bindex,ivdrop,
          llinfo,rep,gridmethod,vfinterpmethod,hmaxpind,
          vfll,utll,rllrun,risave1,risave2,llqhh);
  } else {
    llijc(poldtf,cosave,rgridsave,indexes,tunitsbig,panidbig,brindexbig,
          hhinds,brsize,obshhinds,ivbig,vfsave,cdraws,ic,iiv,ldraws,
          packsize,initx,initb,inits,ngrid,rcutmatout,rivcoefout,rivvariout,rdetivout,
          pgrid,qpts,gvari,gmean,rgrid,impdist,bwmatrix,pindsbig,
	  bstates,revarray,badmissable,bindex,ivdrop2,
          llinfo,rep,gridmethod,vfinterpmethod,maxpind,
          vfll,utll,rllrun,risave1,risave2,llqhh);
  }





  Rcpp::List ret;
  ret["poldtf"] = poldtf;
  ret["pnewtf"] = pnewtf;
  ret["llbrand"] = llbrand;
  ret["riv"] = rivout;
  ret["rivcoef"] = rivcoefout;
  ret["rivvari"] = rivvariout;
  ret["rdetiv"] = rdetivout;
  ret["llqhh"] = llqhh;
  ret["vf"] = vf;
  ret["pstart"] = pstart;

  Rcpp::NumericVector vfout(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  for(int hh=0;hh<nhhs;hh++) {
    for(int j=0;j<nrgrid;j++) {
      for(int i=0;i<nistates;i++) {
        for(int b=0;b<nbsmall;b++) {
          vfout(hh + nhhs*(j+nrgrid*(i+nistates*b))) = vf(hh*vfblocksize2 + j+nrgrid*(i+nistates*b));
        }
      }
    }
  }

  Rcpp::IntegerVector dims(4);
  dims(0) = nhhs;
  dims(1) = nrgrid;
  dims(2) = nistates;
  dims(3) = nbsmall;

  vfout.attr("dim") = dims;

  ret["vfsmall"] = vfout;

  Rcpp::IntegerVector dim2(3);
  dim2[0] = llinfo.nobs;
  dim2[1] = 1+llinfo.capj;
  dim2[2] = nsize;

  vfll.attr("dim") = dim2;
  utll.attr("dim") = dim2;

  ret["vfll"] = vfll;
  ret["utll"] = utll;

  ret["risave2"] = risave2;

  if(llinfo.debug == 2) {
    rllrun.attr("dim") = dim2;
    ret["llbig"] = rllrun;
  }

  return ret;


  END_RCPP

}


// check inventory states

// [[Rcpp::export]]
Rcpp::List checkinvstate(Rcpp::NumericVector hhcrate, Rcpp::List data, Rcpp::List info,
			 Rcpp::IntegerVector hhinds, Rcpp::List hhinitmerge,
			 Rcpp::IntegerVector stvisitinit) {

  BEGIN_RCPP
  
    
  int nobs = Rf_length(VECTOR_ELT(data,0));
  
  int initnobs = Rf_length(VECTOR_ELT(hhinitmerge,0));

  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];

  int maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];

  SEXP sizechoiceinit1 = hhinitmerge["sizechoice"];
  Rcpp::IntegerVector sizechoiceinit(sizechoiceinit1);

  SEXP nbottlesinit1 = hhinitmerge["nbottles"];
  Rcpp::IntegerVector nbottlesinit(nbottlesinit1);

  SEXP tunitstemp(data["totunits"]);
  Rcpp::NumericVector tunits(tunitstemp);

  SEXP brindex1(data["brindex"]);
  Rcpp::NumericVector brindex(brindex1);

  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);
  
  int capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];
  
  int nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];

  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];

  int ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);
  
  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);
  
  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);
  
  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);
  int linitx = initx.length();
  
  int bprime = 0;
  double invprime = 0;
  
  double inv = 0;
  int b=0;

  int iobs = -1;
  int iobsbig = 0;

  int nhh = hhinds(1)-hhinds(0)+1;
  int ninitobsbig = nhh*ninitt;

  Rcpp::NumericVector iinvout(ninitobsbig);
  Rcpp::IntegerVector ibout(ninitobsbig);

  double initinv[nhhs];
  double initb[nhhs];

  int nextbflag=0;
  
  for(int n=hhinds(0)-1;n<hhinds(1);n++) {
    inv = 0;
    b = 0;
    double rcnew = hhcrate(n);
    for(int t=0;t<ninitt;t++) {
      if(stvisitinit(iobsbig)) {
	iobs++;
	nextinv(nbottlesinit(iobs),sizechoiceinit(iobs)-1,b,inv,bstates.begin(),
		nb,maxbottles,packsize.begin(),
		rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                &nextbflag,0);
      } else {
	nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                &nextbflag,0);
      }

      iinvout(iobsbig) = inv;
      ibout(iobsbig) = b+1;

      iobsbig++;
      inv = invprime;
      b = bprime;

    }

    initinv[n] = inv;
    initb[n] = b;

  }
  
  Rcpp::NumericVector invout(nobs);
  Rcpp::IntegerVector bout(nobs);

  for(int hh=hhinds(0)-1;hh<hhinds(1);hh++) {

    inv = initinv[hh];
    b = initb[hh];

    double rcnew = hhcrate(hh);
    
    int obsend;
    if(hh == nhhs-1) {
      obsend = nobs;
    } else {
      obsend = obshhinds(hh+1)-1;
    }

    for(int n=obshhinds(hh)-1;n<obsend;n++) {

      int units = (int)tunits(n);

      if(units < 0) {
	nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                &nextbflag,0);
      } else {
	nextinv(tunits(n),brsize(brindex(n)-1+nbrand)-1,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                &nextbflag,0);
      }

      invout(n) = inv;
      bout(n) = b+1;

      inv = invprime;
      b = bprime;

    }

  }

  Rcpp::List ret;

  ret["initinv"] = iinvout;
  ret["initb"] = ibout;
  
  ret["inv"] = invout;
  ret["b"] = bout;

  return ret;
  
  END_RCPP

}



// ols wrapper for debugging

// [[Rcpp::export]]
Rcpp::List olswrap(Rcpp::NumericMatrix Xin, Rcpp::NumericVector yin) {

  BEGIN_RCPP

  arma::mat X = Rcpp::as<arma::mat>(Xin);
  arma::colvec y = Rcpp::as<arma::colvec>(yin);

  int ncoef = Xin.ncol();
  int nobs = Xin.nrow();

  arma::colvec coef(ncoef);
  arma::colvec resid(nobs);

  double cnum=0;

  ols(X,y,coef,resid,&cnum,0);

  Rcpp::NumericVector outcoef(ncoef);
  Rcpp::NumericVector outresid(nobs);

  coef.print();

  std::copy(coef.begin(),coef.end(),outcoef.begin());
  std::copy(resid.begin(),resid.end(),outresid.begin());

  Rcpp::List res;
  res["coef"] = outcoef;
  res["resid"] = outresid;

  return res;

  END_RCPP

}


// fitiv wrapper - for checking code

// [[Rcpp::export]]
Rcpp::List fitivwrap(const Rcpp::NumericMatrix &co,  const Rcpp::NumericVector &panid,
                     const Rcpp::NumericMatrix &pricemat,
                     const Rcpp::IntegerVector &hhinds, const Rcpp::IntegerVector &brsize,
                     const Rcpp::IntegerVector &obshhinds, const int nco, const int nobs,
                     const int nbrand, const int nsize, const int nhhs, const int sptype,
                     const int ncutiv) {


  BEGIN_RCPP
    
    Rcpp::NumericMatrix riv(nobs,nsize);

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;

  Rcpp::NumericVector rcutmat(cutblocksize*nhhs);

  int nrowiv=0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivcoefblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  Rcpp::NumericVector rivcoef(ivcoefblocksize*nhhs);
  
  int ivvariblocksize = nsize*nsize;
  Rcpp::NumericVector rivvari(ivvariblocksize*nhhs);
  
  Rcpp::NumericVector rdetiv(nhhs);

  fitiv(co,panid,pricemat,hhinds,brsize,obshhinds,nco,nobs,
        nbrand,nsize,nhhs,sptype,ncutiv,riv,
        rcutmat,rivcoef,rivvari,rdetiv,0);

  Rcpp::List res;

  res["iv"] = riv;
  res["cutmat"] = rcutmat;
  res["ivcoef"] = rivcoef;
  res["ivvari"] = rivvari;
  res["detiv"] = rdetiv;

  return res;
    
  END_RCPP
    
    }


// code to do counterfactuals at MCMC draws
// cfdata should be an alternative price matrix

// parallel worker for initial choices

struct initchoicepar : public RcPar::Worker
{

  //input variables
  const RcPar::RMatrix<double> poldtf;
  const RcPar::RVector<int> iobsind;
  const RcPar::RVector<int> iobsindbig;
  const RcPar::RMatrix<double> rgrid;
  const RcPar::RVector<int> brsize;
  const RcPar::RVector<double> rcdraws;
  const RcPar::RVector<double> initx;
  const RcPar::RVector<int> stvisitinit;
  const RcPar::RMatrix<double> cfinit;
  const RcPar::RVector<double> cosave;
  const RcPar::RVector<int> indexes;
  const RcPar::RMatrix<double> bwmatrix;
  const RcPar::RVector<double> rivcoefout;
  const RcPar::RVector<double> rcutmatout;
  const RcPar::RVector<double> rivvariout;
  const RcPar::RVector<double> rdetivout;
  const RcPar::RVector<double> impdist;
  const RcPar::RVector<double> packsize;
  const RcPar::RVector<int> bstates;
  const RcPar::RVector<int> revarray;
  const RcPar::RVector<int> badmissable;
  const RcPar::RVector<int> bindex;
  const RcPar::RVector<double> pgrid;
  const RcPar::RVector<double> vfsave;
  const RcPar::RVector<double> rgridsave;
  const RcPar::RVector<double> gvari;
  const RcPar::RVector<double> gmean;
  const RcPar::RMatrix<double> logiterror;
  const RcPar::RVector<int> ivdrop;

  //const int nsave;
  //const int nbsmall;
  const double * istates;
  const int rep;
  const int nistates;
  const int maxpind;
  const int nd;
  const int linitx;
  const int vfinterpmethod;
  const int ivblocksize;
  const int cutblocksize;
  const int gridmethod;
  const llinputs llinfo;

  //output variables
  RcPar::RVector<int> ichoicej;
  RcPar::RVector<int> ichoices;
  RcPar::RVector<double> initinv;
  RcPar::RVector<int> initb;
  RcPar::RVector<double> risave1;
  RcPar::RVector<int> ibsave;
  
  // initialize with source and destination
  initchoicepar(const Rcpp::NumericMatrix poldtf, const Rcpp::IntegerVector iobsind,
		const Rcpp::IntegerVector iobsindbig, const Rcpp::NumericMatrix rgrid,
		const Rcpp::IntegerVector brsize, const Rcpp::NumericVector rcdraws,
		const Rcpp::NumericVector initx, const Rcpp::IntegerVector stvisitinit,
		const Rcpp::NumericMatrix cfinit, const Rcpp::NumericVector cosave,
		const Rcpp::IntegerVector indexes, const Rcpp::NumericMatrix bwmatrix,
		const Rcpp::NumericVector rivcoefout, const Rcpp::NumericVector rcutmatout,
		const Rcpp::NumericVector rivvariout, const Rcpp::NumericVector rdetivout,
		const Rcpp::NumericVector impdist, const Rcpp::NumericVector packsize,
		const Rcpp::IntegerVector bstates, const Rcpp::IntegerVector revarray,
		const Rcpp::IntegerVector badmissable, const Rcpp::IntegerVector bindex,
		const Rcpp::NumericVector pgrid, const Rcpp::NumericVector vfsave,
		const Rcpp::NumericVector rgridsave, const Rcpp::NumericVector gvari,
		const Rcpp::NumericVector gmean, const Rcpp::NumericMatrix logiterror,
		const Rcpp::IntegerVector ivdrop,
	        const double * istates, const int rep,
		const int nistates, const int maxpind, const int nd, const int linitx,
		const int vfinterpmethod, const int ivblocksize, const int cutblocksize,
		const int gridmethod, const llinputs llinfo,
		Rcpp::IntegerVector ichoicej, Rcpp::IntegerVector ichoices,
		Rcpp::NumericVector initinv, Rcpp::IntegerVector initb,
		Rcpp::NumericVector risave1, Rcpp::IntegerVector ibsave) :
    poldtf(poldtf), iobsind(iobsind), iobsindbig(iobsindbig), rgrid(rgrid), brsize(brsize),
    rcdraws(rcdraws), initx(initx), stvisitinit(stvisitinit), cfinit(cfinit), cosave(cosave),
    indexes(indexes), bwmatrix(bwmatrix), rivcoefout(rivcoefout), rcutmatout(rcutmatout),
    rivvariout(rivvariout), rdetivout(rdetivout), impdist(impdist), packsize(packsize),
    bstates(bstates), revarray(revarray), badmissable(badmissable), bindex(bindex),
    pgrid(pgrid), vfsave(vfsave), rgridsave(rgridsave), gvari(gvari), gmean(gmean),
    logiterror(logiterror), ivdrop(ivdrop), istates(istates), rep(rep),
    nistates(nistates), maxpind(maxpind), nd(nd), linitx(linitx), vfinterpmethod(vfinterpmethod),
    ivblocksize(ivblocksize), cutblocksize(cutblocksize), gridmethod(gridmethod), llinfo(llinfo),
    ichoicej(ichoicej), ichoices(ichoices), initinv(initinv), initb(initb), risave1(risave1),
    ibsave(ibsave) {}

  void operator()(std::size_t begin, std::size_t end) {
    
    int begin1 = (int)begin;
    int end1 = (int)end;

    int nsave = llinfo.nsave;
    int nbsmall = llinfo.nbsmall;
    int nbrand = llinfo.nbrand;
    int nrgrid = llinfo.nrgrid;
    int nco = llinfo.nco;
    int nsize = llinfo.nsize;
    int ninitt = llinfo.ninitt;
    double invub = llinfo.maxinv;
    double maxinv = llinfo.maxinv;
    int capj = llinfo.capj;
    int myopic = llinfo.myopic;
    int nhhs = llinfo.nhhs;
    int cmodel = llinfo.cmodel;
    int ncut = llinfo.ncut;
    int inttype = llinfo.inttype;
    int sptype = llinfo.sptype;
    int invmodel = llinfo.invmodel;
    int maxbottles = llinfo.maxbottles;
    int first20 = llinfo.first20;
    double dg = llinfo.dg;
    int nb = llinfo.nb;
    int hinvbound = llinfo.hinvbound;
    double inv = 0, invprime = 0;
    int b=0, bprime=0;
    //memset(initb,0,nhhs*sizeof(int));
    int coorder[nsave];
    memset(coorder,0,nsave*sizeof(int));
    int iobs = iobsind[begin1];
    int iobsbig = iobsindbig[begin1];
    double ll = 0;

    double ut[(1+capj)*nsize];
    int drop[(1+capj)*nsize];
    double vfhat[(1+capj)*nsize];
    double cerror[(1+capj)*nsize];
    int bchoice[(1+capj)*nsize];
    double ichoice[(1+capj)*nsize];

    double ivout[nsize];
    double ivgrid[(1+capj)*nsize];
    double ivtprobs[nrgrid*maxpind];
    double rgridhh[(nd-1)*llinfo.nrgrid];
    double vftemp[nbsmall*nistates*maxpind];  
    int vffilled[nbsmall*nistates*maxpind];
    for(int i=0;i<nbsmall*nistates*maxpind;i++) {
      vffilled[i] = -1;
    }
    int ivtpfilled[maxpind];
    for(int i=0;i<maxpind;i++) {
      ivtpfilled[i] = -1;
    }

    int nrep = 0, nreptilde=0;
    if(rep <= nsave) {
      nrep = rep-1;
    } else {
      nrep = nsave;
    }
    int ntilde = llinfo.ntilde;
    if(nrep < ntilde) {
      nreptilde = rep-1;
    } else {
      nreptilde = ntilde;
    }

    int prcheck = 0;

    //int capj = llinfo.capj;
    int lomega = llinfo.lomega;
    int hmodel = llinfo.hmodel;
    double clb, cub, alpha, gamma, nu, beta, ccost;
    double xi[nbrand];
    double omega[lomega];

    double kweights[nsave];
    int vfind=0;
      
    for(int n=begin1;n<end1;n++) {
      //Rcpp::Rcout << "hh: " << n << std::endl;
      inv = 0;
      b = 0;
      for(int j=0;j<lomega;j++) {
	omega[j] = poldtf(nbrand+9+j,n);
      }
      clb = poldtf(nbrand,n);
      cub = poldtf(nbrand+1,n);

      alpha = poldtf(nbrand+2,n);
      gamma = poldtf(nbrand+3,n);
      nu = poldtf(nbrand+4,n);
      beta = poldtf(nbrand+5,n);
      ccost = poldtf(nbrand+7,n);

      vfind = n*nrgrid*nbsmall*nistates*nsave;
      
      if(gridmethod == 2) {
	for(int p=0;p<nrgrid;p++) {
	  indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);
	  
	  computeiv(poldtf.begin()+nco*n,rgrid.begin()+nbrand*p,brsize.begin(),
		    nsize,nbrand,rgridhh+(nd-1)*p);

	}
      }
      
      for(int t=0;t<ninitt;t++) {

	//if(rep > 0 && n == 191) {Rcpp::Rcout << "tperiod " << t << std::endl;}

	double chosenq = 0;
	double rcnew;
	//if(usecdraws) {
	indcheck(n+nhhs*t,0,llinfo.lenicdraws,"rcnew offset bug (initpar)");
	rcnew = rcdraws[n+nhhs*t];
	//} else {
	//  rcnew = Rcpp::runif(1)[0];
	//}
	if(rcnew < 0.5) {
	  rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
	} else {
	  rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
	}
	indcheck(t+ninitt*n,0,linitx,"initx: n: %d; t %d",n,t);
	chosenq = initx[t+ninitt*n];  // since this vector is ordered by id, then week, this should be the right indexing
	// chosenq will get overwritten if we're simulating initial choices
    
	invub = hinvbound ? omega[0] : maxinv;
      
	double ivinit[nsize];

	indcheck(iobsbig,0,llinfo.nobsinit,"iobsbig offset bug (initpar)");
	if(stvisitinit[iobsbig]) {
	  iobs++;
	}

	if(stvisitinit[iobsbig]) {

	  // If we are resimulating choices we should create new inclusive values

	  int nprint = -1; //282;
	  int tprint = -1; //1;
	
	  int nbravail = 0;
          int drop1[nbrand];
	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << "****** n: " << n << ", t: " << t << std::endl;
	  }
	
	  for(int size=0;size<nsize;size++) {
	    double umax = -10000000;
	    for(int j=0;j<nbrand;j++) {
	      double price = cfinit(iobs,j);
	      if(price >= 999 || brsize[j+nbrand] != size+1) {
		drop1[j] = 1;
	      } else {
		nbravail++;
		ut[j] = poldtf(j,n) + alpha*price;
		drop1[j] = 0;
		if(ut[j] > umax)
		  umax = ut[j];
		if(n == nprint && t <= tprint) {
		  Rcpp::Rcout << " size: " << size << " ut[" << j << "]: " << ut[j] << "price: " << price << ". ";
		}
	      }
	    }
	    if(n == nprint && t <= tprint) {
	      Rcpp::Rcout << std::endl;
	    }
	    double s = 0;
	    for(int j=0;j<nbrand;j++) {
	      if(!drop1[j]) {
		s += exp(ut[j]-umax);
	      }
	    }
	    s=log(s)+umax;
	    ivinit[size] = s;
	  
	    for(int j=0;j<1+capj;j++) {
	      ivgrid[j+size*(1+capj)] = ivinit[size]*((double)j)/((double)capj);
	    }
	  }

	  //Note that the code below doesn't make use of price states.  It's not worth the trouble
	  //since this part only gets run once
	  //if(rep > 0) {Rcpp::Rcout << "aaa " << t << std::endl;}

	  if(!myopic) {
	    if(vfinterpmethod == 1) {
	      getclosestind(poldtf.begin() + nco*n, cosave.begin() + nsave*nco*n, indexes.begin(), nco, nrep, nsave, rep, coorder);
	    } else {
	      if(!(hmodel && n > begin1)) {
		getkweights(poldtf.begin() + nco*n, cosave.begin() + nsave*nco*n,
			    indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
			    rep, coorder, kweights,0);
	      }
	    }
	  }
	  //for(int j=0;j<nsave;j++) {
	  //  Rcpp::Rcout << "coorder[" << j << "]: " << coorder[j] << std::endl;
	  //}
	  if(!myopic && beta > 0) {
	    ivpred(ivinit, rivcoefout.begin() + ivblocksize*n, rcutmatout.begin() + cutblocksize*n, nd, ncut*(ivdrop[n]==0), sptype, inttype, capj, ivout);
	    if(gridmethod == 2) {
	      for(int p=0;p<nrgrid;p++) {
		ivtprobs[p] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,rivvariout.begin()+(nd-1)*(nd-1)*n,
				     rdetivout[n],impdist[p]);
		if(n==nprint && t <= tprint && p < 10) {
		
		  Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p] << ", " << impdist[p] <<
		    ", " << rdetivout[n] << std::endl;
		  for(int ii=0;ii<nd-1;ii++) {
		    Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
		  }
		  Rcpp::Rcout << std::endl;
		}
	      }
	    }

	    int printflag = 0;

	    //if(n == nprint && t <= tprint) {
	    //  Rcpp::Rcout << "vfin: " << vfin[0] << ", " << vfin[1] << ", " << vfin[2] << std::endl;
	    //}

            int nextbflag = 0;
	    
	    for(int size=0;size<nsize;size++) {
	      for(int j=0;j<capj+1;j++) {
		if( (j == 0 && size == 0) || j > 0 ) {
		  int bprime1 = 0;
		  if(invmodel == 1) {
		    invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		    invprime = invprime > invub ? invub : invprime;
		  } else {
		    nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			    rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		    bprime1 = bprime;
		    bprime = bindex[bprime]-1;
		    indcheck(bprime,0,nbsmall,"bprime (init sim) error");
		  }

		  int iind = itoiind(invprime,istates,nistates,first20);

		  // uncomment if not using price states
		  int vfindx1 = iind + nistates*bprime;
		  int ok = iind < nistates; //check to make sure the state above this one is admissable (I think the second part should always be true)
		  if(invmodel == 2) {
		    ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
		  }
		  indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vf index error cf init 1");
		  if(ok) {
		    if(vffilled[vfindx1] < 0) {
		      vffilled[vfindx1] = 1;
		      if(gridmethod==2) {
			vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						      nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*n,
						      gvari.begin(), rdetivout[n], dg, gmean.begin(),
						      nd-1, istates, nistates, first20, bprime, nbsmall, n == nprint && t <= tprint,0,vfinterpmethod);

		      } else {
			vftemp[vfindx1] = vfavgijchh(ivout, iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						     nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*n,
						     gvari.begin(), rdetivout[n], dg, gmean.begin(),
						     nd-1, istates, nistates, first20, bprime, nbsmall, printflag);
		      }
		    
		      if(n == nprint && t <= tprint) {
			Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		      }

		    }
		  }

		  vfindx1 = iind-1 + nistates*bprime;
		  indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vf index error cf init 2");
		  if(vffilled[vfindx1] < 0) {
		    vffilled[vfindx1] = 1;
		    if(gridmethod == 2) {
		      vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						    nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*n,
						    gvari.begin(), rdetivout[n], dg, gmean.begin(),
						    nd-1, istates, nistates, first20, bprime, nbsmall, n == nprint && t <= tprint, 0,vfinterpmethod);
		    } else {
		      vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						   nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*n,
						   gvari.begin(), rdetivout[n], dg, gmean.begin(),
						   nd-1, istates, nistates, first20, bprime, nbsmall, 0);
		    }

		    if(n == nprint && t <= tprint) {
		      Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		    }

		  }

		}
	      }
	    }
	  }
	  //throw std::range_error("stopping");
	  int first = 1;

	  int vfoffset = 0;
	  double maxutil = -1000000;
          int nextbflag = 0;
          
	  for(int size=0;size<nsize;size++) {
	    for(int j=0;j<capj+1;j++) {
	      int indx = j+(1+capj)*size;
	      drop[indx] = 1;
	      if( (j == 0 && size == 0) || j > 0 ) {
		drop[indx] = 0;
		if(invmodel == 1) {
		  invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		  invprime = invprime > invub ? invub : invprime;
		  ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize[size], rcnew, gamma, nu, ccost, omega, cmodel, 0);
		} else {
		  nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		  bchoice[indx] = bprime;
		  indcheck(bprime,0,nb,"bprime (init sim) error");
		  ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				   nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
		  bprime = bindex[bprime]-1;		
		  ichoice[indx] = invprime;
		}
              
		if(!myopic && beta > 0) {
		  vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20,vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
		} else {
		  vfhat[indx] = 0.0;
		}
		if(n == nprint && t <= tprint) {
		  Rcpp::Rcout << "ut[" << indx << "]: " << ut[indx] << ", vf: " << vfhat[indx] << ", ivgrid: " << ivgrid[indx] <<std::endl;
		}
		maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);

	      }
	    }
	  }

	  double s1 = 0.0;

	  int schoice=0;
	  int jchoice=0;

	  double choiceumax = -1000000;

	  for(int size=0;size<nsize;size++) {
	    for(int j=0;j<capj+1;j++) {
	      int indx = j+(1+capj)*size;
	      if(!drop[indx]) {
		ut[indx] -= maxutil;
		indcheck(n+nhhs*(t+ninitt*(size+nsize*j)),0,nhhs*ninitt*(1+capj)*nsize,"logiterror offset bug (initpar)");
		cerror[indx] = logiterror(n+nhhs*t,size+nsize*j);
		s1 += exp(ut[indx]+vfhat[indx]);
		if(ut[indx]+vfhat[indx]+cerror[indx] > choiceumax) {
		  choiceumax = ut[indx]+vfhat[indx]+cerror[indx];
		  if(invmodel==2) {
		    bprime = bchoice[indx];
		    invprime = ichoice[indx];
		  }
		  if(j > 0) {
		    schoice = size+1;
		    jchoice = j;
		  }
		}
	      }
	    }
	  }

	  if(n == nprint && t <= tprint) {
	    Rcpp::Rcout << "DEBUG: b: " << b << ", inv : " << inv <<
	      ", bprime: " << bprime << ", invprime; " << invprime <<
	      ", jchoice: " << jchoice << ", schoice: " << schoice << std::endl;
	  
	    int tempb; double tempi;
            if(invmodel == 2) {
              nextinv(jchoice,schoice-1,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
                      rcnew,revarray.begin(),nsize,badmissable.begin(),&tempb,&tempi,&nextbflag,0);
            }
	    //Rcpp::Rcout << "TEST: " << tempb << ", " << tempi << std::endl;
	    for(int size=0;size<nsize;size++) {
	      for(int j=0;j<capj+1;j++) {
		int indx = j+(1+capj)*size;
		if(!drop[indx]) {
		  Rcpp::Rcout << "choice: " << indx << " (" << j << ", " << size << "): utility: " << ut[indx]
			      << " vf: " << vfhat[indx] << " error: " << cerror[indx]
			      << " bchoice: " << bchoice[indx] << " ichoice: " << ichoice[indx] << std::endl;		
		}
	      }
	    }
	    //if(t==tprint) {throw std::range_error("stopping");}
	  }

	  double obspack = 0;

	  if(nbravail > 0 && schoice > 0) {
	    chosenq = (double)(packsize[schoice-1]*jchoice);
	    //ichoicej(t+ninitt*(n+nhhs*rep)) = jchoice;
	    //ichoices(t+ninitt*(n+nhhs*rep)) = schoice;
	    indcheck(t+ninitt*n,0,linitx,"initchoice offset bug (initpar)");
	    ichoicej[t+ninitt*n] = jchoice;
	    ichoices[t+ninitt*n] = schoice;
	  } else {
	    chosenq = 0;
	    //ichoicej(t+ninitt*(n+nhhs*rep)) = 0;
	    //ichoices(t+ninitt*(n+nhhs*rep)) = 0;
	    indcheck(t+ninitt*n,0,linitx,"initchoice offset bug (initpar)");
	    ichoicej[t+ninitt*n] = 0;
	    ichoices[t+ninitt*n] = 0;
            if(invmodel == 2) {
              bprime = bchoice[0];
              invprime = ichoice[0];
            }
	  }

	  // reset vffilled
	  memset(vffilled,-1,nb*nistates*maxpind*sizeof(int));

	} else if(invmodel == 2) {
          int nextbflag = 0;
	  nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
	}

	//ichoiceout(t+ninitt*n) = chosenq;
	indcheck(n+nhhs*t,0,llinfo.nsim*nhhs*llinfo.ninitt,"risave1 offset bug (initpar)");
	risave1[n+nhhs*t] = inv;
	//ibsave(t+ninitt*(n+nhhs*rep)) = b;
	indcheck(t+ninitt*n,0,linitx,"initb offset bug (initpar)");
	ibsave[t+ninitt*n] = b;
	if(invmodel == 1) {
	  inv = dmax(inv + chosenq - rcnew,0);
	  inv = inv > invub ? invub : inv;
	} else {
	  inv = invprime;
	  b = bprime;
	}
      
	iobsbig++;

      
      }
      initinv[n] = inv;
      if(invmodel==2) {
	initb[n] = b;
      }

      //if(rep>0) {Rcpp::Rcout << "done init hh " << n << std::endl;}

    }

    
  }

  
};

// parallel worker for estimation sample choices

struct datachoicepar : public RcPar::Worker
{

  //input variables
  const RcPar::RMatrix<double> poldtf;
  const RcPar::RVector<int> obshhindsmerge;
  const RcPar::RVector<int> obshhinds;
  const RcPar::RMatrix<double> rgrid;
  const RcPar::RVector<int> brsize;
  const RcPar::RVector<double> cdraws;
  const RcPar::RVector<double> initx;
  const RcPar::RVector<int> stvisitinit;
  const RcPar::RMatrix<double> cfdata;
  const RcPar::RVector<double> cosave;
  const RcPar::RVector<int> indexes;
  const RcPar::RMatrix<double> bwmatrix;
  const RcPar::RVector<double> rivcoefout;
  const RcPar::RVector<double> rcutmatout;
  const RcPar::RVector<double> rivvariout;
  const RcPar::RVector<double> rdetivout;
  const RcPar::RVector<double> impdist;
  const RcPar::RVector<double> packsize;
  const RcPar::RVector<int> bstates;
  const RcPar::RVector<int> revarray;
  const RcPar::RVector<int> badmissable;
  const RcPar::RVector<int> bindex;
  const RcPar::RVector<double> pgrid;
  const RcPar::RVector<double> vfsave;
  const RcPar::RVector<double> rgridsave;
  const RcPar::RVector<double> gvari;
  const RcPar::RVector<double> gmean;
  const RcPar::RMatrix<double> logiterror;
  const RcPar::RVector<int> ivdrop2;
  // these below are new to this routine
  const RcPar::RVector<double> tunits;
  const RcPar::RMatrix<double> ivbig;
  const RcPar::RVector<double> initinv;
  const RcPar::RVector<int> initb;
  const RcPar::RVector<int> pinds;
  const RcPar::RMatrix<double> logiterrorb;
  const RcPar::RVector<double> panid;
  const RcPar::RVector<int> hhinds;

  //const int nsave;
  //const int nbsmall;
  const double * istates;
  const int rep;
  const int nistates;
  const int maxpind;
  const int nd;
  const int linitx;
  const int vfinterpmethod;
  const int ivblocksize;
  const int cutblocksize;
  const int gridmethod;
  const llinputs llinfo;
  const int burnin;
  const int cfkeep;
  const int keep;
  const int nobsdata;
  const int leniv;

  //output variables
  RcPar::RVector<int> simchoice;
  RcPar::RMatrix<double> vfout;
  RcPar::RMatrix<double> utout;
  //RcPar::RVector<int> savechoice;
  RcPar::RVector<double> risave2;
  
  // initialize with source and destination
  datachoicepar(const Rcpp::NumericMatrix poldtf, const Rcpp::IntegerVector obshhindsmerge,
		const Rcpp::IntegerVector obshhinds, const Rcpp::NumericMatrix rgrid,
		const Rcpp::IntegerVector brsize, const Rcpp::NumericVector cdraws,
		const Rcpp::NumericVector initx, const Rcpp::IntegerVector stvisitinit,
		const Rcpp::NumericMatrix cfdata, const Rcpp::NumericVector cosave,
		const Rcpp::IntegerVector indexes, const Rcpp::NumericMatrix bwmatrix,
		const Rcpp::NumericVector rivcoefout, const Rcpp::NumericVector rcutmatout,
		const Rcpp::NumericVector rivvariout, const Rcpp::NumericVector rdetivout,
		const Rcpp::NumericVector impdist, const Rcpp::NumericVector packsize,
		const Rcpp::IntegerVector bstates, const Rcpp::IntegerVector revarray,
		const Rcpp::IntegerVector badmissable, const Rcpp::IntegerVector bindex,
		const Rcpp::NumericVector pgrid, const Rcpp::NumericVector vfsave,
		const Rcpp::NumericVector rgridsave, const Rcpp::NumericVector gvari,
		const Rcpp::NumericVector gmean, const Rcpp::NumericMatrix logiterror,
		const Rcpp::IntegerVector ivdrop2, const Rcpp::NumericVector tunits,
		const Rcpp::NumericMatrix ivbig, const Rcpp::NumericVector initinv,
		const Rcpp::IntegerVector initb, const Rcpp::IntegerVector pinds,
		const Rcpp::NumericMatrix logiterrorb, Rcpp::NumericVector panid,
		const Rcpp::IntegerVector hhinds,
	        const double * istates, const int rep, 
		const int nistates, const int maxpind, const int nd, const int linitx,
		const int vfinterpmethod, const int ivblocksize, const int cutblocksize,
		const int gridmethod, const llinputs llinfo, const int burnin,
		const int cfkeep, const int keep, const int nobsdata, const int leniv,
		Rcpp::IntegerVector simchoice, Rcpp::NumericMatrix vfout,
		Rcpp::NumericMatrix utout, Rcpp::NumericVector risave2) :
    poldtf(poldtf), obshhindsmerge(obshhindsmerge), obshhinds(obshhinds), rgrid(rgrid), brsize(brsize),
    cdraws(cdraws), initx(initx), stvisitinit(stvisitinit), cfdata(cfdata), cosave(cosave),
    indexes(indexes), bwmatrix(bwmatrix), rivcoefout(rivcoefout), rcutmatout(rcutmatout),
    rivvariout(rivvariout), rdetivout(rdetivout), impdist(impdist), packsize(packsize),
    bstates(bstates), revarray(revarray), badmissable(badmissable), bindex(bindex),
    pgrid(pgrid), vfsave(vfsave), rgridsave(rgridsave), gvari(gvari), gmean(gmean),
    logiterror(logiterror), ivdrop2(ivdrop2), tunits(tunits), ivbig(ivbig), initinv(initinv),
    initb(initb), pinds(pinds), logiterrorb(logiterrorb), panid(panid),
    hhinds(hhinds), istates(istates), rep(rep),
    nistates(nistates), maxpind(maxpind), nd(nd), linitx(linitx), vfinterpmethod(vfinterpmethod),
    ivblocksize(ivblocksize), cutblocksize(cutblocksize), gridmethod(gridmethod), llinfo(llinfo),
    burnin(burnin), cfkeep(cfkeep), keep(keep), nobsdata(nobsdata), leniv(leniv),
    simchoice(simchoice), vfout(vfout), utout(utout), risave2(risave2) {}

  void operator()(std::size_t begin, std::size_t end) {
    
    int begin1 = (int)begin;
    int end1 = (int)end;

    int nsave = llinfo.nsave;
    int nbsmall = llinfo.nbsmall;
    int nbrand = llinfo.nbrand;
    int nrgrid = llinfo.nrgrid;
    int nco = llinfo.nco;
    int nsize = llinfo.nsize;
    int ninitt = llinfo.ninitt;
    double invub = llinfo.maxinv;
    double maxinv = llinfo.maxinv;
    int capj = llinfo.capj;
    int myopic = llinfo.myopic;
    int nhhs = llinfo.nhhs;
    int cmodel = llinfo.cmodel;
    int ncut = llinfo.ncut;
    int inttype = llinfo.inttype;
    int sptype = llinfo.sptype;
    int invmodel = llinfo.invmodel;
    int maxbottles = llinfo.maxbottles;
    int first20 = llinfo.first20;
    double dg = llinfo.dg;
    int nb = llinfo.nb;
    int hinvbound = llinfo.hinvbound;
    int debug = llinfo.debug;
    
    double inv = 0, invprime = 0;
    int b=0, bprime=0;
    //memset(initb,0,nhhs*sizeof(int));
    int coorder[nsave];
    memset(coorder,0,nsave*sizeof(int));
    //int iobs = iobsind[begin1];
    //int iobsbig = iobsindbig[begin1];
    double ll = 0;
    //double rllrun = 0;

    double ut[(1+capj)*nsize];
    int drop[(1+capj)*nsize];
    double vfhat[(1+capj)*nsize];
    double cerror[(1+capj)*nsize];
    int bchoice[(1+capj)*nsize];
    double ichoice[(1+capj)*nsize];

    double ivout[nsize];
    double ivgrid[(1+capj)*nsize];
    double ivtprobs[nrgrid*maxpind];
    double rgridhh[(nd-1)*llinfo.nrgrid];
    double vftemp[nbsmall*nistates*maxpind];  
    int vffilled[nbsmall*nistates*maxpind];
    for(int i=0;i<nbsmall*nistates*maxpind;i++) {
      vffilled[i] = -1;
    }
    int ivtpfilled[maxpind];
    for(int i=0;i<maxpind;i++) {
      ivtpfilled[i] = -1;
    }

    int nrep = 0, nreptilde=0;
    if(rep <= nsave) {
      nrep = rep-1;
    } else {
      nrep = nsave;
    }
    int ntilde = llinfo.ntilde;
    if(nrep < ntilde) {
      nreptilde = rep-1;
    } else {
      nreptilde = ntilde;
    }

    int prcheck = 0;

    //int capj = llinfo.capj;
    int lomega = llinfo.lomega;
    int hmodel = llinfo.hmodel;
    double clb, cub, alpha, gamma, nu, beta, ccost;
    double xi[nbrand];
    double omega[lomega];

    double kweights[nsave];
    int vfind=0;

    
    
    int hh = hhinds[0]-2;
    double llrun = 0;
    double llhh[nhhs];
    double llhhbr[nhhs];

    //double invub = maxinv;


    // importance probs - for checking
    //Rcpp::NumericVector tprobout(obsend*nrgrid);

    int n2 = obshhindsmerge[hhinds[0]-1]-1;

    int cfdataind = -1;

    //prcheck = rep == 32;

    for(int hh=begin1;hh<end1;hh++) {
      cfdataind = obshhindsmerge[hh]-2;  //there's an increment below so need this
      if(prcheck) {Rcpp::Rcout << "hh: " << hh << std::endl;}

      if(!hmodel) {
        memset(vffilled,-1,nbsmall*nistates*maxpind*sizeof(int));
        memset(ivtpfilled,-1,maxpind*sizeof(int));
      }

      indcheck(hh,0,nhhs,"initinv bug (datachoicepar). hh %d",hh);
      inv = initinv[hh];
      b = initb[hh];
      llrun = 0;
      llhh[hh] = 0;
      llhhbr[hh] = 0;
      vfind = hh*nrgrid*nbsmall*nistates*nsave;

      for(int j=0;j<nbrand;j++) {
        xi[j] = poldtf(j,hh);
      }

      clb = poldtf(nbrand,hh);
      cub = poldtf(nbrand+1,hh);

      alpha = poldtf(nbrand+2,hh);
      gamma = poldtf(nbrand+3,hh);
      nu = poldtf(nbrand+4,hh);
      beta = poldtf(nbrand+5,hh);
      ccost = poldtf(nbrand+7,hh);

      for(int j=0;j<lomega;j++) {
	omega[j] = poldtf(nbrand+9+j,hh);
      }

      invub = hinvbound ? omega[0] : maxinv;

      if(!myopic) {
	if(vfinterpmethod == 1) {
	  getclosestind(poldtf.begin() + nco*hh, cosave.begin() + nsave*nco*hh, indexes.begin(), nco, nrep, nsave, rep, coorder);
	} else {
	  if(!(hmodel && hh > hhinds[0]-1)) {
	    getkweights(poldtf.begin() + nco*hh, cosave.begin() + nsave*nco*hh,
			indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
			rep, coorder, kweights,0);
	  }
	}
      }

      int obsend;
      if(hh == nhhs-1) {
	obsend = llinfo.nobs;
      } else {
	obsend = obshhinds[hh+1]-1;
      }

      if(gridmethod == 2) {
	for(int p=0;p<nrgrid;p++) {
	  indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);

	  computeiv(poldtf.begin()+nco*hh,rgrid.begin()+nbrand*p,brsize.begin(),
		    nsize,nbrand,rgridhh+(nd-1)*p);

	}
      }
      //Rcpp::Rcout << "abc" << std::endl;
      for(int n=obshhinds[hh]-1;n<obsend;n++) {
	if(prcheck) {Rcpp::Rcout << "1" << std::endl;}
	indcheck(n,0,llinfo.nobs,"tunits bug (datachoicepar). n: %d; hh %d",n,hh);
	int units = (int)tunits[n];

	double rcnew;
	indcheck(n,0,llinfo.lencdraws,"cdraws. n: %d",n);
	//if(usecdraws) {
	rcnew = cdraws[n];
	//} else {
	//  rcnew = R::runif(0,1);
	//}

	if(rcnew < 0.5) {
	  rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
	} else {
	  rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
	}

        int nextbflag = 0;

	if(units < 0) {
	  // household didn't go to the store - there's nothing to model here.
	  risave2[n] = inv;
	  if(invmodel==1) {
	    inv = dmax(inv - rcnew,0);
	  } else {
	    nextinv(0,0,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		    rcnew,revarray.begin(),nsize,badmissable.begin(),&b,&inv,
                    &nextbflag,0);
	  }
	  //if(debug) {
	  //  rllrun(n) = 1;
	  //}
	  indcheck(n+nobsdata*2,0,llinfo.nobs*3,"simchoice bug (datachoicepar). n: %d; hh %d",n,hh);
	  simchoice[n+nobsdata*0] = 0;
	  simchoice[n+nobsdata*1] = -1;
	  simchoice[n+nobsdata*2] = -1;
	} else {
	  cfdataind++;
	  indcheck(n,0,leniv,"iv. n: %d",n);
	  //Rcpp::Rcout << "1" << std::endl;
	  double maxutil = -1000000;

	  // compute avg predicted inclusive value, needed for importance weighting
	  double rgrid1[nsize];
	  //indcheck(n2,0,nrowivbig,"ivbig: n: %d",n);
	  for(int size=0;size<nsize;size++) {
	    rgrid1[size] = ivbig(n,size);
	  }

	  //int printobs = n <= 51141 ? 51141 : 52607;
	  int printobs = -1; //38943; //n <= 158 ? 158 : 17378;;
	  if(prcheck) {Rcpp::Rcout << "aaa" << std::endl;}
	  if(!myopic && beta > 0) {
	    ivpred(rgrid1, rivcoefout.begin() + ivblocksize*hh, rcutmatout.begin() + cutblocksize*hh, nd, ncut*(ivdrop2[hh]==0), sptype, inttype, capj, ivout);
	    if(gridmethod == 2) {
	      if(n==printobs) {
		Rcpp::Rcout << "pinds[" << n << "]: " << pinds[n] << "; ivptfilled[" << n << "]: " << ivtpfilled[pinds[n]-1] << std::endl;
		for(int ii=0;ii<nsize;ii++) {
		  Rcpp::Rcout << "rgrid1[" << ii << "]: " << rgrid1[ii] << ", ";
		}
		Rcpp::Rcout << std::endl;
		for(int ii=0;ii<nd-1;ii++) {
		  Rcpp::Rcout << "ivout[" << ii << "]: " << ivout[ii] << ", ";
		}
		Rcpp::Rcout << std::endl;
		Rcpp::Rcout << "ivcoef:" << std::endl;
		for(int ii =0;ii<ivblocksize;ii++) {
		  Rcpp::Rcout << "[" << ii << "]: " << *(rivcoefout.begin() + ivblocksize*hh + ii) << "; ";
		}
		Rcpp::Rcout << std::endl;
	      }
	      indcheck(pinds[n]-1,0,maxpind,"ivtpfilled bug (datachoicepar). n: %d; hh %d",n,hh);
	      if(ivtpfilled[pinds[n]-1] < 0) {
		ivtpfilled[pinds[n]-1] = n;
		for(int p=0;p<nrgrid;p++) {
		  ivtprobs[p+nrgrid*(pinds[n]-1)] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,rivvariout.begin()+(nd-1)*(nd-1)*hh,
							   rdetivout[hh],impdist[p]);
		  //tprobout(n+obsend*p) = ivtprobs[p];
		  if(n==printobs && p < 10) {
		  
		    Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p+nrgrid*(pinds[n]-1)] << ", " << impdist[p] <<
		      ", " << rdetivout[hh] << std::endl;
		    for(int ii=0;ii<nd-1;ii++) {
		      Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
		    }
		    Rcpp::Rcout << std::endl;
		  }
		}
	      }
	    }
	  }
	  if(prcheck) {Rcpp::Rcout << "2" << std::endl;}
        
	  double iv0[nsize];
	  for(int size=0;size<nsize;size++) {
	    iv0[size] = ivbig(n,size);
	    for(int j=0;j<1+capj;j++) {
	      ivgrid[j+size*(1+capj)] = iv0[size]*((double)j)/((double)capj);
	      if(n==printobs) {
		Rcpp::Rcout << "iv(" << n << "," << size << "): " << ivbig(n,size) <<
		  " iv0[" << size << "]: " << iv0[size] << " ivgrid[" <<
		  j+size*(1+capj) << "]: " << ivgrid[j+size*(1+capj)] << std::endl;
	      }
	    }
	  }
	  if(prcheck) {Rcpp::Rcout << "3" << std::endl;}
	  // compute value function at points that get visited - could fold this into next loop
	  int hhprint = -1, repprint = -1;
	  int printflag = hh == hhprint && rep == repprint;

	  //int pobs = 0; //n == 1 || n == obshhinds[hh+1]-2;

	  if(n==printobs) {
	    Rcpp::Rcout << "inv: " << inv << ", b: " << b << ", init inv: " << initinv[hh] << ", init b: " << initb[hh] << std::endl;
	  }
	
	  if(!myopic && beta > 0) {
	    for(int size=0;size<nsize;size++) {
	      for(int j=0;j<capj+1;j++) {
		if( (j == 0 && size == 0) || j > 0 ) {
		  int bprime1 = 0;
		  if(invmodel == 1) {
		    invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		    //drop[j+(1+capj)*size] = invprime > maxinv;
		    invprime = invprime > invub ? invub : invprime;
		  } else {
		    nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			    rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		    bprime1 = bprime;
		    bprime = bindex[bprime]-1;
		    indcheck(bprime,0,nb,"bprime (sim) error");
		  }

		  int iind = itoiind(invprime,istates,nistates,first20);
		  /*int iind = 0;
		    double dif = 0;
		    if(invprime > istates[first20-1]) {
		    dif = istates[first20] - istates[first20-1];
		    iind = first20 + (int)((invprime-istates[first20-1])/dif);
		    } else {
		    dif = istates[1] - istates[0];
		    iind = 1 + (int)((invprime-istates[0])/dif);
		    }*/
		  int vfindx1 = iind + nistates*(bprime+nbsmall*(pinds[n]-1));
                

		  // uncomment if not using price states
		  //vfindx1 = iind + nistates*bprime;
		  //if(vffilled[vfindx1] != n) {
		  int ok = iind < nistates;
		  if(invmodel == 2) {
		    ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
		  }
		  if(ok) {
		    indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
		    if(vffilled[vfindx1] < 0) {
		      vffilled[vfindx1] = n;
		      if(gridmethod==2) {
			vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						      nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*hh,
						      gvari.begin(), rdetivout[hh], dg, gmean.begin(),
						      nd-1, istates, nistates, first20, bprime, nbsmall, printflag,0,vfinterpmethod);

			if(n==printobs) {
			  //  Rcpp::Rcout << "nreptilde " << nreptilde << std::endl;
			  Rcpp::Rcout << "ivtprobs: " << *(ivtprobs+nrgrid*(pinds[n]-1)) << std::endl;
			  Rcpp::Rcout << "kweights: " << kweights[0] << std::endl;
			  Rcpp::Rcout << "vfin: " << vfsave[0] << std::endl;
			  Rcpp::Rcout << "vf output (" << n << ", " << iind << "): " << vftemp[iind+nistates*bprime] << std::endl;
			}

		      } else {
			vftemp[vfindx1] = vfavgijchh(ivout, iind, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						     nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*hh,
						     gvari.begin(), rdetivout[hh], dg, gmean.begin(),
						     nd-1, istates, nistates, first20, bprime, nbsmall, printflag);
		    }
		      if(n == printobs) {
			Rcpp::Rcout << " hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
			  " i': " << invprime << "b': " << bprime << " vf: " << vftemp[vfindx1] << std::endl;
		      }

		      //vfout(n+obsend*iind) = vftemp[vfindx1];

		    }
		  }

		  vfindx1 = iind-1 + nistates*(bprime+nbsmall*(pinds[n]-1));
		  indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug. hh %d, n %d",hh,n);
		  //vfindx1 = iind-1 + nistates*bprime;
		  //if(vffilled[vfindx1] != n) {
		  if(vffilled[vfindx1] < 0) {
		    vffilled[vfindx1] = n;
		    if(gridmethod == 2) {
		      vftemp[vfindx1] = vfavgijchh2(ivtprobs+nrgrid*(pinds[n]-1), iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder, kweights,
						    nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*hh,
						    gvari.begin(), rdetivout[hh], dg, gmean.begin(),
						    nd-1, istates, nistates, first20, bprime, nbsmall, 0, 0,vfinterpmethod);
		    } else {
		      vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfsave.begin()+vfind, rgridsave.begin(), indexes.begin(), coorder,
						   nreptilde, nsave, nrgrid, rivvariout.begin()+(nd-1)*(nd-1)*hh,
						   gvari.begin(), rdetivout[hh], dg, gmean.begin(),
						   nd-1, istates, nistates, first20, bprime, nbsmall, 0);
		    }

		    if(n == printobs) {
		      Rcpp::Rcout << " hh: " << hh+1 << " n: " << n << //" pinds[n]: " << pinds[n] <<
			" i': " << invprime << "b': " << bprime << " vf: " << vftemp[vfindx1] << std::endl;
		    }

		    //vfout(n+obsend*(iind-1)) = vftemp[vfindx1];

		  }

		}
	      }
	    }
	  }
	  if(prcheck) {Rcpp::Rcout << "4" << std::endl;}
	  int first = 1;

	  int vfoffset = nistates*nbsmall*(pinds[n]-1);

	  for(int size=0;size<nsize;size++) {
	    for(int j=0;j<capj+1;j++) {
	      int indx = j+(1+capj)*size;
	      drop[indx] = 1;
	      if( (j == 0 && size == 0) || j > 0 ) {
		drop[indx] = 0;
		if(invmodel == 1) {
		  invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
		  invprime = invprime > invub ? invub : invprime;
		  ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize[size], rcnew, gamma, nu, ccost, omega, cmodel, 0);
		} else {
		  nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			  rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,
                           &nextbflag,0);
		  indcheck(bprime,0,nb,"bprime (sim) error");
		  ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				   nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
		  bchoice[indx] = bprime;
		  bprime = bindex[bprime]-1;
		  ichoice[indx] = invprime;
		}
		if(!myopic && beta > 0) {
		  vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20,vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
		} else {
		  vfhat[indx] = 0.0;
		}
		indcheck(n+obsend*(indx),0,obsend*nsize*(1+capj),"vfout bug (datachoicepar). n: %d; hh %d",n,hh);
		vfout(n,indx) = vfhat[indx];
		utout(n,indx) = ut[indx];
		maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);

		if(n==printobs) {
		  Rcpp::Rcout << " ivgrid[" << indx << "]: " << ivgrid[indx] << " ut[" << indx << "]: " << ut[indx] << " vfhat[" << indx << "]: " << vfhat[indx] << std::endl;
		}

	      }
	    }
	  }
	  if(prcheck) {Rcpp::Rcout << "4.5" << std::endl;}
	  //if(printflag)
	  //  Rprintf("\n");

	  //ut[0] -= maxutil;
	  //cerror[0] = -log(-log(R::runif(0,1)));
	  double s1 = 0.0; //exp(ut[0]+vfhat[0]);

	  double choiceumax = -1000000;  //ut[0]+vfhat[0] + cerror[0];
	  indcheck(n+nobsdata*2,0,llinfo.nobs*3,"simchoice bug (2, datachoicepar). n: %d; hh %d",n,hh);
	  for(int j=0;j<3;j++) {
	    //simchoice(n+nobsdata*(j+3*rep)) = 0;
	    simchoice[n+nobsdata*j] = 0;
	  }

	  for(int size=0;size<nsize;size++) {
	    for(int j=0;j<capj+1;j++) {
	      int indx = j+(1+capj)*size;
	      if(!drop[indx]) {
		ut[indx] -= maxutil;
		//cerror[indx] = -log(-log(R::runif(0,1)));
		indcheck(n+nobsdata*(size+nsize*j),0,nobsdata*(1+capj)*nsize,"logiterror bug (datachoicepar). n: %d; hh %d",n,hh);
		cerror[indx] = logiterror(n,size+nsize*j);
		s1 += exp(ut[indx]+vfhat[indx]);
		if(ut[indx]+vfhat[indx]+cerror[indx] > choiceumax) {
		  choiceumax = ut[indx]+vfhat[indx]+cerror[indx];
                  if(invmodel==2) {
                    bprime = bchoice[indx];
                    invprime = ichoice[indx];
                  }
		  if(j > 0) {
		    //simchoice(n+nobsdata*(0+3*rep)) = size+1;
		    //simchoice(n+nobsdata*(1+3*rep)) = j;
		    simchoice[n+nobsdata*0] = size+1;
		    simchoice[n+nobsdata*1] = j;
		  }
		}
	      }
	    }
	  }
	  if(prcheck){Rcpp::Rcout << "4.6" << std::endl;}
	  // simulate brand choice if a size choice > 0 is chosen

	  //if(simchoice(n+nobsdata*(0+3*rep)) > 0) {
	  if(simchoice[n+nobsdata*0] > 0) {
	    choiceumax = -1000000;
	    double ut1;
	    double s=0;
	    double utnoerr = 0;
	    for(int j=0;j<nbrand;j++) {
	      //Rcpp::Rcout << "column index " << prindex << ", " << j << ", " << n << std::endl;
	      //double price = REAL(VECTOR_ELT(data,prindex-1+j))[n];
	      double price = cfdata(cfdataind,j);
	      if(n==printobs) {
		Rcpp::Rcout << "cfdataind: " << cfdataind << ", price: " << price << ", simchoice: "
			    << simchoice[n] << ", brsize: " << brsize[j+nbrand] << std::endl;
	      }
	      //Rcpp::Rcout << "y" << std::endl;
	      if(!(price >= 999 || simchoice[n] != brsize[j+nbrand])) {
		//Rcpp::Rcout << "z" << std::endl;
		//ut1 = xi[j] + alpha*price - log(-log(R::runif(0,1)));
		indcheck(n+nobsdata*j,0,nobsdata*nbrand,"logiterrorb bug (datachoicepar). n: %d; hh %d",n,hh);
		ut1 = xi[j] + alpha*price + logiterrorb(n,j);
		if(n==printobs) {
		  Rcpp::Rcout << "simchoice " << j << ": xi[j]: " << xi[j] << ", alpha: " << alpha <<
		    ", price: " << price << "error: " << logiterrorb(n,j) << std::endl;
		}
		if(ut1 > choiceumax) {
		  choiceumax = ut1;
		  utnoerr = xi[j] + alpha*price;
		  //simchoice(n+nobsdata*(2+3*rep)) = j+1;
		  simchoice[n+nobsdata*2] = j+1;
		}
		s+=exp(xi[j] + alpha*price);
	      }
	    }
	    llhhbr[hh] += log(exp(utnoerr)/s);
	  
	  }

	  if(n==printobs) {
	    Rcpp::Rcout << "simulated choice:" << std::endl;
	    Rcpp::Rcout << simchoice[n+nobsdata*0] << ", " << simchoice[n+nobsdata*1] << ", " <<
	      simchoice[n+nobsdata*2] << std::endl;
	  }
	  
	  if(prcheck) {Rcpp::Rcout << "4.7" << std::endl;}
	  double obspack = 0;
	  int obschoice = 0;

	  /*if(simchoice(n+nobsdata*(0+3*rep)) > 0) {
	    obspack = packsize(simchoice(n+nobsdata*(0+3*rep))-1);
	    obschoice = simchoice(n+nobsdata*(1+3*rep)) + (1+capj)*(simchoice(n+nobsdata*(0+3*rep))-1);
	    }*/

	  if(simchoice[n+nobsdata*0] > 0) {
	    obspack = packsize[simchoice[n+nobsdata*0]-1];
	    obschoice = simchoice[n+nobsdata*1] + (1+capj)*(simchoice[n+nobsdata*0]-1);
	  }
	  if(prcheck) {Rcpp::Rcout << "4.8" << std::endl;}
	  
	  if(prcheck) {Rcpp::Rcout << "4.85" << std::endl;}
	  
	  if(prcheck) {Rcpp::Rcout << "4.8" << std::endl;}
	  //indcheck(obschoice,0,(1+capj)*nsize,"ut. tunits: %d",((int)simchoice(n+nobsdata*(1+3*rep))));
	  indcheck(obschoice,0,(1+capj)*nsize,"ut. tunits: %d",((int)simchoice[n+nobsdata*1]));
	  //if( dmax(inv + obspack*((double)simchoice[n+nobsdata]) - rcnew,0) > invub ) {
	  //  llrun -= 1000000;
	  //} else {
	  //llrun += ut[obschoice] + vfhat[obschoice] - log(s1);
	  //}

	  llrun += log(exp(ut[obschoice] + vfhat[obschoice])/s1);

	  if(debug) {
	    //if( dmax(inv + obspack*((double)simchoice[n+nobsdata]) - rcnew,0) > invub ) {
	    //  rllrun[n] = -1000000;
	    //} else {
	    //rllrun[n] = exp(ut[obschoice]+vfhat[obschoice])/s1;
	    //}
	  }
	  risave2[n] = inv;
	  //bsave2(n) = b;
	  if(invmodel == 1) {
	    //inv = dmax(inv + obspack*((double)simchoice(n+nobsdata*(1+3*rep))) - rcnew,0);
	    inv = dmax(inv + obspack*((double)simchoice[n+nobsdata*1]) - rcnew,0);
	    inv = inv > invub ? invub : inv;
	  } else {
	    inv = invprime;
	    b = bprime;
	  }
	  n2++;
	}
	if(n==(llinfo.nobs-1) || panid[n] != panid[n+1]) {
	  llhh[hh] = llrun;  // I changed to return log likelihood rather than likelihood
	  ll += log(llhh[hh]);
	}
	if(prcheck) {Rcpp::Rcout << "5" << std::endl;}
      }
    }

  }

};

// [[Rcpp::export]]
Rcpp::List MCMCcounterfactual(Rcpp::NumericVector betadraw, Rcpp::List Mcmc, Rcpp::List info,
			      Rcpp::List hhbig, Rcpp::List hhmerge, Rcpp::List hhmergebr,
			      Rcpp::NumericMatrix cfdata, Rcpp::NumericMatrix cfinit, Rcpp::IntegerVector hhinclude, Rcpp::IntegerVector stvisitinit) {

  BEGIN_RCPP

  int maxrep = Rcpp::as< std::vector<int> >(Mcmc["R"])[0];
  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];
  int crfix = Rcpp::as< std::vector<int> >(info["crfix"])[0];
  int datacrate = Rcpp::as< std::vector<int> >(info["data.crate"])[0];
  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  int usecdraws = Rcpp::as< std::vector<int> >(info["usecdraws"])[0];
  int brandonly = Rcpp::as< std::vector<int> >(info["brandonly"])[0];
  int sizeshifter = 0;
  Rcpp::IntegerVector sizebrand(nbrand);
  if(info.containsElementNamed("sizeshifter")) {
    sizeshifter = Rcpp::as< std::vector<int> >(info["sizeshifter"])[0];
    SEXP sb1 = info["sizebrand"];
    Rcpp::IntegerVector sb2(sb1);
    sizebrand = Rcpp::clone(sb2);
  }
  int cfsim = 1;
  if(info.containsElementNamed("cfsim")) {
    cfsim = Rcpp::as< std::vector<int> >(info["cfsim"])[0];
  }

  int pexpect = 0;   // if this is 1, we'll use the original data to compute price expectations
  if(info.containsElementNamed("pexpect")) {
    pexpect = Rcpp::as< std::vector<int> >(info["pexpect"])[0];
  }

  
  Rcpp::NumericVector crate(nhhs);
  if(datacrate) {
    SEXP temp(info["crate"]);
    Rcpp::NumericVector temp1(temp);
    crate = Rcpp::clone(temp1);
  }

  SEXP tunitstemp1(hhbig["totunits"]);
  Rcpp::NumericVector tunits(tunitstemp1);

  SEXP panidtemp(hhbig["PANID"]);
  Rcpp::NumericVector panid(panidtemp);


  SEXP xfull1 = info["xfull"];
  Rcpp::NumericVector xfull(xfull1);

  int nhhinclude = hhinclude.size();
  Rcpp::IntegerVector includeflag(nhhs);

  for(int i=0;i<nhhinclude;i++) {
    includeflag(hhinclude(i)-1) = 1;
  }

  int npbig = xfull.size();

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = Rcpp::clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdrop2(ivdrop1);
    ivdrop = Rcpp::clone(ivdrop2);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = Rcpp::clone(ivdrop);

  SEXP tf1 = info["tform"];
  Rcpp::NumericVector tf(tf1);
  SEXP fixed1 = info["fixed"];
  Rcpp::IntegerVector fixed(fixed1);
  SEXP pequal1 = info["paramequal"];
  Rcpp::IntegerVector paramequal(pequal1);
  SEXP lbounds1 = info["lbounds"];
  Rcpp::NumericVector lbounds(lbounds1);
  SEXP ubounds1 = info["ubounds"];
  Rcpp::NumericVector ubounds(ubounds1);
  int crhhlb = Rcpp::as< std::vector<int> >(info["crhhlb"])[0];

  if(crhhlb) {
    if(info.containsElementNamed("crate")) {
      SEXP temp(info["crate"]);
      Rcpp::NumericVector temp1(temp);
      crate = Rcpp::clone(temp1);
    } else {
      throw std::range_error("crhhlb is TRUE but crate is undefined");
    }
  }

  //if(!inputs.containsElementNamed("pstart"))
  //  throw std::range_error("pstart is NULL in inputs");

  int npsmall = 0;
  for(int i = 0;i < npbig;i++) {
    npsmall += !fixed(i);
  }

  //if(datacrate) {
  //  pstart(nbrand,Rcpp::_) = crate;
  //}

  Rcpp::NumericMatrix pold(npsmall,nhhs);

  Rcpp::NumericMatrix poldtf(npbig,nhhs);

  // this will not copy the object

  SEXP tunitstemp(hhmergebr["totunits"]);
  Rcpp::NumericVector tunitsbr(tunitstemp);

  SEXP panidtemp1(hhmergebr["PANID"]);
  Rcpp::NumericVector panidbr(panidtemp1);

  SEXP brindextemp(hhmergebr["brindex"]);
  Rcpp::IntegerVector brindexbr(brindextemp);

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  int pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  int retindll = Rcpp::as< std::vector<int> >(info["retindll"])[0];
  int varflag = 0;
  if(info.containsElementNamed("varflag"))
    varflag = Rcpp::as< std::vector<int> >(info["varflag"])[0];

  SEXP panidmtemp(hhmerge["PANID"]);
  Rcpp::NumericVector panidmerge(panidmtemp);

  Rcpp::IntegerVector hhinds(2);
  hhinds(0) = 1;
  hhinds(1) = nhhs;

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  SEXP pindsmerge1(info["idpriceindmerge"]);
  Rcpp::IntegerVector pindsmerge(pindsmerge1);

  SEXP pindsbig1(info["idpriceindbig"]);
  Rcpp::IntegerVector pindsbig(pindsbig1);

  Rcpp::IntegerVector pinds(pindsbig.size());

  int hmaxpind = Rcpp::as< std::vector<int> >(info["hmaxpind"])[0];
  SEXP hpindsbig1(info["hpriceindbig"]);
  Rcpp::IntegerVector hpindsbig(hpindsbig1);

  int maxpind = Rcpp::as< std::vector<int> >(info["maxpind"])[0];
  int nobsbr = tunitsbr.size();
  int nobsmerge = panidmerge.size();

  int prindex = Rcpp::as< std::vector<int> >(info["prindex"])[0];
  Rcpp::NumericMatrix pricematbr(nobsbr,nbrand);
  SEXP ptemp = hhmergebr[prindex-1];
  Rcpp::NumericVector ptemp1(ptemp);
  for(int i=0;i<nbrand;i++) {
    ptemp = hhmergebr[prindex-1+i];
    ptemp1 = ptemp;
    pricematbr(Rcpp::_,i) = ptemp1;
  }

  // this copies - should just work with numericmatrix directly

  Rcpp::NumericMatrix pricematmerge(nobsmerge,nbrand);
  SEXP ptemp2 = hhmerge[prindex-1];
  Rcpp::NumericVector ptemp3(ptemp2);
  for(int i=0;i<nbrand;i++) {
    ptemp2 = hhmerge[prindex-1+i];
    ptemp3 = ptemp2;
    pricematmerge(Rcpp::_,i) = ptemp3;
  }


  // define inputs for quantity log likelihood
  llinputs llinfo;

  SEXP hhbigcol1 = hhbig[0];
  Rcpp::NumericVector hhbigc1(hhbigcol1);
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];

  llinfo.nb = nb;
  llinfo.nbsmall = nbsmall;
  llinfo.nco = npbig;
  llinfo.nobs = hhbigc1.size();
  llinfo.nbrand = nbrand;
  llinfo.nsize = nsize;
  llinfo.nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  llinfo.sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];
  llinfo.ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];
  llinfo.inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];
  int nsave = 10;
  if(Mcmc.containsElementNamed("nsave")) {
    llinfo.nsave = Rcpp::as< std::vector<int> >(Mcmc["nsave"])[0];
    nsave = llinfo.nsave;
  } else {
    llinfo.nsave = 10;
  }
  if(Mcmc.containsElementNamed("ntilde")) {
    llinfo.ntilde = Rcpp::as< std::vector<int> >(Mcmc["ntilde"])[0];
  } else {
    llinfo.ntilde = 3;
  }
  llinfo.nobsinit = stvisitinit.length();

  Rcpp::NumericMatrix pricematbig(llinfo.nobs,nbrand);
  SEXP ptemp4 = hhbig[prindex-1];
  Rcpp::NumericVector ptemp5(ptemp4);
  for(int i=0;i<nbrand;i++) {
    ptemp4 = hhbig[prindex-1+i];
    ptemp5 = ptemp4;
    pricematbig(Rcpp::_,i) = ptemp5;
  }

  int nitervf = 50; // number of iterations for value function update
  int niterstop = 500; // rep at which we switch to only a single update
  if(Mcmc.containsElementNamed("nitervf")) {
    nitervf = Rcpp::as< std::vector<int> >(Mcmc["nitervf"])[0];
  }

  if(Mcmc.containsElementNamed("niterstop")) {
    niterstop = Rcpp::as< std::vector<int> >(Mcmc["niterstop"])[0];
  }

  llinfo.nrgrid = Rcpp::as< std::vector<int> >(info["nrgrid"])[0];
  llinfo.necoef = Rcpp::as< std::vector<int> >(info["necoef"])[0];
  llinfo.capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];
  llinfo.nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];
  llinfo.ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];
  llinfo.retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];
  llinfo.myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];
  llinfo.cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];
  llinfo.debug = Rcpp::as< std::vector<int> >(info["debug"])[0];
  llinfo.usevf = 1;
  llinfo.idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];
  llinfo.first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];
  llinfo.hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];
  llinfo.genrnd = Rcpp::as< std::vector<int> >(info["genrnd"])[0];
  llinfo.nd = Rf_length(info["ngrid"]);
  llinfo.nq = INTEGER(Rf_getAttrib(info["qpts"], R_DimSymbol))[0];
  llinfo.ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];
  llinfo.pflag = Rcpp::as< std::vector<int> >(info["print"])[0];
  llinfo.hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  llinfo.repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0]-1;
  llinfo.ncpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[1];
  llinfo.rcpgrid = INTEGER(Rf_getAttrib(info["pgrid"], R_DimSymbol))[0];
  llinfo.ncrate = Rcpp::as< std::vector<int> >(info["ncost"])[0];
  llinfo.dg = Rcpp::as< std::vector<double> >(info["dg"])[0];
  llinfo.maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];
  llinfo.hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];
  llinfo.invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  llinfo.maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];
  llinfo.usecdraws = usecdraws;
  llinfo.initvf = 0;
  if(llinfo.invmodel == 1) {
    llinfo.lomega = 2;
  } else {
    llinfo.lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }
  int leniv = llinfo.nobs;
  int capj = llinfo.capj;
  int lomega = llinfo.lomega;
  int hmodel = llinfo.hmodel;
  double clb, cub, alpha, gamma, nu, beta, ccost;
  double xi[nbrand];
  double omega[lomega];
  int myopic = llinfo.myopic;
  double invub = llinfo.maxinv;
  double dg = llinfo.dg;
  int cmodel = llinfo.cmodel;
  Rcpp::NumericVector llhh(nhhs);
  Rcpp::NumericVector llhhbr(nhhs);
  int debug = llinfo.debug;
  if(hmodel) {
    pinds = Rcpp::clone(hpindsbig);
  } else {
    pinds = Rcpp::clone(pindsbig);
  }
  //throw std::range_error("stopping");
  double kweights[nsave];

  int vfind=0;
  int nco = npbig;

  int nd = Rf_length(info["ngrid"]);

  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);

  int linitx = initx.length();
  int ninitt = llinfo.ninitt;
  double maxinv = llinfo.maxinv;
  int hinvbound = llinfo.hinvbound;
  
  
  //Rcpp::NumericVector vfll(llinfo.nobs*nsize*(1+llinfo.capj));
  //Rcpp::NumericVector utll(llinfo.nobs*nsize*(1+llinfo.capj));


  if(llinfo.debug && llinfo.retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  Rcpp::NumericVector cosave(npbig*nhhs*llinfo.nsave);
  Rcpp::IntegerVector indexes(llinfo.nsave);

  for(int i=0;i<llinfo.nsave;i++) {
    indexes(i) = i+1;
  }

  SEXP tunitsbigtemp(hhbig["totunits"]);
  Rcpp::NumericVector tunitsbig(tunitsbigtemp);

  SEXP panidbigtemp(hhbig["PANID"]);
  Rcpp::NumericVector panidbig(panidbigtemp);

  SEXP brindexbigtemp(hhbig["brindex"]);
  Rcpp::IntegerVector brindexbig(brindexbigtemp);

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP expandbig1(info["expandbig"]);
  Rcpp::IntegerVector expandbig(expandbig1);

  Rcpp::NumericMatrix ivbig(llinfo.nobs,nsize);
  Rcpp::NumericMatrix ivbig1(llinfo.nobs,nsize);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);

  int maxgridlength = ngrid[ngrid.size()-1];

  int vfblocksize1 = llinfo.nsave*llinfo.nrgrid*nbsmall*maxgridlength;
  int vfblocksize2 = llinfo.nrgrid*nbsmall*maxgridlength;

  Rcpp::NumericVector vfsave(vfblocksize1*nhhs);

  if(info.containsElementNamed("vfstart")) {
    SEXP tempvfinfo(info["vfstart"]);
    Rcpp::NumericVector tempvf(tempvfinfo);
    vfsave = Rcpp::clone(tempvf);
    llinfo.initvf = 1;
  }

  // if initvf is true, copy in supplied cosave
  if(llinfo.initvf) {
    SEXP tempcosave1(info["cosave"]);
    Rcpp::NumericVector tempcosave(tempcosave1);
    cosave = Rcpp::clone(tempcosave);
  }

  //Rcpp::Rcout << "aaa" << std::endl;

  Rcpp::NumericVector vf(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);
  //throw std::range_error("stopping");
  // cdraws needs to be reddrawn so we redefine it here
  Rcpp::NumericVector cdraws(llinfo.nobs);

  // icdraws won't be used anymore except for debugging
  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);

  SEXP cdrawtemp(info["cdraws"]);
  Rcpp::NumericVector cdrawinfo(cdrawtemp);

  SEXP iiv1(info["initivdraws"]);
  Rcpp::NumericVector iiv(iiv1);

  SEXP ldraws1(info["logitdraws"]);
  Rcpp::NumericVector ldraws(ldraws1);

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);


  //SEXP initb1(info["initb"]);
  //Rcpp::IntegerVector initb(initb1);

  //SEXP inits1(info["inits"]);
  //Rcpp::IntegerVector inits(inits1);

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  SEXP qpts1(info["qpts"]);
  Rcpp::NumericVector qpts(qpts1);

  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);

  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);

  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);
  //throw std::range_error("stopping");
  int llbrandsize = nhhs;
  if(varflag)
    llbrandsize = nobsbr;

  // declarations of variables that get returned/work variables

  Rcpp::NumericVector llbrand(llbrandsize);
  Rcpp::NumericVector llbrand1(llbrandsize);
  Rcpp::NumericVector llbrand1a(llbrandsize);

  // the iv blocks that have 1 will just be dummies for when we compute expectations using the actual data, rather than the counterfactual data
  Rcpp::NumericMatrix rivout(nobsmerge,nsize);
  Rcpp::NumericMatrix rivout1(nobsmerge,nsize);

  //int cutblocksize = ncutiv == 0 ? 1 : (ncutiv+1)*nsize;
  int cutblocksize = (ncutiv+1)*nsize;
  Rcpp::NumericVector rcutmatout(cutblocksize*nhhs);
  Rcpp::NumericVector rcutmatout1(cutblocksize*nhhs);

  int nrowiv=0;
  if(sptype == 1) {
    nrowiv = 2*ncutiv+2;
  } else {
    if(ncutiv == 0) {
      nrowiv = nsize+1;
    } else {
      nrowiv = ncutiv+4;
    }
  }

  int ivcoefblocksize = ncutiv == 0 ? nrowiv*nsize : nrowiv*nsize*nsize;
  Rcpp::NumericVector rivcoefout(ivcoefblocksize*nhhs);
  Rcpp::NumericVector rivcoefout1(ivcoefblocksize*nhhs);

  int ivvariblocksize = nsize*nsize;
  Rcpp::NumericVector rivvariout(ivvariblocksize*nhhs);
  Rcpp::NumericVector rivvariout1(ivvariblocksize*nhhs);

  Rcpp::NumericVector rdetivout(nhhs);
  Rcpp::NumericVector rdetivout1(nhhs);

  // setup of grids, starting values for B and W,

  Rcpp::IntegerVector fixedpop(npsmall);
  if(Mcmc.containsElementNamed("fixedcoef")) {
    SEXP tempfixed(Mcmc["fixedcoef"]);
    Rcpp::IntegerVector tempfixed1(tempfixed);
    fixedpop = Rcpp::clone(tempfixed1);
  }
  //Rcpp::print(fixedpop);
  //Rcpp::Rcout << "abc" << std::endl;
  int nfixedpop = 0;
  for(int i=0;i<npsmall;i++) {
    nfixedpop += fixedpop(i);
  }
  //Rcpp::Rcout << "def" << std::endl;
  int allfixed = nfixedpop == npsmall;
  int allvarying = nfixedpop == 0;

  if(llinfo.hmodel && !allfixed) {
    Rcpp::Rcout << "Number of population fixed params: "  << nfixedpop << std::endl;
    throw std::runtime_error("if hmodel is true all parameters must be fixed");
  }


  int fixedmbound = nfixedpop;
  if(allvarying)
    fixedmbound = 1;

  // extract and define proposals

  arma::mat proposal(fixedmbound,fixedmbound);
  proposal.eye();
  if(Mcmc.containsElementNamed("proposal")) {
    SEXP tempprop1(Mcmc["proposal"]);
    Rcpp::NumericVector tempprop(tempprop1);
    if(tempprop.length() != nfixedpop*nfixedpop) {
      Rcpp::Rcout << "length of proposal: " << tempprop.length()
                 << "num fixed coefs: " << nfixedpop << std::endl;
      throw std::range_error("Incorrect dimensions for proposal.  Make sure the dimensions correspond to the number of population fixed coefficient.");
    }
    std::copy(tempprop.begin(),tempprop.end(),proposal.begin());
  }

  arma::mat cholprop = chol(proposal);

  // split up the proposal into 2 pieces

  int nbrfixed = 0;
  int ii1 = 0;
  for(int i=0;i<nbrand;i++) {
    if(!fixed(i)) {
      if(fixedpop(ii1)) {
        nbrfixed++;
      }
      ii1++;
    }
  }

  //Rcpp::Rcout << "nbrfixed: " << nbrfixed << " nfixedpop: " << nfixedpop << std::endl;

  int bend1 = nbrfixed;
  int bstart2 = nbrfixed;
  int bend2 = nfixedpop;

  if(nbrfixed == 0 || allvarying)
    bend1 = 1;

  arma::mat cpropbrand = cholprop.submat(0,0,bend1-1,bend1-1);

  if( nfixedpop == nbrfixed ) {
    bstart2 = nfixedpop-1;
    bend2 = nfixedpop;
  }

  if(allvarying) {
    bstart2 = 0;
    bend2 = 1;
  }

  arma::mat cpropother = cholprop.submat(bstart2,bstart2,bend2-1,bend2-1);

  arma::mat gvari1a(gvari.begin(),nsize,nsize,false);
  arma::mat gvar = arma::inv_sympd(gvari1a);
  arma::mat gchol = arma::chol(gvar);

  int nrgrid = llinfo.nrgrid;
  double rgridhh[(nd-1)*llinfo.nrgrid];
  
  // gridmethod == 1: old random grid of inclusive values
  // gridmethod == 2: fixed grid of price vectors (randomly drawn in advance)
  int gridmethod = 1;
  if(Mcmc.containsElementNamed("gridmethod")) {
    gridmethod = Rcpp::as< std::vector<int> >(Mcmc["gridmethod"])[0];
  }

  int vfinterpmethod = 1;
  if(Mcmc.containsElementNamed("vfinterpmethod")) {
    vfinterpmethod = Rcpp::as< std::vector<int> >(Mcmc["vfinterpmethod"])[0];
  }

  Rcpp::NumericMatrix bwmatrix(npbig,npbig);

  double hthumb = pow(4/(3*( (double)nsave )),0.2);
  for(int i=0;i<npbig;i++) {
    bwmatrix(i,i) = 1/hthumb;
  }

  if(Mcmc.containsElementNamed("bwmatrix")) {
    SEXP bwtemp1(Mcmc["bwmatrix"]);
    Rcpp::NumericMatrix bwtemp(bwtemp1);
    bwmatrix = Rcpp::clone(bwtemp);
  }

  int rgridbdbig = gridmethod == 2 ? nbrand*llinfo.nrgrid*llinfo.nsave : nsize*llinfo.nrgrid*llinfo.nsave;

  int rgridbd = gridmethod == 2 ? nbrand : nsize;

  Rcpp::NumericMatrix rgrid(rgridbd,nrgrid);
  Rcpp::NumericVector rgridsave(rgridbdbig);


  if(gridmethod == 1) {
    for(int i=0;i<nsize;i++) {
      rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
    }
  } else if (gridmethod == 2) {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  Rcpp::NumericVector impdist(nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }

  arma::mat rgrida(rgrid.begin(),rgridbd,nrgrid,false);
  arma::colvec gmeana(gmean.begin(),gmean.size(),false);

  arma::mat tempm(nsize,nsize);
  if(gridmethod != 2) {
    tempm = rgrida.t()*gchol;
    rgrida = tempm.t();
    rgrida.each_col() += gmeana;
  }

  int nistates = ngrid[llinfo.nd-1];

  int nkeep = Rcpp::as< std::vector<int> >(Mcmc["keep"])[0];
  int nwrite = Rcpp::as< std::vector<int> >(Mcmc["nwrite"])[0];

  int cfkeep = 10;
  if(Mcmc.containsElementNamed("cfkeep")) {
    cfkeep = Rcpp::as< std::vector<int> >(Mcmc["cfkeep"])[0];
  }
  int burnin = (int)maxrep/2;
  if(Mcmc.containsElementNamed("burnin")) {
    burnin = Rcpp::as< std::vector<int> >(Mcmc["burnin"])[0];
  }

  Rcpp::Environment base("package:base");
  Rcpp::Function save = base["save"];

  int ndrawkeep = maxrep/nkeep;

  // compute number of demographic vars

  SEXP dvarname1(info["dvarnames"]);
  Rcpp::StringVector dvarnames(dvarname1);

  int ndvars = dvarnames.size();

  Rcpp::NumericMatrix bdraw(ndrawkeep,npsmall*ndvars);

  Rcpp::NumericMatrix Wdraw(ndrawkeep,npsmall*npsmall);

  Rcpp::NumericVector llbrsave(ndrawkeep);
  Rcpp::NumericVector llqsave(ndrawkeep);

  // make these TRANSFORMED beta draws

  Rcpp::NumericMatrix W(npsmall,npsmall);


  for(int i=0;i<npsmall;i++) {
    if(!fixedpop[i])
      W(i,i) = 1;
  }
  //Rcpp::Rcout << "ndvars: " << ndvars << std::endl;

  /*Rcpp::NumericMatrix b(npsmall,ndvars);

  for(int i=0;i<npsmall;i++) {
    b(i,0) = Rcpp::mean(pold(i,Rcpp::_));
    }*/

  // for speed, make a matrix of the demographic variables

  SEXP dvar1(info["dvars"]);
  Rcpp::NumericMatrix dvars(dvar1); //check that dims are right
  //Rcpp::Rcout << "dvar rows: " << dvars.nrow() << ". cols: " << dvars.ncol() << std::endl;

  arma::mat W1(npsmall-nfixedpop,npsmall-nfixedpop);
  W1.zeros();

  arma::mat W1chol(npsmall-nfixedpop,npsmall-nfixedpop);
  W1chol.zeros();

  arma::mat W1inv(npsmall-nfixedpop,npsmall-nfixedpop);
  W1inv.zeros();

  arma::mat bnotfixed(npsmall-nfixedpop,nhhs);
  bnotfixed.zeros();

  arma::mat temppnew(npsmall-nfixedpop,nhhs);
  temppnew.zeros();

  arma::mat temppold(npsmall-nfixedpop,nhhinclude);
  temppold.zeros();

  // this is only needed for conformability
  arma::mat temppold1(npsmall-nfixedpop,nhhs);
  temppold1.zeros();

  Rcpp::NumericVector bdiff0(nhhs);
  Rcpp::NumericVector bdiff1(nhhs);

  arma::mat bdiff0a(bdiff0.begin(),bdiff0.size(),false);
  arma::mat bdiff1a(bdiff1.begin(),bdiff1.size(),false);

  arma::mat betamean(npsmall-nfixedpop,1);
  betamean.zeros();
  Rcpp::Rcout << "ndvars: " << ndvars << std::endl;
  arma::mat b1(npsmall-nfixedpop,ndvars);
  b1.zeros();
  Rcpp::Rcout << b1.n_rows << ", " << b1.n_cols << std::endl;
  arma::mat b2(nfixedpop,1);
  b2.zeros();

  arma::mat b2a(nbrfixed,1);
  b2a.zeros();

  arma::mat b2b(nfixedpop-nbrfixed,1);
  b2b.zeros();
  Rcpp::Rcout << "b2b size " << nfixedpop-nbrfixed << std::endl;

  arma::mat btempdiff(npsmall-nfixedpop,1);
  btempdiff.zeros();

  arma::mat bmeandiff(npsmall-nfixedpop,1);
  bmeandiff.zeros();

  Rcpp::NumericVector pnewdrawvec((npsmall-nfixedpop)*nhhs);
  arma::mat pnewdrawmat(pnewdrawvec.begin(),nhhs,npsmall-nfixedpop,false);

  Rcpp::NumericVector llqhh(nhhs);
  Rcpp::NumericVector llqhh1(nhhs);
  Rcpp::NumericVector llqhh1a(nhhs);
  Rcpp::NumericVector rllrun(llinfo.nsim*llinfo.nobs);
  Rcpp::NumericVector risave1(llinfo.nsim*nhhs*llinfo.ninitt);
  Rcpp::NumericVector risave2(llinfo.nsim*llinfo.nobs);

  Rcpp::NumericMatrix pnew(npsmall,nhhs);
  Rcpp::NumericMatrix pnewtf(npbig,nhhs);

  Rcpp::NumericVector r(nhhs);
  Rcpp::NumericVector u(nhhs);

  arma::mat npsmallnf(npsmall-nfixedpop,npsmall-nfixedpop);
  npsmallnf.eye();
  npsmallnf *= npsmall-nfixedpop;

  arma::mat S(npsmall-nfixedpop,npsmall-nfixedpop);
  S.zeros();

  arma::mat Schol(npsmall-nfixedpop,npsmall-nfixedpop);
  Schol.zeros();

  arma::mat XX(npsmall-nfixedpop,npsmall-nfixedpop);
  XX.zeros();

  arma::mat tempb1(1,npsmall - nfixedpop);
  tempb1.zeros();

  arma::mat tempb2(1,nfixedpop);
  tempb2.zeros();

  arma::mat tempb2a(1,nbrfixed);
  tempb2a.zeros();

  arma::mat tempb2b(1,nfixedpop-nbrfixed);
  tempb2b.zeros();

  // varying coefficient rho (1st el is brands, and second is other params)
  Rcpp::NumericVector rho(2);
  if(Mcmc.containsElementNamed("rho")) {
    SEXP temprho(Mcmc["rho"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho in Mcmc must be 2");
    rho = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho(i) = 0.1;
  }

  // fixed coefficients rho1
  Rcpp::NumericVector rho1(2);
  if(Mcmc.containsElementNamed("rho1")) {
    SEXP temprho(Mcmc["rho1"]);
    Rcpp::NumericVector temprho1(temprho);
    if(temprho1.size() < 2)
      throw std::range_error("length of rho1 in Mcmc must be 2");
    rho1 = Rcpp::clone(temprho1);
  } else {
    for(int i=0;i<2;i++)
      rho1(i) = 0.1;
  }

  // saves of rho1 - so we can see when it stabilizes
  Rcpp::NumericMatrix rho1save(ndrawkeep,2);

  int splitfixed = 1;
  if(Mcmc.containsElementNamed("splitfixed")) {
    splitfixed = Rcpp::as< std::vector<int> >(Mcmc["splitfixed"])[0];
  }


  int npropsave = 100;
  if(Mcmc.containsElementNamed("npropsave")) {
    npropsave = Rcpp::as< std::vector<int> >(Mcmc["npropsave"])[0];
  }

  int nprint = 10;
  if(Mcmc.containsElementNamed("nprint")) {
    nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  }

  std::vector<std::string> pnamevec(npbig,"V");
  for(int i=0;i<npbig;i++) {
    std::stringstream ss;
    ss << i+1;
    pnamevec[i] += ss.str();
  }

  Rcpp::CharacterVector pnames = Rcpp::wrap(pnamevec);

  if(Mcmc.containsElementNamed("pnames")) {
    nprint = Rcpp::as< std::vector<int> >(Mcmc["nprint"])[0];
  }

  Rcpp::IntegerMatrix xaccept(npropsave,2);
  Rcpp::IntegerVector acceptflag(nhhs);

  int wmethod = 2;  //inverse wishart drawing method - all seem about the same
  // 1: Rossi code (with demographics)
  // 2: my code (derived from ken train's R code)
  // 3: inverse gamma

  // priors

  // variance matrix prior

  arma::mat S0 = arma::eye(npsmall-nfixedpop,npsmall-nfixedpop);
  //int nu = npsmall-nfixedpop+3;
  if(Mcmc.containsElementNamed("S0")) {
    SEXP temp(Mcmc["S0"]);
    Rcpp::NumericMatrix tempS(temp);
    for(int i=0;i<npsmall-nfixedpop;i++) {
      for(int j=0;j<npsmall-nfixedpop;j++) {
        S0(i,j) = tempS(i,j);
      }
    }
  }
  // not used here
  /*if(Mcmc.containsElementNamed("nu")) {
    nu = Rcpp::as< std::vector<int> >(Mcmc["nu"])[0];
  }
  S0 = S0*( (double)nu );*/

  //arma::mat Stemp(npsmall-nfixedpop,npsmall-nfixedpop+nhhs);
  arma::mat Stemp(npsmall-nfixedpop,nu+nhhinclude);
  Stemp.zeros();
  
  // mean parameter priors

  arma::mat bA = arma::eye(ndvars,ndvars);

  bA = bA*1000;
  if(Mcmc.containsElementNamed("bA")) {
    SEXP temp(Mcmc["bA"]);
    Rcpp::NumericMatrix tempbA(temp);
    bA = Rcpp::as<arma::mat>(tempbA);  // will clone
  }

  // precomputing to save time
  arma::mat RA = arma::chol(bA);
  arma::mat dvarm = Rcpp::as<arma::mat>(dvars);
  arma::mat dvarminc(nhhinclude,ndvars);
  ii1 = 0;
  for(int i=0;i<nhhs;i++) {
    if(includeflag(i)) {
      dvarminc.row(ii1) = dvarm.row(i);
      ii1++;
    }
  }
  arma::mat WW = arma::join_cols(dvarminc,RA);

  arma::mat betabar(npsmall-nfixedpop,ndvars);
  betabar.zeros();

  if(Mcmc.containsElementNamed("betabar")) {
    SEXP temp(Mcmc["betabar"]);
    Rcpp::NumericMatrix tempbb(temp);
    betabar = Rcpp::as<arma::mat>(tempbb);
  }
  arma::mat betabarvec = betabar;
  betabarvec.reshape(ndvars*(npsmall-nfixedpop),1);
  
  //this is a fudge when we have no demographic vars
  //We need to do it right with demographics
  double betaVA = 10.0;
  arma::mat betaA(npsmall-nfixedpop,npsmall-nfixedpop);
  betaA.eye();
  betaA = betaA*betaVA;
  if(Mcmc.containsElementNamed("betaA")) {
    SEXP temp(Mcmc["betaA"]);
    Rcpp::NumericMatrix tempbetaA(temp);
    betaA = Rcpp::as<arma::mat>(tempbetaA);  // will clone
  }
  arma::mat betaS0I = arma::inv_sympd(betaA);

  //nu += nhhs;
  nu += nhhinclude;

  int naccept = 0;
  
  int useparamstart = 0;
  Rcpp::NumericMatrix paramstart(npbig,nhhs);

  tform(pold,xfull,tf,fixed,paramequal,paramstart,lbounds,ubounds,crate,useparamstart,
  	nbrand,npsmall,npbig,nhhs,crfix,datacrate,crhhlb,sizeshifter,nsize,brsize.begin(),
        sizebrand.begin(),poldtf);

  //Rcpp::NumericVector ichoiceout(initx.length());
  
  // choices at initial values
  /*Rcpp::IntegerVector ichoicej(initx.length()*maxrep);
  Rcpp::IntegerVector ichoices(initx.length()*maxrep);
  Rcpp::IntegerVector ibsave(initx.length()*maxrep);*/
  Rcpp::IntegerVector ichoicej(initx.length());
  Rcpp::IntegerVector ichoices(initx.length());
  Rcpp::IntegerVector ibsave(initx.length());
  
  // choices at estimation data
  int obsend;
  if(hhinds(1)==nhhs) {
    obsend = llinfo.nobs;
  } else {
    obsend = obshhinds(hhinds(1));
  }

  int nobsdata = obsend - (obshhinds(hhinds(0)-1)-1);
  
  //Rcpp::IntegerVector simchoice(nobsdata*3*maxrep);
  Rcpp::IntegerVector simchoice(nobsdata*3);
  
  //Rcpp::NumericVector risave2(nobsdata);
  Rcpp::IntegerVector bsave2(nobsdata);

  // output summary statistics

  int cfrepkeep = (int)((maxrep-burnin)/cfkeep);
  Rcpp::IntegerVector savechoice(nobsdata*3*cfrepkeep);
  Rcpp::NumericMatrix saveinv(llinfo.nobs,cfrepkeep);

  llinfo.cfrepkeep = cfrepkeep;

  Rcpp::NumericMatrix revenues(maxrep,nbrand);
  Rcpp::NumericMatrix quantities(maxrep,nbrand);
  Rcpp::NumericMatrix nunits(maxrep,1+capj);
  
  Rcpp::Rcout << "starting simulation loop" <<std::endl;

  
  int maxbottles = llinfo.maxbottles;
  double ut[(1+capj)*nsize];
  int drop[(1+capj)*nsize];
  double vfhat[(1+capj)*nsize];
  double cerror[(1+capj)*nsize];
  int bchoice[(1+capj)*nsize];
  double ichoice[(1+capj)*nsize];

  
  double ivout[nsize];
  double ivgrid[(1+capj)*nsize];
  double vftemp[nbsmall*nistates*maxpind];
  int vffilled[nbsmall*nistates*maxpind];

  int nrep, nreptilde;
  
  Rcpp::NumericVector initinv(nhhs);
  Rcpp::IntegerVector initb(nhhs);

  
  int coorder[nsave];
  
  
  int ncut = llinfo.ncut;
  int inttype = llinfo.inttype;
  int ivblocksize;
  if(sptype == 1) {
    if(ncutiv == 0) {
      ivblocksize = (2*ncutiv+2)*nsize;
    } else {
      ivblocksize = (2*ncutiv+2)*nsize*nsize;
    }
  } else {
    if(ncutiv == 0) {
      ivblocksize = (nsize + 1)*nsize;
    } else {
      ivblocksize = (ncutiv+4)*nsize*nsize;
    }
  }

  double ivtprobs[nrgrid*maxpind];
  int ivtpfilled[maxpind];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }
  int first20 = llinfo.first20;

  //Rcpp::NumericVector detiv(nhhs);
  int invmodel = llinfo.invmodel;

  int keep=-1;

  // predefine initial indexes for choice simulations
  Rcpp::IntegerVector iobsind(nhhs);
  Rcpp::IntegerVector iobsindbig(nhhs);
  int iobs = -1;
  int iobsbig = 0;
  for(int n=hhinds(0)-1;n<hhinds(1);n++) {
    iobsind(n) = iobs;
    iobsindbig(n) = iobsbig;
    for(int t=0;t<ninitt;t++) {
      if(stvisitinit(iobsbig)) {
	iobs++;
      }

      iobsbig++;
    }
  }

  Rcpp::NumericMatrix logiterror(nhhs*ninitt,(1+capj)*nsize);
  Rcpp::NumericMatrix logiterrorobs(nobsdata,(1+capj)*nsize);
  Rcpp::NumericMatrix logiterrorb(nobsdata,nbrand);

  // value function, utilities we save to look at
  Rcpp::NumericMatrix vfout(obsend,nsize*(1+capj));
  Rcpp::NumericMatrix utout(obsend,nsize*(1+capj));

  if(usecdraws) {
    Rcpp::Rcout << "usecdraws is TRUE. cdraws and icdraws will be copied from info." << std::endl;
  } else {
    Rcpp::Rcout << "usecdraws is FALSE. cdraws and icdraws will be drawn." << std::endl;
  }

  Rcpp::NumericVector rcdraws(nhhs*ninitt);

  if(usecdraws) {
    cdraws = Rcpp::clone(cdrawinfo);
    rcdraws = Rcpp::clone(ic);
    llinfo.lencdraws = cdraws.length();
  } else {
    llinfo.lencdraws = llinfo.nobs;
  }
  llinfo.lenicdraws = rcdraws.length();
  
  //throw std::range_error("stopping");
  for(int rep=0;rep<maxrep;rep++) {
    if(rep >= burnin && rep%cfkeep == 0) {
      keep++;
    }
    // initializations
    //int keep=0;

    Rcpp::Timer timer;

    for(int cc=0;cc<cfsim;cc++) {
      double inv = 0, invprime = 0;
      int b=0, bprime=0;
      for(int i=0;i<nhhs;i++) {
	initb(i) = 0;
      }
      memset(coorder,0,nsave*sizeof(int));
    
      double ll = 0;

      //memset(vffilled,-1,nbsmall*nistates*maxpind*sizeof(int));
      //memset(ivtpfilled,-1,maxpind*sizeof(int));
      for(int i=0;i<nbsmall*nistates*maxpind;i++) {
	vffilled[i] = -1;
      }
      for(int i=0;i<maxpind;i++) {
	ivtpfilled[i] = -1;
      }
    
      if(rep <= nsave) {
	nrep = rep-1;
      } else {
	nrep = nsave;
      }
      int ntilde = llinfo.ntilde;
      if(nrep < ntilde) {
	nreptilde = rep-1;
      } else {
	nreptilde = ntilde;
      }

      int prcheck = 0;

      //debugprint = rep == 1334;
      //if(rep==0) {
      if(!usecdraws) {
	//cdraws = Rcpp::runif(llinfo.nobs);
	for(int i=0;i<llinfo.nobs;i++) {
	  cdraws(i) = R::runif(0,1);
	}
      }
      //}


      // initial likelihood calculations
      if(prcheck)
	Rcpp::Rcout << "rep " << rep << std::endl;

      int j=0;
      for(int i=0;i<npbig;i++) {
	if(!fixed(i)) {
	  for(int hh=0;hh<nhhs;hh++) {
	    poldtf(i,hh) = betadraw(rep+maxrep*(j+npsmall*hh));
	    if(sizeshifter && i < nbrand && brsize(i+nbrand) > 1 && sizebrand(i)) {
	      poldtf(i,hh) -= betadraw(rep+maxrep*(npsmall-1-nsize+brsize(i+nbrand)+npsmall*hh));
	    }
	  }
	  j++;
	} else if(paramequal(i) > 0) {
	  if(paramequal(i)-1 >= i)
            throw std::runtime_error("Error (counterfactual) paramequal[i] must be less than i");
	  for(int hh=0;hh<nhhs;hh++) {
	    poldtf(i,hh) = poldtf(paramequal(i)-1,hh);
	  }
	}
	if(sizeshifter && i < nbrand && brsize(i+nbrand) > 1 && sizebrand(i)) {
	  for(int hh=0;hh<nhhs;hh++) {
	    poldtf(i,hh) += betadraw(rep+maxrep*(npsmall-1-nsize+brsize(i+nbrand)+npsmall*hh));
	  }
	}
      }

      // compute inclusive values

      if(hmodel) {
	fitiv(poldtf,  panidmerge, cfdata, hhinds, brsize,
	      obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
	      ncutiv, rivout, rcutmatout, rivcoefout, rivvariout, rdetivout);
	if(pexpect) {
	  fitiv(poldtf,  panidmerge, pricematmerge, hhinds, brsize,
		obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
		ncutiv, rivout1, rcutmatout, rivcoefout, rivvariout, rdetivout);
	}
      } else {
	fitivhh(poldtf,  panidmerge, cfdata, hhinds, brsize,
		obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
		ncutiv, dropmat, ivdrop, rivout, rcutmatout, rivcoefout,
		rivvariout, rdetivout,ivdrop2);
	if(pexpect) {
	  fitivhh(poldtf,  panidmerge, cfdata, hhinds, brsize,
		  obshhindsmerge, npbig, nobsmerge, nbrand, nsize, nhhs, sptype,
		  ncutiv, dropmat, ivdrop, rivout1, rcutmatout, rivcoefout,
		  rivvariout, rdetivout,ivdrop2);
	}
      }
      if(prcheck) {Rcpp::Rcout << "iv computed" << std::endl;}
      //throw std::range_error("stopping");
      for(int i=0;i<llinfo.nobs;i++) {
	for(int j=0;j<nsize;j++) {
	  ivbig(i,j) = rivout(expandbig[i]-1,j);
	}
      }

      // simulate initial choices
      // too slow, need to index price observations here

      if (usecdraws) {
	//rcdraws = Rcpp::runif(nhhs*ninitt);
	for(int i=0;i<nhhs*ninitt;i++) {
	  rcdraws(i) = R::runif(0,1);
	}
      }
    
    
      if(prcheck) {Rcpp::Rcout << "abc" << std::endl;}
      //logiterror = -log(-log(Rcpp::runif(nhhs*ninitt*(1+capj)*nsize)));  //I had problems with this.
      //Rcpp::Rcout << "size: " << nhhs*ninitt*(1+capj)*nsize << std::endl;
      for(int i=0;i<nhhs*ninitt;i++) {
	for(int j=0;j<(1+capj)*nsize;j++) {
	  logiterror(i,j) = -log(-log(R::runif(0,1)));
	}
	//Rcpp::Rcout << " " << i;
      }
      if(prcheck) {Rcpp::Rcout << "cde" << std::endl;}
      initchoicepar initpar(poldtf, iobsind, iobsindbig, rgrid, brsize, rcdraws,
			    initx,  stvisitinit, cfinit, cosave, indexes, bwmatrix,
			    rivcoefout, rcutmatout, rivvariout, rdetivout,
			    impdist, packsize, bstates, revarray, badmissable, bindex,
			    pgrid, vfsave, rgridsave, gvari, gmean, logiterror, ivdrop2,
			    istates, rep,
			    nistates, maxpind, nd, linitx, vfinterpmethod, ivblocksize,
			    cutblocksize, gridmethod, llinfo, ichoicej, ichoices,
			    initinv, initb, risave1, ibsave);

      parallelFor(hhinds[0]-1,hhinds[1],initpar);
      //initpar(hhinds[0]-1,hhinds[1]);
    
      // simulate choices at vf estimate
      // update value function

      //Rcpp::Rcout << "finished initial simulation" << std::endl;
      //timer.step("Iteration");
      //Rcpp::NumericVector tt1(timer);

      //tt1[0] /= 1e9;
      //Rcpp::Rcout << "Elapsed Seconds: "  << tt1 << std::endl;
      vfind = 0;

      int hh = hhinds[0]-2;
      double llrun = 0;

      //double invub = maxinv;

      for(int i=0;i<nobsdata*nsize;i++) {
	for(int j=0;j<(capj+1);j++) {
	  logiterrorobs(i,j) = -log(-log(R::runif(0,1)));
	}
      }

      for(int i=0;i<nobsdata;i++) {
	for(int j=0;j<nbrand;j++) {
	  logiterrorb(i,j) = -log(-log(R::runif(0,1)));
	}
      }

      datachoicepar datapar(poldtf,obshhindsmerge,obshhinds,rgrid,brsize,
			    cdraws,initx,stvisitinit,cfdata,cosave,indexes,
			    bwmatrix,rivcoefout,rcutmatout,rivvariout,
			    rdetivout,impdist,packsize,bstates,revarray,
			    badmissable,bindex,pgrid,vfsave,rgridsave,
			    gvari,gmean,logiterrorobs,ivdrop2,tunits,
			    ivbig,initinv,initb,pinds,logiterrorb,
			    panid,hhinds,istates,
			    rep,nistates,maxpind,nd,linitx,vfinterpmethod,
			    ivblocksize,cutblocksize,gridmethod,llinfo,
			    burnin,cfkeep,keep,nobsdata,leniv,simchoice,
			    vfout,utout,risave2);
			  
      parallelFor(hhinds[0]-1,hhinds[1],datapar);
      //datapar(hhinds[0]-1,hhinds[1]);

      int cfdataind = -1;

      for(int n=0;n<nobsdata;n++) {

	int units = (int)tunits(n);

	if(units >= 0) {

	  cfdataind++;

	  if(simchoice(n+nobsdata*2) > 0) {
	    nunits(rep,simchoice(n+nobsdata*1)) += 1.0/((double)cfsim);
	    quantities(rep,simchoice(n+nobsdata*2)-1) += simchoice(n+nobsdata*1)/((double)cfsim);
	    revenues(rep,simchoice(n+nobsdata*2)-1) += cfdata(cfdataind,simchoice(n+nobsdata*2)-1)*((double)simchoice(n+nobsdata*1))/((double)cfsim);
	  } else {
	    nunits(rep,0) += 1.0/((double)cfsim);
	  }

	}
      
      }

    }

    // save output
    if(rep >= burnin && rep%cfkeep == 0) {
      for(int n=0;n<nobsdata;n++) {
        //indcheck(n+nobsdata*(2+3*keep),0,nobsdata*3*llinfo.cfrepkeep,"(datachoicepar) savechoice bug, n %d, keep %d",n,keep);
        for(int i=0;i<3;i++) {
          savechoice(n+nobsdata*(i+3*keep)) = simchoice(n+nobsdata*i);
        }
      }
      for(int n=0;n<llinfo.nobs;n++) {
	saveinv(n,keep) = risave2(n);
      }
    }
    
    // importance probs - for checking
    //Rcpp::NumericVector tprobout(obsend*nrgrid);

    // update the value function

    int indx1 = rep%nsave;

    if(!llinfo.myopic && !brandonly) {
      if(rep == niterstop)
        nitervf = 1;
      //Rcpp::Rcout << "Updating VF" << std::endl;
      vfupdateijc(hhinds,bstates,pgrid,ngrid,poldtf,vfsave,rgrid,cosave,rgridsave,
                  packsize,rcutmatout,rivcoefout,rivvariout,rdetivout,gvari,gmean,
                  indexes,brsize,impdist,bwmatrix,revarray,badmissable,bindex,ivdrop2,
		  llinfo,rep,nitervf,gridmethod,
                  vfinterpmethod,vf);


      // see if can speed this up w std::copy - will need to reindex
      for(int hh=0;hh<nhhs;hh++) {
        for(int i=0;i<vfblocksize2;i++) {
          indcheck(hh*vfblocksize1 + indx1 + nsave*i,0,vfsave.length(),"vfsave copying bug: hh: %d; rep: %d; i: %d; vfblocksize1: %d; nsave: %d",hh,rep,i,
                   vfblocksize1,nsave);
          indcheck(hh*vfblocksize2+i,0,vf.length(),"vf copying bug: hh: %d; rep: %d; i: %d",hh,rep,i);
          vfsave[hh*vfblocksize1 + indx1 + nsave*i] = vf[hh*vfblocksize2+i];
        }
      }

    }

    for(int s=0;s<rgridbd;s++) {
      for(int j=0;j<nrgrid;j++) {
        indcheck(s+rgridbd*(j+nrgrid*indx1),0,rgridsave.length(),"rgridsave copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        indcheck(s+rgridbd*j,0,rgridsave.length(),"rgrid copying bug: s: %d; rep: %d; j: %d",s,rep,j);
        rgridsave[s+rgridbd*(j+nrgrid*indx1)] = rgrid[s+rgridbd*j];
      }
    }

    for(int hh=0;hh<nhhs;hh++) {
      int pind = hh*npbig*nsave;
      for(int j=0;j<npbig;j++) {
        indcheck(j+pind+indx1*npbig,0,cosave.length(),"cosave copying bug: hh: %d; rep: %d; j: %d",hh,rep,j);
        indcheck(j+npbig*hh,0,poldtf.length(),"poldtf copying bug: hh: %d; rep: %d; j: %d",hh,rep,j);
        cosave[j+pind+indx1*npbig] = poldtf[j+npbig*hh];
      }
    }

    if(gridmethod != 2) {
      for(int i=0;i<rgridbd;i++) {
	rgrid(i,Rcpp::_) = Rcpp::rnorm(nrgrid);
      }

      tempm = rgrida.t()*gchol;
      rgrida = tempm.t();
      rgrida.each_col() += gmeana;

    }

    //if(((rep+1)%nprint == 0 && nprint > 1) || nprint == 1) {
      Rcpp::Rcout << "Finished counterfactual iteration " << rep+1 << ":" << std::endl;

      timer.step("Iteration");

      Rcpp::NumericVector tt(timer);
      tt[0] /= 1e9;
      //tt[1] /= 1e9;
      Rcpp::Rcout << "Elapsed Seconds: "  << tt << std::endl;
      Rcpp::Rcout << std::endl;

      Rcpp::Rcout << "Average Value Function:" << std::accumulate(vf.begin(),vf.end(),0.0)/((double)vf.size()) << std::endl;
      Rcpp::Rcout << "Inventory stats: " << std::endl;
      Rcpp::Rcout << "Min: " << Rcpp::min(risave2) << ", median: " << Rcpp::median(risave2) <<
	", mean: " << Rcpp::mean(risave2) << "max: " << Rcpp::max(risave2) << ", sd: " <<
	Rcpp::sd(risave2) << std::endl;
      Rcpp::Rcout << "Total quantity: " << Rcpp::sum(quantities(rep,Rcpp::_)) << std::endl;
      Rcpp::Rcout << "Total revenue: " << Rcpp::sum(revenues(rep,Rcpp::_)) << std::endl;

      Rcpp::Rcout << std::endl;

      /*base["initinv"] = initinv;
      base["vf"] = vf;
      base["quantities"] = quantities;
      base["revenues"] = revenues;
      base["utout"] = utout;
      base["vfout"] = vfout;
      base["logiterrorobs"] = logiterrorobs;
      base["logiterrorb"] = logiterrorb;
      base["simchoice"] = simchoice;
      base["initinv"] = initinv;
      base["risave2"] = risave2;
      base["poldtf"] = poldtf;

      throw std::runtime_error("stopping.");*/

      //}

  }


  Rcpp::List ret;
  ret["poldtf"] = poldtf;
  //ret["pnewtf"] = pnewtf;
  ret["llbrand"] = llhhbr;
  //ret["llbrand1"] = llbrand1;
  ret["riv"] = rivout;
  ret["rivcoef"] = rivcoefout;
  ret["rivvari"] = rivvariout;
  ret["rdetiv"] = rdetivout;
  ret["llqhh"] = llhh;
  //ret["llqhh1"] = llqhh1;
  //ret["naccept"] = naccept;
  //ret["r"] = r;
  //ret["b"] = b;
  //ret["W"] = W;
  //ret["acceptflag"] = acceptflag;
  //ret["llqhh1a"] = llqhh1a;
  //ret["llbrand1a"] = llbrand1a;
  //ret["xaccept"] = xaccept;
  ret["vf"] = vf;
  //ret["bdraw"] = bdraw;
  //ret["Wdraw"] = Wdraw;
  //ret["betadraw"] = betadraw;
  //ret["rho1"] = rho1save;
  //ret["pstart"] = pstart;

  Rcpp::NumericVector vfout1(llinfo.nrgrid*nbsmall*maxgridlength*nhhs);

  for(int hh=0;hh<nhhs;hh++) {
    for(int j=0;j<nrgrid;j++) {
      for(int i=0;i<nistates;i++) {
        for(int b=0;b<nbsmall;b++) {
          vfout1(hh + nhhs*(j+nrgrid*(i+nistates*b))) = vf(hh*vfblocksize2 + j+nrgrid*(i+nistates*b));
        }
      }
    }
  }

  Rcpp::IntegerVector dims(4);
  dims(0) = nhhs;
  dims(1) = nrgrid;
  dims(2) = nistates;
  dims(3) = nbsmall;

  vfout1.attr("dim") = dims;

  ret["vfsmall"] = vfout1;
  //ret["llqsave"] = llqsave;
  //ret["llbrsave"] = llbrsave;

  //Rcpp::IntegerVector dim2(3);
  //dim2[0] = llinfo.nobs;
  //dim2[1] = 1+llinfo.capj;
  //dim2[2] = nsize;

  //vfll.attr("dim") = dim2;
  //utll.attr("dim") = dim2;

  //ret["vfll"] = vfll;
  //ret["utll"] = utll;

  ret["risave2"] = risave2;

  Rcpp::IntegerVector dim2(3);
  dim2[0] = nobsdata;
  dim2[1] = 3;
  dim2[2] = cfrepkeep;
  savechoice.attr("dim") = dim2;
  ret["savechoice"] = savechoice;
  ret["revenues"] = revenues;
  ret["quantities"] = quantities;
  ret["nunits"] = nunits;
  ret["saveinv"] = saveinv;



  //Rcpp::List ret;
  //ret["a"] = 0;
  
  return ret;

  END_RCPP

}


// calculate choice probabilities at a vector of prices and inventories (based on choicesimijc)
// this is for the identification plots using the data
// when we do this, we should pre-solve for the value function 
// and inclusive value process using ll or computevf
// here, hhinds will index which households we want to compute everything for

// [[Rcpp::export]]
Rcpp::List getcprob(Rcpp::NumericMatrix param, Rcpp::NumericMatrix prices, 
		    Rcpp::NumericVector inventories, Rcpp::List info,
		    Rcpp::List ivvars, Rcpp::NumericVector vfin, Rcpp::NumericVector oldparam,
		    Rcpp::NumericVector oldrgrid, int nrgrid, Rcpp::IntegerVector indexes,
		    Rcpp::IntegerVector hhinds, int gridmethod,
		    Rcpp::IntegerVector vfimethod, Rcpp::NumericMatrix bwmatrix) {

  BEGIN_RCPP

  int rep = 1;

  int nobs = prices.nrow();

  int nbrand = Rcpp::as< std::vector<int> >(info["nbrand"])[0];
  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];
  int invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  int maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];
  
  int nsize = Rcpp::as< std::vector<int> >(info["nsize"])[0];

  int inttype = Rcpp::as< std::vector<int> >(info["inttype"])[0];

  double * iv[nsize];

  SEXP obshhinds1(info["obshhinds"]);
  Rcpp::IntegerVector obshhinds(obshhinds1);

  SEXP obshhinds2(info["obshhindsmerge"]);
  Rcpp::IntegerVector obshhindsmerge(obshhinds2);

  int hmodel = Rcpp::as< std::vector<int> >(info["h.model"])[0];

  int maxpind = 0;

  if(hmodel) {
    maxpind = Rcpp::as< std::vector<int> >(info["hmaxpind"])[0];
  } else {
    maxpind = Rcpp::as< std::vector<int> >(info["maxpind"])[0];
  }

  SEXP pindsbig1(info["idpriceindbig"]);
  Rcpp::IntegerVector pindsbig(pindsbig1);

  SEXP hpindsbig1(info["hpriceindbig"]);
  Rcpp::IntegerVector hpindsbig(hpindsbig1);

  Rcpp::IntegerVector pinds(pindsbig.size());
  if(hmodel) {
    pinds = Rcpp::clone(hpindsbig);
  } else {
    pinds = Rcpp::clone(pindsbig);
  }

  int rgridbd = gridmethod == 2 ? nbrand : nsize;

  Rcpp::NumericMatrix rgrid(rgridbd,nrgrid);

  if (gridmethod == 2) {
    if(info.containsElementNamed("rgrid")) {
      SEXP temp1rgrid = info["rgrid"];
      Rcpp::NumericMatrix temp2rgrid(temp1rgrid);
      rgrid = Rcpp::clone(temp2rgrid);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$rgrid must be specified.");
    }
  }

  int vfinterpmethod = vfimethod(0);

  Rcpp::NumericVector impdist(nrgrid);

  if(gridmethod == 2) {
    if(info.containsElementNamed("impdist")) {
      SEXP temp1imp = info["impdist"];
      Rcpp::NumericVector temp2imp(temp1imp);
      impdist = Rcpp::clone(temp2imp);
    } else {
      throw std::runtime_error("If gridmethod is 2, info$impdist must be specified.");
    }
  }

  


  int nsave = 1;
  int ntilde = 1;


  SEXP cdrawtemp(info["cdraws"]);
  Rcpp::NumericVector cdraws(cdrawtemp);

  int usecdraws = Rcpp::as< std::vector<int> >(info["usecdraws"])[0];

  int capj = Rcpp::as< std::vector<int> >(info["bigJ"])[0];

  int nsim = Rcpp::as< std::vector<int> >(info["nsim"])[0];

  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];

  int ninitt = Rcpp::as< std::vector<int> >(info["ninitt"])[0];

  int ncutiv = Rcpp::as< std::vector<int> >(info["ncutiv"])[0];

  int retinv = Rcpp::as< std::vector<int> >(info["retinv"])[0];

  Rcpp::IntegerMatrix dropmat(nhhs,nsize);
  if(info.containsElementNamed("dropmat")) {
    SEXP dropmat1 = info["dropmat"];
    Rcpp::IntegerMatrix dropmat2(dropmat1);
    dropmat = clone(dropmat2);
  }

  Rcpp::IntegerVector ivdrop(nhhs);
  if(info.containsElementNamed("ivdrop")) {
    SEXP ivdrop1 = info["ivdrop"];
    Rcpp::IntegerVector ivdropb(ivdrop1);
    ivdrop = clone(ivdropb);
  }

  Rcpp::IntegerVector ivdrop2(nhhs);
  ivdrop2 = clone(ivdrop);


  SEXP ic1(info["initcdraws"]);
  Rcpp::NumericVector ic(ic1);
  
  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);
  
  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP rarray1(info["revarray"]);
  Rcpp::IntegerVector revarray(rarray1);

  SEXP badmissable1(info["badmissable"]);
  Rcpp::IntegerVector badmissable(badmissable1);

  SEXP bindex1(info["bindex"]);
  Rcpp::IntegerVector bindex(bindex1);

  
  double maxinv = Rcpp::as< std::vector<double> >(info["maxinv"])[0];

  int myopic = Rcpp::as< std::vector<int> >(info["myopic"])[0];

  int cmodel = Rcpp::as< std::vector<int> >(info["contmodel"])[0];

  int debug = Rcpp::as< std::vector<int> >(info["debug"])[0];

  int sptype = Rcpp::as< std::vector<int> >(info["splinetype"])[0];

  int idrawtype = Rcpp::as< std::vector<int> >(info["idrawtype"])[0];

  SEXP brsizetemp(info["brsize"]);
  Rcpp::IntegerVector brsize(brsizetemp);

  int first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];

  int hinvbound = Rcpp::as< std::vector<int> >(info["h.invbound"])[0];

  SEXP initx1(info["initx"]);
  Rcpp::NumericVector initx(initx1);
  int linitx = initx.length();

  if(debug && retinv)
    throw std::range_error("Cannot have both debug and retinv set to TRUE.");

  double ll = 0;

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  int siminitchoices = Rcpp::as< std::vector<int> >(info["siminitchoices"])[0];

  SEXP ivcoef1(ivvars[0]);
  Rcpp::NumericVector ivcoef(ivcoef1);
  SEXP cutmat1(ivvars[1]);
  Rcpp::NumericVector cutmat(cutmat1);
  SEXP ivvari1(ivvars[2]);
  Rcpp::NumericVector ivvari(ivvari1);
  SEXP detiv1(ivvars[3]);
  Rcpp::NumericVector detiv(detiv1);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);
  int nd = Rf_length(info["ngrid"]);
  int ncut = Rcpp::as< std::vector<int> >(info["ncut"])[0];

  
  SEXP gvari1(info["gvari"]);
  Rcpp::NumericVector gvari(gvari1);
  SEXP gmean1(info["gmean"]);
  Rcpp::NumericVector gmean(gmean1);
  double dg = Rcpp::as< std::vector<int> >(info["dg"])[0];

  int nistates = ngrid(nd-1);

  double rgridhh[(nd-1)*nrgrid];

  double ivtprobs[nrgrid*maxpind];
  int ivtpfilled[maxpind];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }

  int ivblocksize;
  if(sptype == 1) {
    if(ncutiv == 0) {
      ivblocksize = (2*ncutiv+2)*nsize;
    } else {
      ivblocksize = (2*ncutiv+2)*nsize*nsize;
    }
  } else {
    if(ncutiv == 0) {
      ivblocksize = (nsize + 1)*nsize;
    } else {
      ivblocksize = (ncutiv+4)*nsize*nsize;
    }
  }

  int cutblocksize = (ncutiv+1)*nsize;

  int nrep, nreptilde;

  if(rep <= nsave) {
    nrep = rep-1;
  } else {
    nrep = nsave;
  }

  if(nrep < ntilde) {
    nreptilde = rep-1;
  } else {
    nreptilde = ntilde;
  }

  double initinv[nhhs];
  int initb[nhhs];

  memset(initb,0,nhhs*sizeof(int));

  int nco = INTEGER(Rf_getAttrib(param, R_DimSymbol))[0];

  Rcpp::NumericVector llhh(nhhs);
  Rcpp::NumericVector llhhbr(nhhs);

  double ut[(1+capj)*nsize];
  int drop[(1+capj)*nsize];
  double vfhat[(1+capj)*nsize];
  double cerror[(1+capj)*nsize];
  int bchoice[(1+capj)*nsize];
  double ichoice[(1+capj)*nsize];

  int b=0;
  double ivout[nsize];
  double ivgrid[(1+capj)*nsize];
  double vftemp[nbsmall*nistates*maxpind];
  int vffilled[nbsmall*nistates*maxpind];
  for(int i=0;i<nbsmall*nistates*maxpind;i++) {
    vffilled[i] = -1;
  }
  for(int i=0;i<maxpind;i++) {
    ivtpfilled[i] = -1;
  }
  double invprime = 0;
  int coorder[nsave];

  int lencdraws = Rf_length(info["cdraws"]);

  if(lencdraws == 0 && usecdraws) {
    throw std::range_error("cdraws is length zero but usecdraws is true");
  }

  Rcpp::NumericVector rllrun(nsim*nobs);
  Rcpp::NumericVector risave1(nsim*nhhs*ninitt);
  Rcpp::List invsave(2);

  int hhprint = Rcpp::as< std::vector<int> >(info["hhprint"])[0]-1;
  int repprint = Rcpp::as< std::vector<int> >(info["repprint"])[0];

  int bprime = 0;

  // compute initial inventory
  // we might want to approach this a bit differently - probably should have an AR step??


  double clb, cub, alpha, gamma, nu, beta, ccost;
  double xi[nbrand];

  int lomega = 2;

  if(invmodel==2) {
    lomega = Rcpp::as< std::vector<int> >(info["lomega"])[0];
  }

  double omega[lomega];
  double invub = maxinv;

  double inv = 0;
  
  // simulate initial choices

  int iobs = -1;
  int iobsbig = 0;

  double kweights[nsave];

  int vfind=0;

  int nhhc = hhinds(1)-hhinds(0)+1;

  Rcpp::NumericVector cprobout(nhhc*nobs*(1+capj)*nsize);
  //Rcpp::IntegerVector ichoicej(initx.length());
  //Rcpp::IntegerVector ichoices(initx.length());
  //Rcpp::IntegerVector ibsave(initx.length());
  //Rcpp::IntegerVector inbr(initx.length());

  //Rcpp::Rcout << "Starting calculation" << std::endl;
  
  for(int n=hhinds(0)-1;n<hhinds(1);n++) {
    inv = 0;
    b = 0;
    for(int j=0;j<lomega;j++) {
      omega[j] = param(nbrand+9+j,n);
    }
    clb = param(nbrand,n);
    cub = param(nbrand+1,n);

    alpha = param(nbrand+2,n);
    gamma = param(nbrand+3,n);
    nu = param(nbrand+4,n);
    beta = param(nbrand+5,n);
    ccost = param(nbrand+7,n);

    vfind = n*nrgrid*nbsmall*nistates*nsave;

    if(gridmethod == 2) {
      for(int p=0;p<nrgrid;p++) {
        indcheck(nbrand*p,0,rgrid.length(),"rgrid offset bug p: %d",p);

        computeiv(param.begin()+nco*n,rgrid.begin()+nbrand*p,brsize.begin(),
                  nsize,nbrand,rgridhh+(nd-1)*p);

      }
    }
    //Rcpp::Rcout << "aaa" << std::endl;
    for(int t=0;t<nobs;t++) {

      double chosenq = 0;
      double rcnew;
      if(usecdraws) {
	rcnew = ic(n+nhhs*t);
      } else {
	rcnew = Rcpp::runif(1)[0];
      }
      if(rcnew < 0.5) {
	rcnew = clb + (cub-clb)*sqrt(0.5*rcnew);
      } else {
	rcnew = cub - (cub-clb)*sqrt(0.5*(1.0-rcnew));
      }
      // indcheck(t+ninitt*n,0,linitx,"initx: n: %d; t %d",n,t);
      //chosenq = initx(t+ninitt*n);  // since this vector is ordered by id, then week, this should be the right indexing
      // chosenq will get overwritten if we're simulating initial choices
    
      invub = hinvbound ? omega[0] : maxinv;
      
      double ivinit[nsize];

      inv = inventories(t);

      // If we are resimulating choices we should create new inclusive values

      int nprint = -1; //282;
      int tprint = -1; //1;
	
      int nbravail = 0;
      if(n == nprint && t <= tprint) {
	Rcpp::Rcout << "****** n: " << n << ", t: " << t << std::endl;
      }
      int drop1[nbrand];
	
      for(int size=0;size<nsize;size++) {
	double umax = -10000000;
	for(int j=0;j<nbrand;j++) {
	  double price = prices(t,j);
	  if(price >= 999 || brsize[j+nbrand] != size+1) {
	    drop1[j] = 1;
	  } else {
	    nbravail++;
	    ut[j] = param(j,n) + alpha*price;
	    drop1[j] = 0;
	    if(ut[j] > umax)
	      umax = ut[j];
	    if(n == nprint && t <= tprint) {
	      Rcpp::Rcout << " size: " << size << " ut[" << j << "]: " << ut[j] << " price: " << price << ". ";
	    }
	  }
	}
	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << std::endl;
	}
	double s = 0;
	for(int j=0;j<nbrand;j++) {
	  if(!drop1[j]) {
	    s += exp(ut[j]-umax);
	    if(n == nprint && t <= tprint) {
	      Rcpp::Rcout << "j: " << j << "; ut[j]: " << ut[j] << "; udiff" << ut[j] - umax << "; s: " << s << std::endl;
	    }
	  }
	}
	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << "s: " << s << "; umax: " << umax << std::endl;
	}
	s=log(s)+umax;
	ivinit[size] = s;
	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << "ivinit[" << size << "]: " << ivinit[size] << std::endl;
	}

	for(int j=0;j<1+capj;j++) {
	  ivgrid[j+size*(1+capj)] = ivinit[size]*((double)j)/((double)capj);
	}
      }

      //Note that the code below doesn't make use of price states.  It's not worth the trouble
      //since this part only gets run once

      if(!myopic) {
	if(vfinterpmethod == 1) {
	  getclosestind(param.begin() + nco*n, oldparam.begin() + nsave*nco*n, indexes.begin(), nco, nrep, nsave, rep, coorder);
	} else {
	  if(!(hmodel && n > hhinds[0]-1)) {
	    getkweights(param.begin() + nco*n, oldparam.begin() + nsave*nco*n,
			indexes.begin(), bwmatrix.begin(), nco, nrep, nsave,
			rep, coorder, kweights,0);
	  }
	}
      }

      if(!myopic && beta > 0) {
	ivpred(ivinit, ivcoef.begin() + ivblocksize*n, cutmat.begin() + cutblocksize*n, nd, ncut*(ivdrop(n)==0), sptype, inttype, capj, ivout);
	if(gridmethod == 2) {
	  for(int p=0;p<nrgrid;p++) {
	    ivtprobs[p] = impfn2(ivout,rgridhh+(nd-1)*p,nd-1,ivvari.begin()+(nd-1)*(nd-1)*n,
				 detiv(n),impdist(p));
	    if(n==nprint && t <= tprint && p < 10) {
	      
	      Rcpp::Rcout << "ivtprobs[" << p << "]: " << ivtprobs[p] << ", " << impdist(p) <<
		", " << detiv(n) << std::endl;
	      for(int ii=0;ii<nd-1;ii++) {
		Rcpp::Rcout << "rgridhh[" << ii << "]: " << *(ii+rgridhh+(nd-1)*p) << ", ";
	      }
	      Rcpp::Rcout << std::endl;
	    }
	  }
	}

	int printflag = 0;

	if(n == nprint && t <= tprint) {
	  Rcpp::Rcout << "vfin: " << vfin[0] << ", " << vfin[1] << ", " << vfin[2] << std::endl;
	}

	int nextbflag = 0;
	//Rcpp::Rcout << "n, t: " << n << ", " << t << std::endl;
	for(int size=0;size<nsize;size++) {
	  for(int j=0;j<capj+1;j++) {
	    if( (j == 0 && size == 0) || j > 0 ) {
	      int bprime1 = 0;
	      if(invmodel == 1) {
		invprime = dmax(inv + packsize(size)*((double)j) - rcnew,0);
		invprime = invprime > invub ? invub : invprime;
	      } else {
		nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
			rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
		bprime1 = bprime;
		bprime = bindex[bprime]-1;
		indcheck(bprime,0,nbsmall,"bprime (init sim) error");
	      }

	      int iind = itoiind(invprime,istates,nistates,first20);
	      //int vfindx1 = iind + nistates*(bprime+nb*(pinds[n]-1));
                

	      // uncomment if not using price states
	      int vfindx1 = iind + nistates*bprime;
	      int ok = iind < nistates; //check to make sure the state above this one is admissable (I think the second part should always be true)
	      indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug (choice sim - 1). iind-1 %d, bprime %d",iind-1,bprime);

	      if(invmodel == 2) {
		ok = ok && ((bprime == 0 && iind == 0) || (bprime > 0 && pgrid[nd-1+nd*iind] <= packsize[bstates[0+maxbottles*bprime1]-1]));
	      }
	      if(ok) {
		if(vffilled[vfindx1] < 0) {
		  vffilled[vfindx1] = 1;
		  if(gridmethod==2) {
		    vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
						  ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
						  gvari.begin(), detiv(n), dg, gmean.begin(),
						  nd-1, istates, nistates, first20, bprime, nbsmall, 0,0,vfinterpmethod);

		  } else {
		    vftemp[vfindx1] = vfavgijchh(ivout, iind, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
						 ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
						 gvari.begin(), detiv(n), dg, gmean.begin(),
						 nd-1, istates, nistates, first20, bprime, nbsmall, printflag);
		  }

		  if(n == nprint && t <= tprint) {
		    Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		  }


		  //vfout(n+obsend*iind) = vftemp[vfindx1];

		}
	      }

	      //vfindx1 = iind-1 + nistates*(bprime+nb*(pinds[n]-1));
	      vfindx1 = iind-1 + nistates*bprime;
	      indcheck(vfindx1,0,nbsmall*nistates*maxpind,"vfindx1 1 bug (choice sim - 2). iind-1 %d, bprime %d",iind-1,bprime);
	      if(vffilled[vfindx1] < 0) {
		vffilled[vfindx1] = 1;
		if(gridmethod == 2) {
		  vftemp[vfindx1] = vfavgijchh2(ivtprobs, iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder, kweights,
						ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
						gvari.begin(), detiv(n), dg, gmean.begin(),
						nd-1, istates, nistates, first20, bprime, nbsmall, 0, 0,vfinterpmethod);
		} else {
		  vftemp[vfindx1] = vfavgijchh(ivout, iind-1, vfin.begin()+vfind, oldrgrid.begin(), indexes.begin(), coorder,
					       ntilde, nsave, nrgrid, ivvari.begin()+(nd-1)*(nd-1)*n,
					       gvari.begin(), detiv(n), dg, gmean.begin(),
					       nd-1, istates, nistates, first20, bprime, nbsmall, 0);
		}

		if(n == nprint && t <= tprint) {
		  Rcpp::Rcout << "vftemp[" << vfindx1 << "]: " << vftemp[vfindx1] << std::endl;
		}

		//vfout(n+obsend*(iind-1)) = vftemp[vfindx1];

	      }

	    }
	  }
	}
      }

      int first = 1;

      int vfoffset = 0;
      double maxutil = -1000000;
      int nextbflag = 0;
        
      for(int size=0;size<nsize;size++) {
	for(int j=0;j<capj+1;j++) {
	  int indx = j+(1+capj)*size;
	  drop[indx] = 1;
	  if( (j == 0 && size == 0) || j > 0 ) {
	    drop[indx] = 0;
	    if(invmodel == 1) {
	      invprime = dmax(inv + packsize[size]*((double)j) - rcnew,0);
	      invprime = invprime > invub ? invub : invprime;
	      ut[indx] = utiliv(j, &b, ivgrid[indx], inv, packsize(size), rcnew, gamma, nu, ccost, omega, cmodel, 0);
	    } else {
	      nextinv(j,size,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		      rcnew,revarray.begin(),nsize,badmissable.begin(),&bprime,&invprime,&nextbflag,0);
	      bchoice[indx] = bprime;
	      indcheck(bprime,0,nb,"bprime (init sim) error");
	      ut[indx] = utiliv2(j,size,ivgrid[indx],rcnew, gamma, nu, ccost, omega, lomega, cmodel, inv, b, bprime,
				 nextbflag, bstates.begin(), nb, maxbottles, packsize.begin(), 0);
	      bprime = bindex[bprime]-1;		
	      ichoice[indx] = invprime;
	    }
              
	    if(!myopic && beta > 0) {
	      vfhat[indx] = beta*vfinterplinfast(invprime, vftemp+vfoffset, bprime, bchoice[indx], istates, nistates, first20,vffilled+vfoffset,invmodel,bstates.begin(),nb,maxbottles,packsize.begin());
	    } else {
	      vfhat[indx] = 0.0;
	    }
	    if(n == nprint && t <= tprint) {
	      Rcpp::Rcout << "ut[" << indx << "]: " << ut[indx] << ", vf: " << vfhat[indx] << ", ivgrid: " << ivgrid[indx] <<std::endl;
	    }
	    maxutil = dmax(maxutil,ut[indx]+vfhat[indx]);
	    
	  }
	}
      }

      double s1 = 0.0;

      int schoice=0;
      int jchoice=0;
      
      double choiceumax = -1000000;

      for(int size=0;size<nsize;size++) {
	for(int j=0;j<capj+1;j++) {
	  int indx = j+(1+capj)*size;
	  if(!drop[indx]) {
	    ut[indx] -= maxutil;
	    cerror[indx] = -log(-log(R::runif(0,1)));
	    s1 += exp(ut[indx]+vfhat[indx]);
	    cprobout(n+nhhc*(t+nobs*(j+(1+capj)*size))) = exp(ut[indx]+vfhat[indx]);
	    if(ut[indx]+vfhat[indx]+cerror[indx] > choiceumax) {
	      choiceumax = ut[indx]+vfhat[indx]+cerror[indx];
	      if(invmodel==2) {
		bprime = bchoice[indx];
		invprime = ichoice[indx];
	      }
	      if(j > 0) {
		schoice = size+1;
		jchoice = j;
	      }
	    }
	  }
	}
      }
      
      for(int size=0;size<nsize;size++) {
	for(int j=0;j<capj+1;j++) {
	  int indx = j+(1+capj)*size;
	  if(!drop[indx]) {
	    cprobout(n+nhhc*(t+nobs*(j+(1+capj)*size))) /= s1;
	  }
	}
      }

      if(n == nprint && t <= tprint) {
	Rcpp::Rcout << "DEBUG: b: " << b << ", inv : " << inv <<
	  ", bprime: " << bprime << ", invprime; " << invprime <<
	  ", jchoice: " << jchoice << ", schoice: " << schoice << std::endl;
	
	int tempb; double tempi;
	if(invmodel==2) {	  
	  nextinv(jchoice,schoice-1,b,inv,bstates.begin(),nb,maxbottles,packsize.begin(),
		  rcnew,revarray.begin(),nsize,badmissable.begin(),&tempb,&tempi,&nextbflag,0);
	  Rcpp::Rcout << "TEST: " << tempb << ", " << tempi << std::endl;
	}
	for(int size=0;size<nsize;size++) {
	  for(int j=0;j<capj+1;j++) {
	    int indx = j+(1+capj)*size;
	    if(!drop[indx]) {
	      Rcpp::Rcout << "choice: " << indx << " (" << j << ", " << size << "): utility: " << ut[indx]
			  << " vf: " << vfhat[indx] << " error: " << cerror[indx]
			  << " bchoice: " << bchoice[indx] << " ichoice: " << ichoice[indx]
			  << "choice prob: " << exp(ut[indx]+vfhat[indx])/s1 
			  << "tutil: " << ut[indx]+vfhat[indx]+cerror[indx] << std::endl;		
	    }
	  }
	}
	//if(t==tprint) {throw std::range_error("stopping");}
      }

      // reset vffilled
      memset(vffilled,-1,nb*nistates*maxpind*sizeof(int));

      //Rcpp::Rcout << "done t " << t << std::endl;

    }

  }

  Rcpp::List res;
  Rcpp::IntegerVector dim1(4);
  dim1[0] = nhhc;
  dim1[1] = nobs;
  dim1[2] = 1+capj;
  dim1[3] = nsize;

  cprobout.attr("dim") = dim1;
  res["cprob"] = cprobout;

  return res;

  END_RCPP

}

// vf interpolation wrapper

// [[Rcpp::export]]
Rcpp::List vfinterpw(Rcpp::NumericVector inventories, 
		     Rcpp::IntegerVector bottlestates, 
		     Rcpp::IntegerVector bbigstates, 
		     Rcpp::IntegerVector pstates,
		     Rcpp::List info, Rcpp::NumericVector vfin) {

  BEGIN_RCPP
    ;
  int nobs = inventories.length();

  int nb = Rcpp::as< std::vector<int> >(info["nb"])[0];
  int nbsmall = Rcpp::as< std::vector<int> >(info["nbsmall"])[0];

  int invmodel = Rcpp::as< std::vector<int> >(info["invmodel"])[0];
  int maxbottles = Rcpp::as< std::vector<int> >(info["maxbottles"])[0];
  int first20 = Rcpp::as< std::vector<int> >(info["first20"])[0];
  int nhhs = Rcpp::as< std::vector<int> >(info["nhhs"])[0];
  int nrgrid = Rcpp::as< std::vector<int> >(info["nrgrid"])[0];

  SEXP packsize1(info["packsize"]);
  Rcpp::NumericVector packsize(packsize1);
  
  SEXP bstates1(info["bstates"]);
  Rcpp::IntegerVector bstates(bstates1);

  SEXP pgrid1(info["pgrid"]);
  Rcpp::NumericVector pgrid(pgrid1);

  SEXP ngrid1(info["ngrid"]);
  Rcpp::IntegerVector ngrid(ngrid1);
  int nd = Rf_length(info["ngrid"]);

  int nistates = ngrid(nd-1);

  //double rgridhh[(nd-1)*nrgrid];

  //double ivtprobs[nrgrid*maxpind];
  //int ivtpfilled[maxpind];

  double istates[nistates];
  for(int i=0;i<nistates;i++) {
    istates[i] = pgrid[nd-1+nd*i];
  }

  Rcpp::NumericVector vfout(nobs*nhhs);

  int vffilled[nbsmall*nistates*nrgrid];
  memset(vffilled,1,nbsmall*nistates*nrgrid*sizeof(int));

  double vftemp[nbsmall*nistates];
  memset(vftemp,0,nbsmall*nistates*sizeof(double));
  for(int hh=0;hh<nhhs;hh++) {
    for(int i=0;i<nobs;i++) {

      int p = pstates(i)-1;

      for(int ii=0;ii<nistates;ii++) {
	for(int b=0;b<nbsmall;b++) {
	  vftemp[ii+nistates*b] = vfin(hh+nhhs*(p+nrgrid*(ii+nistates*b)));
	}
      }

      vfout(i+nobs*hh) = vfinterplinfast(inventories(i), vftemp, bottlestates(i)-1, 
					 bbigstates(i)-1, istates,  nistates, first20, 
					 vffilled,  invmodel, bstates.begin(), nb,
					 maxbottles, packsize.begin());
    }
  }

  Rcpp::List ret;

  Rcpp::IntegerVector dims(2);
  dims(0) = nobs;
  dims(1) = nhhs;
  //dims[2] = nrgrid;
  vfout.attr("dim") = dims;

  ret["vf"] = vfout;

  return ret;

  END_RCPP


}

