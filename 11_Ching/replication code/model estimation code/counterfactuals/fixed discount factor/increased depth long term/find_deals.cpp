// this file can be used to identify deals in the IRI store data
// it is based on code I wrote to do the same thing in the Nielsen
// homescan data

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

int quarter(const int month) {
  if(month >=1 && month <=3) {
    return 1;
  } else if (month >=4 && month <=6) {
    return 2;
  } else if (month >=7 && month <=9) {
    return 3;
  } else if (month >=10 && month <=12) {
    return 4;
  } else {
    throw std::range_error("Invalid month value in quarter calculation");
  }
}

// [[Rcpp::export]]
Rcpp::List findmodes(Rcpp::List stores, int nbrand, Rcpp::IntegerVector minmode) {

  BEGIN_RCPP

  SEXP dates1 = stores["WEEK"];
  Rcpp::IntegerVector dates(dates1);

  SEXP month1 = stores["month"];
  Rcpp::IntegerVector month(month1);

  SEXP stcode1 = stores["IRI_KEY"];
  Rcpp::IntegerVector stcode(stcode1);

  //SEXP upc1 = stores["upc"];
  //Rcpp::CharacterVector upc(upc1);
  
  int nobs = dates.size();
  
  Rcpp::NumericMatrix modalprices(nobs,nbrand);
  
  //SEXP prmult1 = stores["prmult"];
  //Rcpp::IntegerVector prmult(prmult1);

  //SEXP units1 = stores["units"];
  //Rcpp::IntegerVector units(units1);

  int minm = minmode(0);

  for(int u=0;u<nbrand;u++) {

    SEXP pvec1 = stores[2+u];
    Rcpp::NumericVector pvec(pvec1);

    int lprice = 100;

    Rcpp::NumericVector prices(lprice);
    Rcpp::IntegerVector pcount(lprice);

    //std::string currentupc = Rcpp::as< std::string>(upc(0));
    int currentstore = stcode(0);
  
    //Rcpp::Rcout << "date[0]: " << dates[0] << std::endl;
    //Rcpp::Rcout << "upc[0]: " << currentupc << std::endl;

    int restart = 1;

    int istart=0;

    int modeprice=-1; //index of maximum mode
    int nmode=-1; //number of obs at current mode

    for(int i=0;i<nobs;i++) {

      if(restart) {

	// new upc or store
	memset(prices.begin(),0,lprice*sizeof(double));
	memset(pcount.begin(),0,lprice*sizeof(int));

	currentstore = stcode(i);

	restart = 0;

	istart = i;

      }

      if(pvec(i) > 0) {

	int indx = -1;
	for(int j=0;j<lprice;j++) {
	  if(pvec(i) == prices(j) || pcount(j) == 0) {
	    indx = j;
	    break;
	  }
	}
	if(indx == -1) {
	  Rcpp::Rcout << "obs " << i+1 << std::endl;
	  for(int j=0;j<lprice;j++) {
	    Rcpp::Rcout << "prices(j): " << prices(j) << "; pcount(j): " << pcount(j) << std::endl;
	  }
	  throw std::range_error("Need to increase the size of lprice.");
	} else {
	  if(pcount(indx) == 0) {
	    prices(indx) = pvec(i);
	    pcount(indx) = 1;
	  } else {
	    pcount(indx)++;
	  }
	  
	}

      }

      if(i < nobs-1) {
	int nextstore = stcode(i+1);

	restart = (currentstore != nextstore) || (quarter(month(i)) != quarter(month(i+1)));

      }

      //if(i==166) {
      //  for(int j=0;j<lprice;j++) {
      //Rcpp::Rcout << "prices(j): " << prices(j) << "; pcount(j): " << pcount(j) << std::endl;
      //  }
      //  Rcpp::Rcout << "modal price (" << modeprice << "): " << prices(modeprice) << std::endl;
      //  Rcpp::Rcout << "restart " << i << ": " << restart << std::endl;
      //  Rcpp::Rcout << 
      //}

      if(restart || i == nobs-1) {


	nmode = pcount(0);
	modeprice = 0;
	for(int j=1;j<lprice;j++) {
	  if(pcount(j) == 0) {
	    break;
	  }
	  if(pcount(j) > nmode) {
	    nmode = pcount(j);
	    modeprice = j;
	  } else if (pcount(j) == nmode) {
	    if(minm) {
	      if(prices(j) < prices(modeprice)) {
		nmode = pcount(j);
		modeprice = j;
	      }
	    } else {
	      if(prices(j) > prices(modeprice)) {
		nmode = pcount(j);
		modeprice = j;
	      }
	    }
	  }
	}

	/*if(indx == modeprice) {
	  nmode = pcount(indx);
	  } else {

	  if(pcount(indx) > nmode) {
	  nmode = pcount(indx);
	  modeprice = indx;
	  } else if (pcount(indx) == nmode) {
	  if(minm) {
	  if(prices(indx) < prices(modeprice)) {
	  nmode = pcount(indx);
	  modeprice = indx;
	  }
	  } else {
	  if(prices(indx) > prices(modeprice)) {
	  nmode = pcount(indx);
	  modeprice = indx;
	  }
	  }
	  }
	  }*/
	
	for(int j=istart;j<=i;j++) {
	  modalprices(j,u) = prices(modeprice);
	}
      }

    }

  }

  Rcpp::List ret;

  ret["modalprice"] = modalprices;

  return ret;
    
  END_RCPP

}
