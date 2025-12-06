#include <Rcpp.h>

// find first el of array where indexes match

int findinarray(int i, int t, int s, const Rcpp::IntegerVector &i1, const Rcpp::IntegerVector &i2, const Rcpp::IntegerVector &i3, int li) {

  int done = 0;
  int index = 0;
  int findex = -1;

  while(!done) {
    done = index == li-1 || (i == i1(index) && t == i2(index) && s == i3(index));
    if(i == i1(index) && t == i2(index) && s == i3(index))
      findex = index;
    index++;
  }

  return(findex);

}

// code to do fast aggregation of prices

// [[Rcpp::export]]
Rcpp::List agprices(Rcpp::List mpricest, Rcpp::List mpricehh, Rcpp::List mfeat, Rcpp::List mdisp, Rcpp::List mpr, SEXP nbrsize, SEXP nweek, SEXP nstore) {

  int nbr = *INTEGER(nbrsize);
  int nw = *INTEGER(nweek);
  int ns = *INTEGER(nstore);


  Rcpp::NumericVector prarray(nbr*nw*ns);
  Rcpp::NumericVector prchange(nbr*nw*ns);
  Rcpp::IntegerVector farray(nbr*nw*ns);
  Rcpp::IntegerVector darray(nbr*nw*ns);
  Rcpp::IntegerVector ptarray(nbr*nw*ns);

  SEXP brindex1(mpricest["brindex"]);
  Rcpp::IntegerVector brindex(brindex1);

  SEXP stindex1(mpricest["stindex"]);
  Rcpp::IntegerVector stindex(stindex1);

  SEXP wkindex1(mpricest["wkindex"]);
  Rcpp::IntegerVector wkindex(wkindex1);

  SEXP ppunit1(mpricest["ppunit"]);
  Rcpp::NumericVector ppunitst(ppunit1);

  SEXP feat1(mfeat["Feat"]);
  Rcpp::IntegerVector feat(feat1);

  SEXP disp1(mdisp["Disp"]);
  Rcpp::IntegerVector disp(disp1);

  SEXP pr1(mpr["PR"]);
  Rcpp::IntegerVector pr(pr1);

  SEXP brindexhh1(mpricehh["brindex"]);
  Rcpp::IntegerVector brindexhh(brindexhh1);

  SEXP stindexhh1(mpricehh["stindex"]);
  Rcpp::IntegerVector stindexhh(stindexhh1);

  SEXP wkindexhh1(mpricehh["wkindex"]);
  Rcpp::IntegerVector wkindexhh(wkindexhh1);

  SEXP ppunithh1(mpricehh["ppunit"]);
  Rcpp::NumericVector ppunithh(ppunithh1);

  int lst = brindex.size();
  int lhh = brindexhh.size();

  for(int i=0;i<nbr;i++) {
    for(int t=0;t<nw;t++) {
      for(int s=0;s<ns;s++) {
        prarray[t+nw*(s+ns*i)] = -1;
        prchange[t+nw*(s+ns*i)] = -1;
        farray[t+nw*(s+ns*i)] = -1;
        darray[t+nw*(s+ns*i)] = -1;
        ptarray[t+nw*(s+ns*i)] = -1;
        int ind1 = findinarray(i+1,t+1,s+1,brindex,wkindex,stindex,lst);
        if(ind1 > -1) {
          prarray[t+nw*(s+ns*i)] = ppunitst[ind1];
          farray[t+nw*(s+ns*i)] = feat[ind1];
          darray[t+nw*(s+ns*i)] = disp[ind1];
          ptarray[t+nw*(s+ns*i)] = pr[ind1];
        }
        ind1 = findinarray(i+1,t+1,s+1,brindexhh,wkindexhh,stindexhh,lhh);
        if(ind1 > -1) {
          prchange[t+nw*(s+ns*i)] = prarray[t+nw*(s+ns*i)]-ppunithh[ind1];
          prarray[t+nw*(s+ns*i)] = ppunithh[ind1];
        }
      }
    }
  }

  Rcpp::List res(5);
  Rcpp::IntegerVector dim(3);

  dim(0) = nw;
  dim(1) = ns;
  dim(2) = nbr;

  prarray.attr("dim") = dim;
  prchange.attr("dim") = dim;
  farray.attr("dim") = dim;
  darray.attr("dim") = dim;
  ptarray.attr("dim") = dim;

  res(0) = prarray;
  res(1) = prchange;
  res(2) = farray;
  res(3) = darray;
  res(4) = ptarray;

  return(res);

}


















