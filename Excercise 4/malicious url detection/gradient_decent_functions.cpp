// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
using namespace Rcpp;
using namespace std;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef MSpMat::InnerIterator InIterMat;
typedef NumericVector NV;
typedef IntegerVector IV;
typedef double db;
typedef vector<int> vi;
#define lop(i,a,b) for (int i=a; i<=b; i++)

//[[Rcpp::export]]
NV sgdC(IV rn, MSpMat X,NV y,db eps,int ite,db lambda){
  int p=X.cols(),r;
  NV beta(p),H(p,0.001),re(ite+p);
  db og,g;
  lop(i,0,ite-1){
    og=0;
    r=rn[i];
    lop(j,0,p-1)og-=X.coeff(r,j)*beta[j];
    og=1/(exp(og)+1);
    if(y[r])re[i]=-log(og);
    else re[i]=-log(1-og);
    lop(j,0,p-1){
      g=(og-y[r])*X.coeff(r,j);
      if(beta[j]>0)g+=lambda;
      else if(beta[j]<0)g-=lambda;
      H[j]+=g*g;
      beta[j]-=eps*g/sqrt(H[j]);
    }
  }
  lop(i,ite,ite+p-1)re[i]=beta[i-ite]; 
  return re;
}
//[[Rcpp::export]]
NV sgdC_sparse(NV rn, MSpMat X,NV y,db eps,int epoch,db lambda){
  int p=X.rows(),r,ite=X.cols()*epoch;
  NV H(p,0.001),re(ite+p),beta(p);
  db og,g;
  lop(i,0,ite-1){
    og=0;
    r=rn[i];
    for(InIterMat it(X,r); it;++it)og-=it.value()*beta[it.row()];
    og=1/(exp(og)+1);
    if(y[r])re[i]=-log(og);
    else re[i]=-log(1-og);
    /*vi updated(p,0);*/
    for(InIterMat it(X,r); it;++it){
      g=(og-y[r])*it.value();
      int j=it.row();
      if(beta[j]>0)g+=lambda;
      else if(beta[j]<0)g-=lambda;
      H[j]+=g*g;
      beta[j]-=eps*g/sqrt(H[j]);
      /*updated[it.row()]=1;*/
    }
    /*lop(j,0,p-1){
      if(!updated[j]){
        g=0;
        if(beta[j]>0)g=lambda;
        else if(beta[j]<0)g=-lambda;
        H[j]+=g*g;
        beta[j]-=eps*g/sqrt(H[j]);
      }
    }*/
  }
  lop(i,ite,ite+p-1)re[i]=beta[i-ite];       
  return re;
}