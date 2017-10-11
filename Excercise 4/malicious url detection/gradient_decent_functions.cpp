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
typedef vector<double> vdb;
#define sz(a) int((a).size()) 
#define lop(i,a,b) for (int i=a; i<=b; i++)
#define vlop(i,v) lop(i,0,sz(v)-1)
#define pb push_back

//[[Rcpp::export]]
NV sgdC(NV rn, MSpMat X,NV y,db eps,int epoch,db lambda){
  int p=X.rows(),r,ite=X.cols()*epoch;
  NV re(ite+p);
  vdb H(p,0.001),beta(p,0.001);
  vi nzbeta;
  db og,g;
  lop(i,0,ite-1){
    og=0;
    r=rn[i];
    for(InIterMat it(X,r); it;++it)og-=it.value()*beta[it.row()];
    og=1/(exp(og)+1);
    if(y[r])re[i]=-log(og);
    else re[i]=-log(1-og);
    for(InIterMat it(X,r); it;++it){
      g=(og-y[r])*it.value();
      int j=it.row();
      if(!beta[j])nzbeta.pb(j);
      H[j]+=g*g;
      beta[j]-=eps*g/sqrt(H[j]);
    }
    vlop(j,nzbeta){
      if(beta[nzbeta[j]]>0)beta[nzbeta[j]]-=lambda;
      else beta[nzbeta[j]]+=lambda;
    }
  }
  lop(i,ite,ite+p-1)re[i]=beta[i-ite];       
  return re;
}
