library(Rcpp)
library(RcppEigen)
library(inline)

cppFunction(depends="RcppEigen",'
  NumericVector sgdC(NumericVector rn, Eigen::MappedSparseMatrix<double> X,NumericVector y,NumericVector beta0,double eps,int ite,double lambda){
NumericVector beta=beta0;
int p=beta0.size();
NumericVector H(p,1.0);
NumericVector re(ite+p,1.0);
double og,g;
for(int i=0;i<ite;i++){
og=0;
int r=rn[i];
for(int j=0; j<p;j++)og-=X.coeff(r,j)*beta[j];
og=1/(exp(og)+1);
if(y[r])re[i]=-log10(og);
else re[i]=-log10(1-og);
for(int j=0;j<p;j++){
g=(og-y[r])*X.coeff(r,j);
if(beta[j]>0)g+=lambda;
else g-=lambda;
H[j]+=g*g;
beta[j]-=eps*g/H[j];
}
}
for(int i=ite;i<ite+p;i++)re[i]=beta[i-ite];       
return re;
  }
')
