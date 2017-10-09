library(Rcpp)
library(RcppEigen)
library(inline)

cppFunction(depends="RcppEigen",'
  NumericVector sgdC(IntegerVector rn, Eigen::MappedSparseMatrix<double> X,NumericVector y,NumericVector beta0,double eps,int ite,double lambda){
int p=beta0.size(),r;
NumericVector beta(p);
NumericVector H(p,1.0);
NumericVector re(ite+p);
double og,g;
for(int i=0;i<ite;i++){
og=0;
r=rn[i];
for(int j=0; j<p;j++)og-=X.coeff(r,j)*beta[j];
og=1/(exp(og)+1);
if(y[r])re[i]=-log(og);
else re[i]=-log(1-og);
for(int j=0;j<p;j++){
g=(og-y[r])*X.coeff(r,j);
if(beta[j]>0)g+=lambda;
else if(beta[j]<0)g-=lambda;
H[j]+=g*g;
beta[j]-=eps*g/sqrt(H[j]);
}
return beta;
}
for(int i=ite;i<ite+p;i++)re[i]=beta[i-ite];       
return re;
  }
')
cppFunction(depends="RcppEigen",'
  NumericVector sgdC_sparse(NumericVector rn, Eigen::MappedSparseMatrix<double> X,NumericVector y,NumericVector beta0,double eps,int ite,double lambda){
            NumericVector beta=beta0;
            int p=beta0.size();
            NumericVector H(p,1.0);
            NumericVector re(ite+p,1.0);
            double og,g;
            for(int i=0;i<ite;i++){
            og=0;
            int r=rn[i];
            for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,r); it;++it)og-=it.value()*beta[it.row()];
            og=1/(exp(og)+1);
            if(y[r])re[i]=-log(og);
            else re[i]=-log(1-og);
NumericVector updated(p,0.0);
            for(Eigen::MappedSparseMatrix<double>::InnerIterator it(X,r); it;++it){
               g=(og-y[r])*it.value();
int j=it.row();
if(beta[j]>0)g+=lambda;
            else g-=lambda;
            H[j]+=g*g;
            beta[j]-=eps*g/sqrt(H[j]);
updated[it.row()]=1;
            }
            for(int j=0;j<p;j++){
if(!updated[j]){
            if(beta[j]>0)g+=lambda;
            else g-=lambda;
            H[j]+=g*g;
            beta[j]-=eps*g/sqrt(H[j]);
}
            }
            }
            for(int i=ite;i<ite+p;i++)re[i]=beta[i-ite];       
            return re;
            }
            ')