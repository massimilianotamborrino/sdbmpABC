#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

//#General Matrix-Vector Multiplication mv_mult
// [[Rcpp::export]]
NumericVector mv_multgen_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  double temp=0;
  for(int i = 0; i < mat.nrow(); i++)
  {
    for(int j = 0; j < vec.size(); j++)
    {
      temp = temp + mat(i,j) * vec[j];
    }
    ret[i]=temp;
    temp=0;
  }
  return ret;
};


//#Matrix-Vector Multiplication when the matrix is cm, a sparse matrix
//#that has non zero in the main diagonal and in cm(4,1), cm(5,2), cm(6,3)
// [[Rcpp::export]]
NumericVector mv_multcm_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  ret[0]=mat(0,0)*vec[0];
  ret[1]=mat(1,1)*vec[1];
  ret[2]=mat(2,2)*vec[2];
  ret[3]=mat(3,3)*vec[3]+mat(3,0)*vec[0];
  ret[4]=mat(4,4)*vec[4]+mat(4,1)*vec[1];
  ret[5]=mat(5,5)*vec[5]+mat(5,2)*vec[2];
  return ret;
};


//#Matrix-Vector Multiplication mv_mult when the matrix is dm,
//#a band matrix that has non zero in the main diagonal and in the
//#3rd diagonal above
// [[Rcpp::export]]
NumericVector mv_multdm_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  ret[0]=mat(0,0)*vec[0]+mat(0,3)*vec[3];
  ret[1]=mat(1,1)*vec[1]+mat(1,4)*vec[4];
  ret[2]=mat(2,2)*vec[2]+mat(2,5)*vec[5];
  ret[3]=mat(3,3)*vec[3]+mat(3,0)*vec[0];
  ret[4]=mat(4,4)*vec[4]+mat(4,1)*vec[1];
  ret[5]=mat(5,5)*vec[5]+mat(5,2)*vec[2];
  return ret;
};


//sigmoid function
// [[Rcpp::export]]
double sigmoid_Cpp_(double x, double vmax, double v0, double r)
{
  double ret=vmax/(1+exp(r*(v0-x)));
  return ret;
};

//Solution of the linear SDE subsystem for general matrices dm and cm
// [[Rcpp::export]]
NumericVector SDE_Cpp_gen_(NumericVector vec, NumericMatrix dm, NumericMatrix cm, NumericVector randvec)
{
  NumericVector ret=mv_multgen_(dm,vec)+mv_multgen_(cm,randvec);
  return ret;
};


//Solution of the linear SDE subsystem when dm and cm are the specifc matrices for the JR-NMM
// [[Rcpp::export]]
NumericVector SDE_Cpp_(NumericVector vec, NumericMatrix dm, NumericMatrix cm, NumericVector randvec)
{
  NumericVector ret=mv_multdm_(dm,vec)+mv_multcm_(cm,randvec);
  return ret;
};

//Solution of the non-linear ODE subsystem
// [[Rcpp::export]]
NumericVector ODE_Cpp_(NumericVector vec, double h, double Aa, double mu, double BbC, double C1, double C2, double C3, double vmax, double v0, double r)
{
  NumericVector help(vec.size());
  help(3)=Aa*sigmoid_Cpp_(vec(1)-vec(2),vmax,v0,r);
  help(4)=Aa*(mu+C2*sigmoid_Cpp_(C1*vec(0),vmax,v0,r));
  help(5)=BbC*sigmoid_Cpp_(C3*vec(0),vmax,v0,r);
  NumericVector ret=vec+h*help;
  return ret;
};

//#Splitting for the JR-NMM
// [[Rcpp::export]]
NumericMatrix Splitting_JRNMM_gen_Cpp_(double h_i,NumericVector startv, NumericVector grid_i, NumericMatrix dm_i, NumericMatrix cm_i, double mu_i, double C_i, double A, double B, double a, double b,double v0, double r, double vmax)
{
  double h=h_i;
  NumericVector start=startv;
  NumericVector grid=grid_i;
  int iter=grid.size();

  NumericMatrix randarr(6,iter);
  randarr(0,_)=rnorm(iter);
  randarr(1,_)=rnorm(iter);
  randarr(2,_)=rnorm(iter);
  randarr(3,_)=rnorm(iter);
  randarr(4,_)=rnorm(iter);
  randarr(5,_)=rnorm(iter);

  NumericMatrix dm=dm_i;
  NumericMatrix cm=cm_i;

  double mu=mu_i;
  double C1=C_i;

  //fixed values: from the original JR-NMM

  double C2=0.8*C1;
  double C3=0.25*C1;
  double C4=C3;
  double Aa=A*a;
  double BbC=B*b*C4;

  NumericMatrix sol(6,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec(6);

  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    newv=SDE_Cpp_gen_(newv,dm,cm,randvec);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    sol(_,i)=newv;
  }

  NumericMatrix ret=sol;

  return sol;
};


//#Splitting for the JR-NMM
// [[Rcpp::export]]
NumericMatrix Splitting_JRNMM_Cpp_(double h_i,NumericVector startv, NumericVector grid_i, NumericMatrix dm_i, NumericMatrix cm_i, double mu_i, double C_i, double A, double B, double a, double b,double v0, double r, double vmax)
{
  double h=h_i;
  NumericVector start=startv;
  NumericVector grid=grid_i;
  int iter=grid.size();

  NumericMatrix randarr(6,iter);
  randarr(0,_)=rnorm(iter);
  randarr(1,_)=rnorm(iter);
  randarr(2,_)=rnorm(iter);
  randarr(3,_)=rnorm(iter);
  randarr(4,_)=rnorm(iter);
  randarr(5,_)=rnorm(iter);

  NumericMatrix dm=dm_i;
  NumericMatrix cm=cm_i;

  double mu=mu_i;
  double C1=C_i;

  //fixed values: from the original JR-NMM
  double C2=0.8*C1;
  double C3=0.25*C1;
  double C4=C3;
  double Aa=A*a;
  double BbC=B*b*C4;

  NumericMatrix sol(6,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec(6);

  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    newv=SDE_Cpp_(newv,dm,cm,randvec);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    sol(_,i)=newv;
  }

  NumericMatrix ret=sol;

  return sol;
};

//#Splitting for the JR-NMM returning only the output y=x[2]-x[3]
// [[Rcpp::export]]
NumericVector Splitting_JRNMM_output_Cpp_(double h_i,NumericVector startv, NumericVector grid_i, NumericMatrix dm_i, NumericMatrix cm_i, double mu_i, double C_i, double A, double B, double a, double b,double v0, double r, double vmax)
{
  double h=h_i;
  NumericVector start=startv;
  NumericVector grid=grid_i;
  int iter=grid.size();

  NumericMatrix randarr(6,iter);
  randarr(0,_)=rnorm(iter);
  randarr(1,_)=rnorm(iter);
  randarr(2,_)=rnorm(iter);
  randarr(3,_)=rnorm(iter);
  randarr(4,_)=rnorm(iter);
  randarr(5,_)=rnorm(iter);

  NumericMatrix dm=dm_i;
  NumericMatrix cm=cm_i;

  double mu=mu_i;
  double C1=C_i;

  //fixed values: from the original JR-NMM
   double C2=0.8*C1;
  double C3=0.25*C1;
  double C4=C3;
  double Aa=A*a;
  double BbC=B*b*C4;

  NumericMatrix sol(6,iter);
  sol(_, 0)=start;
  NumericVector newv=start;
  NumericVector randvec(6);

  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    newv=SDE_Cpp_(newv,dm,cm,randvec);
    newv=ODE_Cpp_(newv,h/2,Aa,mu,BbC,C1,C2,C3,vmax,v0,r);
    sol(_,i)=newv;
  }
  NumericVector ret=sol(1,_)-sol(2,_);
  return ret;
};
