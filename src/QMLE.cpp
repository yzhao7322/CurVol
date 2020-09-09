
//library(RcppArmadillo)
//Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
//#include <Rcpp.h>
//#include <iostream.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @useDynLib CurVol
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double compute_sigma_inp_k_pq_cpp(unsigned int k, arma::mat d, arma::mat As,
                                  arma::mat Bs, arma::mat prev_ys,
                                  arma::mat prev_sigmas,unsigned int p,
                                  unsigned int q,int M){

  double sigma;
  arma::mat A;
  arma::mat B, tmp;

  sigma=d(k,0);

  for(unsigned int i=0;i<q;i++){
    // Rcout << "i is " <<i << std::endl;

    A=As.submat(k+i*M,0,k+i*M,M-1);
    // Rcout << "A is " <<A << std::endl;
    // Rcout << "tmp is " <<tmp << std::endl;
    tmp=(A*prev_ys.submat(0,i,M-1,i));
    sigma=sigma+tmp(0,0);
  }

  for(unsigned int i=0;i<p;i++){
    // Rcout << "i is " <<i << std::endl;
    B=Bs.submat(k+i*M,0,k+i*M,M-1);
    // Rcout << "tmp is  " <<tmp << std::endl;
    tmp=B*prev_sigmas.submat(0,i,M-1,i);
    sigma=sigma+tmp(0,0);
  }

  return sigma;
}


//' @useDynLib CurVol
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double obj_fn_pq_cpp(arma::mat d, arma::mat As,
                     arma::mat Bs, arma::mat ysq,
                     unsigned int p,
                     unsigned int q){

  //initialize
  arma::mat sigma0,y0,prev_sigma,prev_ys,new_sigma, sigmat;
  unsigned int cap=q;
  double obj;

  //ysq=NxM
  unsigned int N=size(ysq)(0);
  unsigned int M=size(ysq)(1);

  // Rcout << "ysq " <<size(ysq) << std::endl;


  sigma0=sum(ysq.submat(0,0,4,M-1),0)/5;
  y0=sigma0;


  //compute initial values to start recursion
  //

  prev_sigma.set_size(p,M);
  prev_ys.set_size(q,M);

  //
  // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
  // Rcout << "sigma0 is " <<sigma0 << std::endl;

  for(unsigned int i=0;i<p;i++)
    prev_sigma.row(i)=sigma0;

  for(unsigned int i=0;i<q;i++)
    prev_ys.row(i)=sigma0;

  // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
  // Rcout << "prev_ys is " <<prev_ys  << std::endl;


  //get initial max(p,q) sigma inner products:

  new_sigma.set_size(1,M);


  if(p>q)
    cap=p;

  sigmat.set_size(cap,M);

  for(unsigned int i=0;i<cap;i++){


    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_pq_cpp(j,d,As,Bs,prev_ys.t(),prev_sigma.t(),p,q,M);
    }
    sigmat.row(i)=new_sigma;

    // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
    // Rcout << "prev_ys is " <<prev_ys  << std::endl;
    // Rcout << "sigmat is " <<sigmat << std::endl;
    //
    if(q>1)
      prev_ys.submat(0,0,q-2,M-1)=prev_ys.submat(1,0,q-1,M-1);
    prev_ys.row(q-1)=ysq.row(i);

    if(p>1)
      prev_sigma.submat(0,0,p-2,M-1)=prev_sigma.submat(1,0,p-1,M-1);
    prev_sigma.row(p-1)=new_sigma;

  }

  // Rcout << "sigmat is " <<sigmat << std::endl;


  for(unsigned int i=cap;i<N;i++){
    //just do this matrix type in case p=1 or q=1
    prev_ys=ysq.submat(i-q,0,(i-1),M-1);
    prev_sigma=sigmat.submat(i-p,0,(i-1),M-1);
    //recursively compute sigma inner products

    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_pq_cpp(j,d,As,Bs,prev_ys.t(),prev_sigma.t(),p,q,M);
    }

    sigmat=join_vert(sigmat,new_sigma);
  }


  // Rcout << "sigmat is " <<sigmat << std::endl;

  obj=mean(sum(ysq/sigmat+log(sigmat),1));


  return obj;
}





/*** R

# As=matrix(1:(25*3),ncol=5,nrow=15)
# Bs=matrix((25*3+1):(25*2),ncol=5,nrow=10)
# # d=matrix(1:5,ncol=1)
# y_vec=matrix(rnorm(15),ncol=3,nrow=5)^2
# y_s=matrix(rnorm(10),ncol=2,nrow=5)^2
#
# d1=matrix(pars$d,ncol=1)
# As2=rbind(pars$As[[1]],pars$As[[2]],pars$As[[3]])
# Bs2=rbind(pars$Bs[[1]],pars$Bs[[2]])
#
# compute_sigma_inp_k_pq_cpp(k=0,d=d1,As=As2,Bs=Bs2,y_vec,y_s,p=2,q=3,M=5)
# # compute_sigma_inp_k_pq(1,pars$d,pars$As,pars$Bs,t(y_vec),t(y_s),2,3)
# obj_fn_pq_cpp(d=d1,As=As2,Bs=Bs2,ysq,p=2,q=3)
#
# obj_fn_pq(pars$d,pars$As,pars$Bs,ysq,2,3)

# pars

# obj_fn_pq_cpp(d,As,Bs,ysq,p,q)

# obj_fn_pq_cpp(d=matrix(pars$d,ncol=1),As=pars$As[[1]],Bs=pars$Bs[[1]],y_vec,p=1,q=1)

# compute_sigma_inp_k_pq(1,d,As,Bs,y_vec,y_s,2,3,5)

*/

