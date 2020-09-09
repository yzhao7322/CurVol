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
double compute_sigma_inp_k_x_cpp(unsigned int k, arma::mat d, arma::mat As,
                                 arma::mat Bs, arma::mat Gs, arma::mat prev_ys,
                                 arma::mat prev_sigmas, arma::mat prev_x, int M){

  double sigma;
  arma::mat A;
  arma::mat B;
  arma::mat G, tmp;

  sigma=d(k,0);

  for(unsigned int i=0;i<1;i++){
    // Rcout << "i is " <<i << std::endl;

    A=As.submat(k+i*M,0,k+i*M,M-1);
    // Rcout << "A is " <<A << std::endl;
    // Rcout << "tmp is " <<tmp << std::endl;
    tmp=(A*prev_ys.submat(0,i,M-1,i));
    sigma=sigma+tmp(0,0);
  }

  for(unsigned int i=0;i<1;i++){
    // Rcout << "i is " <<i << std::endl;
    B=Bs.submat(k+i*M,0,k+i*M,M-1);
    // Rcout << "tmp is  " <<tmp << std::endl;
    tmp=B*prev_sigmas.submat(0,i,M-1,i);
    sigma=sigma+tmp(0,0);
  }

  for(unsigned int i=0;i<1;i++){
    // Rcout << "i is " <<i << std::endl;
    G=Gs.submat(k+i*M,0,k+i*M,M-1);
    // Rcout << "tmp is  " <<tmp << std::endl;
    tmp=G*prev_x.submat(0,i,M-1,i);
    sigma=sigma+tmp(0,0);
  }

  return sigma;
}



//' @useDynLib CurVol
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double obj_fn_x_cpp(arma::mat d, arma::mat As,
                    arma::mat Bs, arma::mat Gs, arma::mat ysq,
                    arma::mat xsq){

  //initialize
  arma::mat sigma0,y0,prev_sigma,prev_ys,prev_x, new_sigma, sigmat;
  unsigned int cap=1;
  double obj;

  //ysq=NxM
  unsigned int N=size(ysq)(0);
  unsigned int M=size(ysq)(1);

  // Rcout << "ysq " <<size(ysq) << std::endl;

  sigma0=sum(ysq.submat(0,0,4,M-1),0)/5;
  y0=sigma0;


  //compute initial values to start recursion
  //

  prev_sigma.set_size(1,M);
  prev_ys.set_size(1,M);
  prev_x.set_size(1,M);
  //
  // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
  // Rcout << "sigma0 is " <<sigma0 << std::endl;

  for(unsigned int i=0;i<1;i++)
    prev_sigma.row(i)=sigma0;

  for(unsigned int i=0;i<1;i++)
    prev_ys.row(i)=sigma0;

  for(unsigned int i=0;i<1;i++)
    prev_x.row(i)=sigma0;

  // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
  // Rcout << "prev_ys is " <<prev_ys  << std::endl;
  // Rcout << "prev_x is " <<prev_x  << std::endl;

  //get initial 1 sigma inner products:

  new_sigma.set_size(1,M);


  cap=1;

  sigmat.set_size(cap,M);

  for(unsigned int i=0;i<cap;i++){

    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_x_cpp(j,d,As,Bs,Gs,prev_ys.t(),prev_sigma.t(),prev_x.t(),M);
    }
    sigmat.row(i)=new_sigma;

    // Rcout << "prev_sigma is " <<prev_sigma << std::endl;
    // Rcout << "prev_ys is " <<prev_ys  << std::endl;
    // Rcout << "sigmat is " <<sigmat << std::endl;
    //
  }

  // Rcout << "sigmat is " <<sigmat << std::endl;


  for(unsigned int i=cap;i<N;i++){
    //just do this matrix type
    prev_ys=ysq.submat(i-1,0,(i-1),M-1);
    prev_sigma=sigmat.submat(i-1,0,(i-1),M-1);
    prev_x=xsq.submat(i-1,0,(i-1),M-1);
    //recursively compute sigma inner products

    for(unsigned int j=0;j<M;j++){
      new_sigma(0,j)=compute_sigma_inp_k_x_cpp(j,d,As,Bs,Gs,prev_ys.t(),prev_sigma.t(),prev_x.t(),M);
    }

    sigmat=join_vert(sigmat,new_sigma);
  }


  // Rcout << "sigmat is " <<sigmat << std::endl;

  obj=mean(sum(ysq/sigmat+log(sigmat),1));


  return obj;
}
