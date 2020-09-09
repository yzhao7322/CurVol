#' Intra-day Return Curves
#'
#' @description intra.return function calculates three types of intra-day return curves: intra-day log return (IDR), cumulative intra-day log return (CIDR), and overnight cumulative intra-day log return (OCIDR).
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N intra-day price curves.
#'
#' @return List of return data:
#' @return idr: the IDR.
#' @return cidr: the CIDR.
#' @return ocidr: the OCIDR.
#'
#' @export
#'
#' @details
#' Suppose that \eqn{P_i(t)} denotes the intra-day price curves, for \eqn{1\leq i \leq N} and \eqn{t\in[0,1]}, we then have,\cr
#' Intra-day log returns: \eqn{x_i(t)=\log P_i(t) - \log P_i(t-\Delta)}, where \eqn{\Delta} is the frequency at intra-day grid points;
#' Cumulative intra-day log returns: \eqn{X_i(t)=\log P_i(t) - \log P_i(0)};
#' Overnight cumulative intra-day log returns: \eqn{x_i(t)=\log P_i(t) - \log P_{i-1}(1)}.
#'
#' @examples
#' # generate intra-day price curve data for the FARCH process.
#' yd = dgp.farch(50,100)
#' yd_arch = yd$arch_mat
#'
#' # calculate discrete data drawn from N intra-day return curves.
#' fcurve = intra.return(yd_arch)
#' idr = fcurve$idr
#' cidr = fcurve$cidr
#' ocidr = fcurve$ocidr
#'
#' @references
#' Gabrys, R., Horvath, L., Kokoszka, P. (2010). Tests for error correlation in the functional linear model. Journal of the American Statistical Association, 105(491), 1113-1125.\cr
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Forecasting Value at Risk via Intra-Day Return Curves. To appear in International Journal of Forecasting.
intra.return <- function(yd){
  grid_point=nrow(yd)
  N=ncol(yd)
  idr_re=matrix(0,grid_point-1,N)
  cidr_re=matrix(0,grid_point,N)
  ocidr_re=matrix(0,grid_point,N-1)

  for(j in c(2:grid_point)){
    for(i in c(1:N)){
      idr_re[j-1,i]=log(yd[j,i])-log(yd[j-1,i])
    }
  }

  for(j in c(1:grid_point)){
    for(i in c(1:N)){
      cidr_re[j,i]=log(yd[j,i])-log(yd[1,i])
    }
  }

  for(j in c(1:grid_point)){
    for(i in c(2:N)){
      ocidr_re[j,i-1]=log(yd[j,i])-log(yd[grid_point,i-1])
    }
  }

  return(list(idr = idr_re, cidr = cidr_re, ocidr = ocidr_re))
}




#' Positive Polynominals
#'
#' @description basis.pp function returns positive vectors drawn from exponential and Bernstein functions.
#' @param grid_point The number of grid point in each curve observation.
#' @param M The number/order of basis functions.
#'
#' @return List of objects:
#' @return exp: A (grid_point) x (M) matrix drawn from M Exponential basis functions.
#' @return bern: A (grid_point) x (M) matrix drawn from M Bernstein basis functions.
#'
#' @export
#' @importFrom stats time
#'
#' @seealso \code{\link{basis.tfpca}} \code{\link{basis.fsnn}} \code{\link{basis.pf}}
#' @examples
#' #generate discrete evaluations of exponential and Bernstein basis functions with order one.
#' ppb = basis.pp(50,1)
#' exp = ppb$exp
#' bern = ppb$bern
basis.pp <- function(grid_point,M){
  time=1:grid_point/grid_point
  exp_matrix=matrix(NA,grid_point,M)
  bernstein_matrix=matrix(NA,grid_point,M)
  for (k in 1:M){
    exp_matrix[,k]=exp(1)^(k*time)
    bernstein_matrix[,k]=(factorial(M-1)/(factorial(k-1)*factorial(M-k)))*(time^(k-1))*((1-time)^(M-k))
  }
  return(list(exp = exp_matrix, bern = bernstein_matrix))
}


#' Truncated Functional Pricipal Component Analysis
#'
#' @description basis.tfpca function extracts truncated non-negative functional principal components from the squared process of conditionally heteroscedastic curves, presumably, so that they may be modelled with a Functional GARCH model subsequently.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N discrete evaluation curves.
#' @param M The number/order of basis functions, by default \eqn{M = 10}.
#'
#' @return List of objects:
#' @return basis: a (grid_point) x (M) matrix containing discrete evaluations of M truncated functional principal components.
#' @return tve: the total variations explained by each truncated functional principal component.
#'
#' @export
#'
#' @import fda
#'
#' @seealso \code{\link{basis.pp}} \code{\link{basis.fsnn}} \code{\link{basis.pf}}
#' @examples
#' # generate discrete evaluations of the FARCH process.
#' fd = dgp.farch(50,100)
#' fd_arch = fd$arch_mat
#'
#' # decompose the first truncated non-negative functional principal component.
#' dt = basis.tfpca(fd_arch,M=1)
#' tbasis = dt$basis
#' tve = dt$tve
#'
#' @references
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of econometrics, 209(2), 353-375.
basis.tfpca <- function(yd,M){
  grid_point=nrow(yd)
  N=ncol(yd)
  yd=yd-rowMeans(yd)

  times=rep(0,grid_point)
  for(i in 1:grid_point){times[i]=i/grid_point}
  squared_y=yd*yd
  basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=32,norder=4)
  y_sq_fdsmooth=smooth.basis(argvals=times,y=squared_y,fdParobj=basis_obj)
  y_sq_fd=y_sq_fdsmooth$fd
  pca_obj=pca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd), centerfns=TRUE)  #centerfns=FALSE.
  eigen_functs=pca_obj$harmonics
  eigen_v=pca_obj$values
  tve=eigen_v/sum(eigen_v)
  ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
  for(j in 1:10){indi_basis=eval.fd(times,eigen_functs[j])
  indi_basis[indi_basis<0]<-0
  ortho_basis_matrix[,j]=indi_basis}
  return(list(basis = ortho_basis_matrix[,1:M], tve = tve))
}


#' Sparse and Non-Negative Functional Principal Component Analysis
#'
#' @description basis.fsnn function extracts sparse and non-negative functional principal components from the squared process of intra-day returns.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N curves.
#' @param M The number/order of basis functions, by default \eqn{M = 10}.
#'
#' @return List of objects:
#' @return basis: a (grid_point) x (M) matrix containing discrete evaluations of M sparse and non-negative functional principal components.
#' @return tve: the total variations explained by each sparse and non-negative functional principal component.
#'
#' @export
#'
#' @importFrom nsprcomp nsprcomp
#' @import fda
#'
#' @seealso \code{\link{basis.pp}} \code{\link{basis.tfpca}} \code{\link{basis.pf}}
#' @examples
#' # generate discrete evaluations of the FARCH process.
#' fd = dgp.farch(50,100)
#' fd_arch = fd$arch_mat
#'
#' # decompose the first sparse and non-negative functional principal component.
#' basis.fsnn(fd_arch,M=1)
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Functional GARCH-X model with an application to forecasting crude oil return curves. Working Paper.
basis.fsnn <- function(yd, M){
  grid_point = nrow(yd)
  N = ncol(yd)
  yd = yd-rowMeans(yd)

  times = rep(0,grid_point)
  for(i in 1:grid_point){times[i] = i/grid_point}
  squared_y = yd*yd

  nspca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
                       centerfns = TRUE,samplesize){
    meanfd = mean.fd(fdobj)
    if (centerfns) {
      fdobj = center.fd(fdobj)
    }

    #  get coefficient matrix and its dimensions

    coef = fdobj$coefs

    coefd = dim(coef) # dim=[samplesize,nbasis]
    ndim  = length(coefd)
    nrep  = coefd[2]
    coefnames = dimnames(coef)
    if (nrep < 2) stop("PCA not possible without replications.")

    basisobj = fdobj$basis
    nbasis = basisobj$nbasis
    type = basisobj$type

    #  set up HARMBASIS

    harmbasis = harmfdPar$fd$basis
    nhbasis = harmbasis$nbasis

    #  set up LFDOBJ and LAMBDA

    Lfdobj = harmfdPar$Lfd
    lambda = harmfdPar$lambda

    #  compute CTEMP whose cross product is needed

    if (ndim == 3) {
      nvar = coefd[3]
      ctemp = matrix(0, nvar * nbasis, nrep)
      for(j in 1:nvar) {
        index = 1:nbasis + (j - 1) * nbasis
        ctemp[index,  ] = coef[,  , j]
      }
    } else {
      nvar = 1
      ctemp = coef
    }

    indi_pc=nsprcomp(t(ctemp),k=nbasis,nneg=TRUE,center = FALSE)
    eigvalc=(indi_pc$sdev)^2
    eigvecc=indi_pc$rotation[, 1:nharm] # check the defination of nsprcomp

    sumvecc = apply(eigvecc, 2, sum)
    eigvecc[,sumvecc < 0] =  - eigvecc[, sumvecc < 0]

    varprop = eigvalc[1:nharm]/sum(eigvalc)

    #  set up harmfd

    if (nvar == 1) {
      harmcoef = eigvecc
    } else {
      harmcoef = array(0, c(nbasis, nharm, nvar))
      for (j in 1:nvar) {
        index = 1:nbasis + (j - 1) * nbasis
        temp = eigvecc[index,  ]
        harmcoef[,  , j] =  temp #slove(Wmat)%*% temp
      }
    }

    harmnames = rep("", nharm)
    for(i in 1:nharm)
      harmnames[i] = paste("PC", i, sep = "")
    if(length(coefd) == 2)
      harmnames = list(coefnames[[1]], harmnames,"values")
    if(length(coefd) == 3)
      harmnames = list(coefnames[[1]], harmnames, coefnames[[3]])
    harmfd = fd(harmcoef, harmbasis, harmnames)

    #  set up harmscr

    if (nvar == 1) {
      harmscr = inprod(fdobj, harmfd)
    } else {
      harmscr = array(0, c(nrep,   nharm, nvar))
      coefarray = fdobj$coefs
      harmcoefarray = harmfd$coefs
      for (j in 1:nvar) {
        fdobjj= fd(as.matrix(    coefarray[,,j]), basisobj)
        harmfdj = fd(as.matrix(harmcoefarray[,,j]), basisobj)
        harmscr[,,j] = inprod(fdobjj, harmfdj)
      }
    }

    pcafd = list(harmfd, eigvalc, harmscr, varprop, meanfd)
    class(pcafd) = "pca.fd"
    names(pcafd) = c("harmonics", "values", "scores", "varprop", "meanfd")

    return(pcafd)
  }

  # smooth data into functional curves
  basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=32,norder=4)
  y_sq_fdsmooth=smooth.basis(argvals=times,y=squared_y,fdParobj=basis_obj)
  y_sq_fd=y_sq_fdsmooth$fd

  # get sparse and non-negative basis
  pca_obj=nspca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd),centerfns = TRUE) # by default to compute the first 10 principal components
  eigen_functs=pca_obj$harmonics
  eigen_v=pca_obj$values
  tve=eigen_v/sum(eigen_v)
  ortho_basis_matrix=matrix(0,nrow=grid_point,ncol=10)
  for(j in 1:10){ortho_basis_matrix[,j]=eval.fd(times,eigen_functs[j])}

  return(list(basis = ortho_basis_matrix[,1:M], tve = tve))
}





#' Predictive Factors
#'
#' @description basis.pf function extracts non-negative predictive factors from the squared process of conditionally heteroscedastic curves, presumably, so that they may be modelled with a Functional GARCH model subsequently.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N discrete evaluation curves.
#' @param M The number/order of basis functions, by default \eqn{M = 10}.
#' @param a The regularization parameter with recommended values of 0.73 - 0.75.
#' @param h The forecasting horizon.
#'
#' @return List of objects:
#' @return basis: a (grid_point) x (M) matrix containing discrete evaluations of M non-negative predictive factors.
#' @return tve: the total variations explained by each non-negative predictive factor.
#'
#' @export
#' @importFrom QZ qz.geigen
#'
#' @seealso \code{\link{basis.pp}} \code{\link{basis.tfpca}} \code{\link{basis.fsnn}}
#' @examples
#' # generate discrete evaluations of the FARCH process.
#' fd = dgp.farch(50,100)
#' fd_arch = fd$arch_mat
#'
#' # decompose the first two non-negative predictive factors.
#' basis.pf(fd_arch,M=2,a=0.75,h=1)
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Functional GARCH-X model with an application to forecasting crude oil return curves. Working Paper.
basis.pf<-function(yd,M,a,h){

  est_cov<-function(yd,h){
    #h: forecast horizon
    yd=yd-rowMeans(yd)
    # dimension
    n = ncol(yd)
    # Calculate Autocovariance
    C = yd%*%t(yd)/n
    # Calcualte Crosscovariance
    C1 = yd[,((h+1):n)]%*%t(yd[,(1:(n-h))])/(n-h)
    return(list(C, C1))
  }

  d = nrow(yd)
  C_est=est_cov(yd,h)
  C=C_est[[1]]
  C1=C_est[[2]]

  # Calculate Ca
  Ca = C+a*diag(d)

  # Find Predictive Factors
  aaa = t(C1)%*%C1

  out = qz.geigen(aaa, B = Ca, only.values = FALSE)
  # eigen values and vectors
  D1 = out$ALPHA
  V1 = out$V

  Dhat1=D1[1:M]
  Vhat1=V1[,1:M]

  VV=sqrt(1./diag(t(Vhat1)%*%Ca%*%Vhat1))
  Vhat1=Vhat1%*%diag(VV)
  R=Vhat1

  F=C1%*%R

  # truncate the negative part
  if(any(F[,1]<0)==TRUE){
    F=-F}
  F[F<0]<-0

  tve=as.vector(rep(NA,length(D1)))
  for (i in 1:length(D1)){
    tve[i]=D1[i]/sum(D1)
  }

  return(list(basis = F, tve = tve))
}
