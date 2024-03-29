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
#' # generate intra-day price curve data for the FGARCH process.
#' yd = dgp.fgarch(50, 100, "garch")
#' yd = yd$garch_mat
#'
#' # calculate discrete data drawn from N intra-day return curves.
#' fcurve = intra.return(yd)
#' idr = fcurve$idr
#' cidr = fcurve$cidr
#' ocidr = fcurve$ocidr
#'
#' @references
#' Gabrys, R., Horvath, L., Kokoszka, P. (2010). Tests for error correlation in the functional linear model. Journal of the American Statistical Association, 105(491), 1113-1125. <doi:10.1198/jasa.2010.tm09794>.\cr
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Forecasting Value at Risk via Intra-Day Return Curves. International Journal of Forecasting. <doi:10.1016/j.ijforecast.2019.10.006>.\cr
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




#' Estimating non-negative basis functions
#'
#' @description basis.est function generates/estimates non-negative basis functions, so that they can be used to estimate a Functional GARCH model subsequently.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N discrete evaluation curves.
#' @param M The number/order of basis functions, by setting \eqn{M <= 10}.
#' @param type A string to switch estimation method: "poly" - exponential and Bernstein functions; "tfpca" - truncated functional principal components; "fsnn" - sparse and non-negative functional principal components; "tpf" - truncated predictive factors.
#'
#' @return List of objects:
#' @return basis: a (grid_point) x (M) matrix containing discrete evaluations of M data-driven bases; when "poly" is applied, the basis contains a list with 'exp' for exponential and 'bern' for Bernstein.
#' @return tve: the total variations explained by data-driven bases, not applied for polynomials.
#'
#' @export
#'
#' @import fda
#' @importFrom stats time
#' @importFrom nsprcomp nsprcomp
#' @importFrom QZ qz.geigen
#'
#' @examples
#' # generate discrete evaluations of the FGARCH process.
#' yd = dgp.fgarch(50, 100, "garch")
#' yd = yd$garch_mat
#'
#' # decompose the first truncated non-negative functional principal component.
#' dt = basis.est(yd, M=1, "tfpca")
#' tbasis = dt$basis
#' tve = dt$tve
#'
#' @references
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
#' Rice, G., Wirjanto, T., Zhao, Y. (2021) Exploring volatility of crude oil intra-day return curves: a functional GARCH-X model. MPRA Paper No. 109231. <https://mpra.ub.uni-muenchen.de/109231>.\cr
basis.est <- function(yd, M, type){
  grid_point=nrow(yd)
  N=ncol(yd)
  yd=yd-rowMeans(yd)

  switch(type,
         poly = {
           time=1:grid_point/grid_point
           exp_matrix=matrix(NA,grid_point,M)
           bernstein_matrix=matrix(NA,grid_point,M)
           for (k in 1:M){
             exp_matrix[,k]=exp(1)^(k*time)
             bernstein_matrix[,k]=(factorial(M-1)/(factorial(k-1)*factorial(M-k)))*(time^(k-1))*((1-time)^(M-k))
           }
           basis1 = list(exp = exp_matrix, bern = bernstein_matrix)
           tve = NULL
         },
         tfpca = {
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
           basis1 = ortho_basis_matrix[,1:M]
         },
         fsnn = {
           times = rep(0,grid_point)
           for(i in 1:grid_point){times[i] = i/grid_point}
           squared_y = yd*yd

           nspca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj), centerfns = TRUE,samplesize){
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
           basis1 = ortho_basis_matrix[,1:M]
         },
         tpf = {
           # the regularization parameter (tuning parameter)
           a = 0.75
           # the forecasting horizon
           h = 1
           est_cov<- function(yd,h){
             if (class(yd) != 'matrix'){stop('Input data type is not matrix')}
             # variance operator
             N = ncol(yd)
             C = (yd %*% t(yd)) / N
             # covaraince operator
             CC = (yd[,(h+1):N] %*% t(yd[,1:(N-h)])) / (N-h)
             return(list(C,CC))}

           K = nrow(yd)
           # variance and corvariance operator
           Temp1 = est_cov(yd,h)
           C = as.matrix(data.frame(Temp1[1]))
           CC = as.matrix(data.frame(Temp1[2]))
           # empirical operator calculation
           Ca = C + a * diag(K)
           # Find the predictive factor
           aaa = CC %*% t(CC)

           out = qz.geigen(aaa, B = Ca, only.values = FALSE)
           # eigen values and vectors
           eigen.value = out$ALPHA
           eigen.vector = out$V

           # order and sort the eigen value
           id = order(eigen.value,decreasing=TRUE)
           eigen.value = sort(eigen.value,decreasing=TRUE)
           eigen.vector.new = eigen.vector[,id]

           # Final estimation
           eigen.value.hat = eigen.value[1:M]
           eigen.vector.hat = eigen.vector.new[,1:M]

           vv = sqrt(1/(diag(t(eigen.vector.hat) %*% Ca %*% eigen.vector.hat)))
           vv.matrix = diag(0,length(vv),length(vv))
           diag(vv.matrix) = vv

           eigen.vector.hat = eigen.vector.hat %*% vv.matrix
           R = eigen.vector.hat
           FF = CC %*% R
           # truncate the negative part
           if(any(FF[,1]<0)==TRUE){
             FF[,1]=-FF[,1]}
           if(any(FF[,2]<0)==TRUE){
             FF[,2]=-FF[,2]}
           FF[FF<0]<-0

           tve=as.vector(rep(NA,length(eigen.value)))
           for (i in 1:length(eigen.value)){
             tve[i]=eigen.value[i]/sum(eigen.value)
           }
           basis1 = FF
         },
         stop("Enter something to switch me!"))
  return(list(basis = basis1, tve = tve))
}

