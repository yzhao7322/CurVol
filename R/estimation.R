#' Estimate Functional ARCH Model
#'
#' @description est.fArch function estimates the Functional ARCH(q) model by using the Quasi-Maximum Likelihood Estimation method.
#'
#' @param fdata The functional data object with N paths.
#' @param basis The M-dimensional basis functions.
#' @param q The order of the depedence on past squared observations. If it is missing, q=1.
#'
#' @return List of model paramters:
#' @return d: d Parameter vector, for intercept function \eqn{\delta}.
#' @return As: A Matrices, for \eqn{\alpha} operators.
#'
#' @export
#'
#' @importFrom nloptr cobyla
#' @importFrom fda Data2fd create.bspline.basis eval.fd
#'
#' @details
#' This function estimates the Functional ARCH(q) model:\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, for \eqn{t \in [0,1]} and \eqn{1\leq i \leq N},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+ \sum_{j=1}^q \int \alpha_j(t,s) x^2_{i-j}(s)ds}.
#'
#' @seealso \code{\link{est.fGarch}} \code{\link{est.fGarchx}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FARCH process and smooth them into a functional data object.
#' yd = dgp.fgarch(grid_point=50, N=200, "arch")
#' yd = yd$garch_mat
#' fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # estimate an FARCH(1) model with basis when M=1.
#' arch1_est = est.fArch(fd, basis_est[,1])
#' }
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
#' Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the ARCH model. Econometric Theory. 29(2), 267-288. <doi:10.1017/S0266466612000345>.\cr
est.fArch=function(fdata, basis, q=1){
  max_eval=10000

  ## projecting the squared process onto given basis functions to get functional scores.
  basis=as.matrix(basis)
  N=dim(fdata$coefs)[2]
  grid_point=dim(basis)[1]
  M=dim(basis)[2]

  int_approx=function(mat){
    temp_n=NROW(mat)
    return((1/temp_n)*sum(mat))}
  fd_eval=eval.fd(seq(0, 1, length.out=grid_point), fdata)

  y_vec=matrix(NA,N,M)
  for(i in 1:N){
    for(j in 1:M){
      y_vec[i,j]=int_approx(fd_eval[,i]^2*basis[,j])
    }
  }
  y_vec = as.matrix(y_vec)

  # Objective function 2
  obj_fn2_q=function(theta,ysq,M,q){

    params=get_dAq_from_theta(theta,M,q)
    d=params$d
    As=params$As

    #for cpp
    As=do.call(rbind,As)
    d=matrix(d,ncol = 1)

    return(obj_fn_q_cpp(d,As,ysq,q))
  }

  # Reformats theta (vector of parameters) into the matrices d, As and Bs.
  get_dAq_from_theta=function(theta,M,q){
    #first M entries are d
    d=theta[1:M]

    #Next p*M^2 entries are A, A_inds are indices of A matrix params
    num_mat_params=(M^2-M*(M-1)/2)
    A_inds=M+(1:(q*(num_mat_params)))

    theta_A=theta[A_inds]

    #This could be changed to search over M^2-M(M-1)/2 paramters instead
    As=list()
    for(i in 1:q){
      curr_A_vals=theta_A[1:num_mat_params+(num_mat_params)*(i-1)]
      A=matrix(0,ncol=M,nrow=M)
      diag(A)=curr_A_vals[1:M]
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
      A=t(A)
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

      As=append(As,list(A))
    }
    return(list(d=d,As=As))
  }

  #If univariate fix type
  if(is.null(dim(y_vec)))
    y_vec=matrix(y_vec,ncol=1)

  #num of basis functions
  M=dim(y_vec)[2]
  #num params in search space
  num_params=q*(M^2-M*(M-1)/2)+M; num_params
  #num of params in B matrices
  num_mat_params=(M^2-M*(M-1)/2)

  #initial values for optimization

  # lower bounds
  d_lower=rep(10^-20,M)
  a_lower=rep(10^-10,q*num_mat_params)

  #constrained optimization

  #formualte objective in terms of theta only
  obj_fn3_q=function(theta){obj_fn2_q(theta=theta,ysq=y_vec,M=M,q=q)}

  #nonlinear inequality constraints on parameters
  eval_g0<-function(theta){
    M = M
    q = q
    h <- numeric(1)
    params=get_dAq_from_theta(theta, M, q)
    upb = 1-(10^-20)
    h[1] <- -(sum(Reduce("+",params$As))-upb)
    return(h)
  }

  #non linear local optimization
  #see nlopt documentation
  #run 10 times to reduce variability in optimization
  for(i in 1:10){
    init=c(runif(num_params,0,1))

    ress=nloptr::cobyla(x0 = init,fn = obj_fn3_q, lower=c(d_lower,a_lower),
                        upper=c(rep(1,num_params)),
                        hin = eval_g0, nl.info = FALSE,
                        control = list(xtol_rel = 1e-18, maxeval=max_eval))
    # print("1")
    if(i==1){
      curr_val=ress$value
      curr_res=ress
      theta=ress$par; #print(res)
      params=get_dAq_from_theta(theta,M,q)

      d=params$d; d
      if(sum(d==0)==M)
        flag=T
      else
        flag=F
    }
    else{
      theta=ress$par; #print(res)
      params=get_dAq_from_theta(theta,M,q)

      d=params$d; d

      if(((ress$value<curr_val)&(sum(d==0)!=M))||flag){
        # print("3")
        flag=FALSE
        curr_val=ress$value
        curr_res=ress
        # print(curr_res)
      }
    }
  }
  # print("4")
  res=curr_res

  #printed in case there was any problems with the optimization
  theta=res$par; print(res)
  params=get_dAq_from_theta(theta,M,q)

  d=params$d; d
  As=params$As; As

  return(params)
}



#' Estimate Functional GARCH Model
#'
#' @description est.fGarch function estimates the Functional GARCH(p,q) model by using the Quasi-Maximum Likelihood Estimation method.
#'
#' @param fdata The functional data object with N paths.
#' @param basis The M-dimensional basis functions.
#' @param p order of the depedence on past volatilities.
#' @param q order of the depedence on past squared observations.
#'
#' @return List of model paramters:
#' @return d: d Parameter vector, for intercept function \eqn{\delta}.
#' @return As: A Matrices, for \eqn{\alpha} operators.
#' @return Bs: B Matrices, for \eqn{\beta} operators.
#'
#' @export
#'
#' @import stats
#' @import sfsmisc
#' @importFrom nloptr cobyla
#' @importFrom fda Data2fd create.bspline.basis eval.fd
#'
#' @seealso \code{\link{est.fArch}} \code{\link{est.fGarchx}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FGARCH process and smooth them into a functional data object.
#' yd = dgp.fgarch(grid_point=50, N=200, "garch")
#' yd = yd$garch_mat
#' fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # estimate an FGARCH(1,1) model with basis when M=1.
#' garch11_est = est.fGarch(fd, basis_est[,1])
#' }
#'
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#' Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
est.fGarch=function(fdata, basis, p=1, q=1){
  max_eval=10000

  ## projecting the squared process onto given basis functions to get functional scores.
  basis=as.matrix(basis)
  N=dim(fdata$coefs)[2]
  grid_point=dim(basis)[1]
  M=dim(basis)[2]

  int_approx=function(mat){
    temp_n=NROW(mat)
    return((1/temp_n)*sum(mat))}
  fd_eval=eval.fd(seq(0, 1, length.out=grid_point), fdata)

  y_vec=matrix(NA,N,M)
  for(i in 1:N){
    for(j in 1:M){
      y_vec[i,j]=int_approx(fd_eval[,i]^2*basis[,j])
    }
  }
  y_vec = as.matrix(y_vec)

  # Objective function 2
  obj_fn2_pq_cpp=function(theta,ysq,M,p,q){

    params=get_dABpq_from_theta(theta,M,p,q)
    d=params$d
    As=params$As
    Bs=params$Bs

    #for cpp
    As=do.call(rbind,As)
    Bs=do.call(rbind,Bs)
    d=matrix(d,ncol = 1)

    #cpp function
    return(obj_fn_pq_cpp(d,As,Bs,ysq,p,q))
  }
  # A function to transform theta

  get_dABpq_from_theta=function(theta,M,p,q){
    #first M entries are d
    d=theta[1:M]

    #Next p*M^2 entries are A, A_inds are indices of A matrix params
    num_mat_params=(M^2-M*(M-1)/2)
    A_inds=M+(1:(q*(num_mat_params)))

    #Rest are Bs..
    B_inds=sfsmisc::last(A_inds)+(1:(p*(num_mat_params)))

    theta_A=theta[A_inds]
    theta_B=theta[B_inds]

    As=list()
    for(i in 1:q){
      curr_A_vals=theta_A[1:num_mat_params+(num_mat_params)*(i-1)]
      A=matrix(0,ncol=M,nrow=M)
      diag(A)=curr_A_vals[1:M]
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
      A=t(A)
      A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

      As=append(As,list(A))
    }

    #see above comment
    Bs=list()
    for(i in 1:p){
      curr_B_vals=theta_B[1:num_mat_params+(num_mat_params)*(i-1)]
      B=matrix(0,ncol=M,nrow=M)
      diag(B)=curr_B_vals[1:M]
      B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
      B=t(B)
      B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]

      Bs=append(Bs,list(B))
    }

    return(list(d=d,Bs=Bs,As=As))
  }

  #If univariate fix type
  if(is.null(dim(y_vec)))
    y_vec=matrix(y_vec,ncol=1)

  #num of basis functions
  M=dim(y_vec)[2]
  #num params in search space
  num_params=(p+q)*(M^2-M*(M-1)/2)+M; num_params
  #num of params in B and A matrices
  num_mat_params=(M^2-M*(M-1)/2)

  # lower bounds
  d_lower=rep(10^-20,M)
  a_lower=rep(10^-10,q*num_mat_params)
  b_lower=rep(10^-20,p*num_mat_params)

  #constrained optimization
  #formualte objective in terms of theta only
  obj_fn3_pq=function(theta,M,p,q){obj_fn2_pq_cpp(theta=theta,ysq=y_vec,M=M,p=p,q=q)}

  #nonlinear inequality constraints on parameters
  eval_g0<-function(theta,M,p,q){
    h <- numeric(2)
    params=get_dABpq_from_theta(theta,M,p,q)
    upb = 1-(10^-20)
    h[1] <- -(sum(Reduce("+",params$As))-upb)
    h[2] <- -(sum(Reduce("+",params$Bs))-upb)
    #h[1] <- -(sum(Reduce("+",params$As))+sum(Reduce("+",params$Bs))-upb)
    return(h)
  }

  #non linear local optimization
  #see nlopt documentation
  #run 10 times to reduce variability in optimization
  for(i in 1:10){
    init=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))

    ress=nloptr::cobyla(x0 = init,fn = obj_fn3_pq, lower=c(d_lower,a_lower,b_lower),
                        upper=c(rep(1,num_params)),
                        hin = eval_g0, nl.info = FALSE,
                        control = list(maxeval=max_eval), M = M, p = p,q = q)
    # print("1")
    if(i==1){
      curr_val=ress$value
      curr_res=ress
      theta=ress$par; #print(res)
      params=get_dABpq_from_theta(theta,M,p,q)

      d=params$d; d
      if(sum(d==0)==M)
        flag=T
      else
        flag=F
    }

    else{
      theta=ress$par; #print(res)
      params=get_dABpq_from_theta(theta,M,p,q)

      d=params$d; d

      if(((ress$value<curr_val)&(sum(d==0)!=M))||flag){
        # print("3")
        flag=FALSE
        curr_val=ress$value
        curr_res=ress
        # print(curr_res)
      }
    }
  }
  # print("4")
  res=curr_res

  #printed in case there was any problems with the optimization
  theta=res$par; print(res)
  params=get_dABpq_from_theta(theta,M,p,q)

  d=params$d; d
  As=params$As; As
  Bs=params$Bs; Bs

  return(params)
}


#' Estimate Functional GARCH-X Model
#'
#' @description est.fGarchx function estimates the Functional GARCH-X model by using the Quasi-Maximum Likelihood Estimation method.
#'
#' @param fdata_y The functional data object with N paths for the objective data.
#' @param fdata_x The functional data object with N paths for the covariate X.
#' @param basis The M-dimensional basis functions.
#'
#'
#' @return List of model paramters:
#' @return d: d Parameter vector, for intercept function \eqn{\delta}.
#' @return As: A Matrices, for \eqn{\alpha} operators.
#' @return Bs: B Matrices, for \eqn{\beta} operators.
#' @return Gs: G Matrices, for \eqn{\gamma} operators.
#'
#' @export
#'
#' @import sfsmisc
#' @importFrom nloptr cobyla
#' @importFrom fda Data2fd create.bspline.basis eval.fd
#'
#' @details
#' This function estimates the Functional GARCH-X model:\cr
#' \eqn{y_i(t)=\sigma_i(t)\varepsilon_i(t)}, for \eqn{t \in [0,1]} and \eqn{1\leq i \leq N},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) y^2_{i-1}(s)ds+\int \beta(t,s) \sigma^2_{i-1}(s)ds + \int \theta(t,s) x^2_{i-1}(s)ds},\cr
#' where \eqn{x_i(t)} is an exogenous variable.
#'
#' @seealso \code{\link{est.fArch}} \code{\link{est.fGarch}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations and smooth them into functional data objects.
#' yd = dgp.fgarch(grid_point=50, N=200, "arch")
#' yd = yd$garch_mat
#' xd = dgp.fgarch(grid_point=50, N=200, "garch")
#' xd = xd$garch_mat
#' fdy = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#' fdx = fda::Data2fd(argvals=seq(0,1,len=50),y=xd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # estimate an FGARCH-X model with basis when M=1.
#' garchx_est = est.fGarchx(fdy, fdx, basis_est[,1])
#' }
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2021) Exploring volatility of crude oil intra-day return curves: a functional GARCH-X model. MPRA Paper No. 109231. <https://mpra.ub.uni-muenchen.de/109231>. Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
#'
est.fGarchx=function(fdata_y, fdata_x, basis){
  max_eval=10000
  ## projecting the squared process onto given basis functions to get functional scores.
  basis=as.matrix(basis)
  N=dim(fdata_y$coefs)[2]
  grid_point=dim(basis)[1]
  M=dim(basis)[2]

  int_approx=function(mat){
    temp_n=NROW(mat)
    return((1/temp_n)*sum(mat))}
  fd_eval_y=eval.fd(seq(0, 1, length.out=grid_point), fdata_y)
  fd_eval_x=eval.fd(seq(0, 1, length.out=grid_point), fdata_x)

  y_vec = x_vec = matrix(NA,N,M)
  for(i in 1:N){
    for(j in 1:M){
      y_vec[i,j]=int_approx(fd_eval_y[,i]^2*basis[,j])
      x_vec[i,j]=int_approx(fd_eval_x[,i]^2*basis[,j])
    }
  }
  y_vec = as.matrix(y_vec)
  x_vec = as.matrix(x_vec)

  # Objective function 2
  obj_fn2=function(theta,ysq,xs,M){

    params=get_dABG_from_theta(theta,M)
    d=params$d
    As=params$As
    Bs=params$Bs
    Gs=params$Gs

    #for cpp
    As=do.call(rbind,As)
    Bs=do.call(rbind,Bs)
    Gs=do.call(rbind,Gs)
    d=matrix(d,ncol = 1)

    return(obj_fn_x_cpp(d,As,Bs,Gs,ysq,xs))
  }

  get_dABG_from_theta=function(theta,M){
    #first M entries are d
    d=theta[1:M]

    #Next p*M^2 entries are A, A_inds are indices of A matrix params
    num_mat_params=(M^2-M*(M-1)/2)
    A_inds=M+(1:(num_mat_params))

    #Rest are Bs..
    B_inds=sfsmisc::last(A_inds)+(1:(num_mat_params))
    G_inds=sfsmisc::last(B_inds)+(1:(num_mat_params))

    theta_A=theta[A_inds]
    theta_B=theta[B_inds]
    theta_G=theta[G_inds]
    #For each M^2 parameters, make As pos def by t(V)%*%V
    #Not sure if As need to be pos def, but it works..
    #So the search space is over the V_i such that Ai=Vi'Vi
    #This could be changed to search over M^2-M(M-1)/2 paramters instead
    curr_A_vals=theta_A[1:num_mat_params]
    A=matrix(0,ncol=M,nrow=M)
    diag(A)=curr_A_vals[1:M]
    A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
    A=t(A)
    A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
    As=list()
    As=append(As,list(A))

    #see above comment
    curr_B_vals=theta_B[1:num_mat_params]
    B=matrix(0,ncol=M,nrow=M)
    diag(B)=curr_B_vals[1:M]
    B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
    B=t(B)
    B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
    Bs=list()
    Bs=append(Bs,list(B))

    curr_G_vals=theta_G[1:num_mat_params]
    G=matrix(0,ncol=M,nrow=M)
    diag(G)=curr_G_vals[1:M]
    G[lower.tri(G)]=curr_G_vals[(M+1):length(curr_G_vals)]
    G=t(G)
    G[lower.tri(G)]=curr_G_vals[(M+1):length(curr_G_vals)]
    Gs=list()
    Gs=append(Gs,list(G))

    return(list(d=d,As=As,Bs=Bs,Gs=Gs))
  }

  #If univariate fix type
  if(is.null(dim(y_vec)))
    y_vec=matrix(y_vec,ncol=1)

  #num of basis functions
  M=dim(y_vec)[2]
  #num params in search space
  num_params=3*(M^2-M*(M-1)/2)+M; num_params
  #num of params in B, G and A matrices
  num_mat_params=(M^2-M*(M-1)/2)

  # lower bounds
  d_lower=rep(10^-20,M)
  ag_lower=rep(10^-10,num_mat_params)
  b_lower=rep(10^-20,num_mat_params)
  #constrained optimization
  #formualte objective in terms of theta only
  obj_fn3=function(theta){obj_fn2(theta=theta,ysq=y_vec,xs=x_vec,M=M)}

  #nonlinear inequality constraints on parameters
  eval_g0<-function(theta){
    M = M
    h <- numeric(2)
    params=get_dABG_from_theta(theta, M)
    upb = 1-(10^-20)
    h[1] <- -(sum(Reduce("+",params$As))-upb)
    h[2] <- -(sum(Reduce("+",params$Bs))-upb)
    return(h)
  }

  #non linear local optimization
  #see nlopt documentation
  #run 10 times to reduce variability in optimization
  for(i in 1:10){
    init=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))

    ress=nloptr::cobyla(x0 = init,fn = obj_fn3, lower=c(d_lower,ag_lower,b_lower,ag_lower),
                        upper=c(rep(1,num_params)),
                        hin = eval_g0, nl.info = FALSE,
                        control = list(xtol_rel = 1e-18, maxeval=max_eval))
    # print("1")
    if(i==1){
      curr_val=ress$value
      curr_res=ress
      theta=ress$par; #print(res)
      params=get_dABG_from_theta(theta,M)

      d=params$d; d
      if(sum(d==0)==M)
        flag=T
      else
        flag=F
    }

    else{
      theta=ress$par; #print(res)
      params=get_dABG_from_theta(theta,M)

      d=params$d; d

      if(((ress$value<curr_val)&(sum(d==0)!=M))||flag){
        # print("3")
        flag=FALSE
        curr_val=ress$value
        curr_res=ress
        # print(curr_res)
      }
    }
  }
  # print("4")
  res=curr_res

  #printed in case there was any problems with the optimization
  theta=res$par; print(res)
  params=get_dABG_from_theta(theta,M)

  return(params)
}





#' Diagnostic information derived from the estimation
#'
#' @description diagnostic.fGarch function provides the estimation parameters that can be used as the inputs for a diagnostic purpose.
#'
#' @param params List of model paramters: d, As, Bs (and Gs).
#' @param basis The M-dimensional basis functions for the projection.
#' @param yd A grid_point x N matrix drawn from N functional curves.
#' @param p The order of the depedence on past squared observations. If it is missing, p=0.
#' @param q The order of the depedence on past volatilities. If it is missing, q=1.
#' @param xd A grid_point x N matrix drawn from N covariate X curves. The default value is NULL.
#'
#' @return List of parameters:
#' @return eps: a grid_point x N matrix containing fitted residuals.
#' @return sigma2: a grid_point x (N+1) matrix that the first N columns are the fitted conditional variance, the N+1 column is the predicted conditional variance.
#' @return yfit: a grid_point x N matrix drawn from N fitted intra-day price curves.
#' @return kernel_op: estimated kernel coefficient operators.
#'
#' @export
#'
#' @seealso \code{\link{est.fArch}} \code{\link{est.fGarch}} \code{\link{est.fGarchx}} \code{\link{gof.fgarch}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FGARCH process and smooth them into a functional data object.
#' yd = dgp.fgarch(grid_point=50, N=200, "garch")
#' yd = yd$garch_mat
#' fdy = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # fit the curve data with an FARCH(1) model.
#' arch_est = est.fArch(fdy, basis_est)
#'
#' # get parameters for diagnostic checking.
#' diag_arch  = diagnostic.fGarch(arch_est, basis_est, yd)
#' }
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020) Forecasting Value at Risk via intra-day return curves. International Journal of Forecasting. <doi:10.1016/j.ijforecast.2019.10.006>.\cr
#'
diagnostic.fGarch=function(params,basis,yd,p=0,q=1,xd=NULL){
  grid_point = nrow(yd)
  N=ncol(yd)
  basis = as.matrix(basis)
  M = ncol(basis)

  d=params$d
  As=params$As
  if(p!=0){
    Bs=params$Bs}
  if( is.null(xd) == FALSE){
    Gs=params$Gs
  }

  int_approx = function(mat){
    temp_n=NROW(mat)
    return((1/temp_n)*sum(mat))}
  positive_c <- function(x) { x[x<0] <- 0; x}

  times = 1:grid_point/grid_point
  delta_hat=0
  for (i in 1:M){
    delta_hat=delta_hat+d[i]*basis[,i]
  }

  alpha_Op_hat=list()
  for (h in 1:q){
    alpha_M=matrix(0,NROW(times),NROW(times))
    A_mat=as.matrix(As[[h]])
    temp_alpha_m=0
    for(t in 1:M){
      for(s in 1:M){
        temp_alpha_m=temp_alpha_m+A_mat[t,s]*basis[,t]%*%t(basis[,s])
      }
    }
    alpha_M=temp_alpha_m
    alpha_Op_hat[[h]]=positive_c(alpha_M)
  }

  if(p!=0){
    beta_Op_hat=list()
    for (h in 1:p){
      beta_M=matrix(0,NROW(times),NROW(times))
      B_mat=as.matrix(Bs[[h]])
      temp_beta_m=0
      for(t in 1:M){
        for(s in 1:M){
          temp_beta_m=temp_beta_m+B_mat[t,s]*basis[,t]%*%t(basis[,s])
        }
      }
      beta_M=temp_beta_m
      beta_Op_hat[[h]]=positive_c(beta_M)
    }
  }

  if(is.null(xd) == FALSE){
    gamma_M=matrix(0,NROW(times),NROW(times))
    G_mat=as.matrix(Gs[[1]])
    temp_gamma_m=0
    for(t in 1:M){
      for(s in 1:M){
        temp_gamma_m=temp_gamma_m+G_mat[t,s]*basis[,t]%*%t(basis[,s])
      }
    }
    gamma_M=temp_gamma_m
    gamma_Op_hat=positive_c(gamma_M)
  }

  sigma2_fit=matrix(NA,grid_point,(N+1))
  for(i in 1:max(p,q)){
    sigma2_fit[,i]=delta_hat
  }

  for(i in (max(p,q)+1):(N+1)){
    #first fill in sigma2 column:
    for(j in 1:grid_point){
      fit_arch_op = fit_garch_op = 0
      for (h in 1:q){
        fit_arch_op = fit_arch_op + int_approx(alpha_Op_hat[[h]][j,] * yd[,i-h]^2)}

      if(p!=0){
        for (h in 1:p){
          fit_garch_op = fit_garch_op + int_approx(beta_Op_hat[[h]][j,] * sigma2_fit[,i-h])}
      }
      sigma2_fit[j,i] = delta_hat[j] + fit_arch_op + fit_garch_op
      if( is.null(xd) ==FALSE){
        sigma2_fit[j,i] = delta_hat[j] + fit_arch_op + fit_garch_op + int_approx(gamma_Op_hat[j,] * xd[,i-1])
      }
    }
  }

  error_fit=yd/sqrt(abs(sigma2_fit[,1:N]))
  error_fit[,1]=0
  # simulate an error following OU process
  error_sim <- function(grid_point, samplesize){
    ti=1:grid_point/grid_point
    comat=matrix(NA,grid_point,grid_point)
    for (i in 1:grid_point){
      comat[i,]=exp(1)^(-ti[i]/2-ti/2)*pmin(exp(1)^(ti[i]),exp(1)^(ti))
    }
    epsilon=MASS::mvrnorm(n = samplesize, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)

    return(t(epsilon))
  }
  fd_fit=sqrt(abs(sigma2_fit[,1:N]))*error_sim(grid_point,N)

  # kernel coefficients
  kernel_coef = list(delta = delta_hat, alpha = alpha_Op_hat)
  if(p!=0){kernel_coef = list(delta = delta_hat, alpha = alpha_Op_hat, beta = beta_Op_hat)}
  if(is.null(xd) == FALSE){kernel_coef = list(delta = delta_hat, alpha = alpha_Op_hat, beta = beta_Op_hat, gamma = gamma_Op_hat)}

  return(list(eps = error_fit, sigma2 = sigma2_fit, yfit = fd_fit, kernel_op = kernel_coef))
}
