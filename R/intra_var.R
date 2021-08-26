#' Forecast Daily/intra-day Value-at-Risk
#'
#' @description var.forecast function forecasts the daily VaR and intra-day VaR curves according to intra-day return curves.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N functional curves.
#' @param sigma_pred The predicted conditional standard deviation curve obtained by using the functional ARCH/GARCH model.
#' @param error_fit The residual obtained from the functional ARCH/GARCH model.
#' @param quantile_v The percentile of the VaR forecasts.
#' @param Method A string to indicate which method will be implemented to estimate the quantile of the error: "normal" - assuming errors following a normal distribution; "bootstrap" - bootstrapping the empirical distribution from the historical data; "evt" - fitting the error with a generalized Pareto distribution.
#'
#' @return List of objects:
#' @return day_VaR: the daily VaR.
#' @return day_ES: the daily expected shortfall.
#' @return intraday_VaR: the intra-day VaR curves.
#' @return vio_seq: the violation curve for the intra-day VaR curves.
#'
#' @export
#'
#' @importFrom POT fitgpd
#' @importFrom POT rgpd
#'
#' @details
#' This function uses a two-step approach to forecast intra-day Value-at-Risk, in formula that\cr
#' \eqn{VaR_{i+1}^\tau(t)=\hat{\sigma}_{i+1}(t)\hat{\varepsilon}^\tau(t)}, for \eqn{1\leq i \leq N}, \eqn{t\in[0,1]}, and the percentile \eqn{\tau \in [0,1]},\cr
#' where the forecats of conditional standard deviation \eqn{\hat{\sigma}_{i+1}(t)} can be obtained by using \code{\link{est.fArch}}, \code{\link{est.fGarch}}, or \code{\link{est.fGarchx}}. Note that when \eqn{t=1}, \eqn{VaR_{i+1}^\tau(1)} is the forecast of daily VaR.
#'
#' @seealso \code{\link{basis.est}} \code{\link{est.fGarch}} \code{\link{diagnostic.fGarch}} \code{\link{var.backtest}} \code{\link{var.vio}}
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FGARCH(1,1) process.
#' grid_point = 50; N = 200
#' yd = dgp.fgarch(grid_point, N, "garch")
#' yd = yd$garch_mat
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # fit the curve data and the conditional volatility by using an FGARCH(1,1) model with M=1.
#' fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))
#' garch11_est = est.fGarch(fd, basis_est[,1])
#' diag_garch = diagnostic.fGarch(garch11_est, basis_est[,1], yd)
#'
#' # get the in-sample fitted conditional variance and error.
#' sigma_fit = diag_garch$sigma2[,1:N]
#' error_fit = diag_garch$eps
#'
#' # get in-sample intra-day VaR curve by assuming a point-wisely Gaussian distributed error term.
#' var_obj = var.forecast(yd, sigma_fit, error_fit, quantile_v=0.01, Method="normal")
#'
#' # the intra-day VaR curves.
#' var_obj$intraday_VaR
#' }
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Forecasting Value at Risk via intra-day return curves. International Journal of Forecasting. <doi:10.1016/j.ijforecast.2019.10.006>.\cr
var.forecast <- function(yd,sigma_pred,error_fit,quantile_v,Method){
  # a function to resample curve data with replacement.
  resmp <- function(data,boot_N,boot_in_N){
    bootres=matrix(NA,1,boot_N)
    for (j in 1:boot_N){
      bootres[,j]=data[sample(c(1:boot_in_N),1,replace=TRUE)]
    }
    return(bootres)}

  point_grid=nrow(yd)
  in_N=ncol(error_fit)
  switch(Method,
         normal = {
           error_var=qnorm(p=quantile_v,mean=0,sd=1)
           daily_var=(sigma_pred[point_grid])*error_var
           ea=yd[point_grid,]
           ealog=(ea<daily_var)
           daily_es=mean(ea[ealog])
           var_curves=sigma_pred*error_var
         },
         bootstrap = {
           bootN=1000
           error_var=matrix(NA,point_grid,1)
           for (j in 1:point_grid){
             bootres=resmp(error_fit[j,],bootN,in_N)
             error_var[j,1]=quantile(bootres,quantile_v)}
           daily_var=sigma_pred[point_grid]*error_var[point_grid,1]
           ea=yd[point_grid,]
           ealog=(ea<daily_var)
           daily_es=mean(ea[ealog])
           var_curves=sigma_pred*error_var
         },
         evt = {
           M=4 # this can be altered, here we stick on the implementation of Rice, Wirjanto, Zhao (2020).
           times=rep(0,point_grid)
           for(i in 1:point_grid){times[i]=i/point_grid}
           basis_obj=create.bspline.basis(rangeval=c(times[1],1),nbasis=32,norder=4)
           y_sq_fdsmooth=smooth.basis(argvals=times,y=error_fit,fdParobj=basis_obj)
           y_sq_fd=y_sq_fdsmooth$fd
           pca_obj=pca.fd(y_sq_fd,nharm=10,harmfdPar=fdPar(y_sq_fd), centerfns=FALSE)  #also try centerfns = TRUE
           eigen_functs=pca_obj$harmonics
           error_basis=matrix(0,nrow=point_grid,ncol=M)
           for(j in 1:M){error_basis[,j]=eval.fd(times,eigen_functs[j])}
           int_approx=function(x){temp_n=NROW(x)
           return((1/temp_n)*sum(x))}
           error_proj=matrix(NA,in_N,M)
           for(i in 1:in_N){
             for(j in 1:M){error_proj[i,j]=int_approx(error_fit[,i]*error_basis[,j])}}

           #### predict error_proj by using the generalised extreme value theory, the threshold parameter mythh is subjectively selected by graphical methods.
           mythh=rep(1,M)
           for (i in 1:M){
             mythh[i]=quantile(ecdf(error_proj[,i]),0.03)}
           gpd_sim=matrix(NA,1000,M)

           for (j in 1:M){
             gev=fitgpd(data=as.vector(error_proj[,j]),threshold=mythh[j],est="mle")
             scale_sigma=gev[[1]][1];shape_xi=gev[[1]][2]
             gpd_sim[,j]=rgpd(1000,loc=mythh[j],scale=scale_sigma,shape=shape_xi)
           }
           #### reverse-KL explansion: from the error_proj (from generalised pareto distribution) to infinite dimensional [0,1] space error_pred
           errorq_pred=matrix(0,point_grid,1000)
           ytq_es=matrix(NA,point_grid,1)
           for (i in 1:1000){
             for(j in 1:M){
               errorq_pred[,i]=errorq_pred[,i]+gpd_sim[i,j]*error_basis[,j] }}

           errorq_var=matrix(NA,point_grid,1)

           for(i in 1:point_grid){
             errorq_var[i,]=quantile(errorq_pred[i,],quantile_v)
           }

           daily_var=(sigma_pred[point_grid])*errorq_var[point_grid,]
           ea=yd[point_grid,]
           ealog=(ea<daily_var)
           daily_es=mean(ea[ealog])
           var_curves=sigma_pred*errorq_var
         },
         stop("Enter something to switch me!")
  )

  return(list(day_VaR = daily_var, day_ES = daily_es, intraday_VaR = var_curves))
}


#' Violation Process
#'
#' @description var.vio function returns a violation process for the intra-day VaR curves.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N functional curves.
#' @param var_curve A (grid_point) x (number of observations) matrix storing the forecasts of intra-day VaR curves.
#'
#' @return A violation process of the intra-day return curves based on the forecasts of intra-day VaR curves.
#' @export
#'
#' @details
#' Given the intra-day return curves \eqn{x_i(t)}, and the forecasts of intra-day VaR curves \eqn{\widehat{VaR}_i^\tau(t)} obtained from \code{\link{var.forecast}}, the violation process \eqn{Z_i^\tau(t)} can be defined as, \cr
#' \eqn{Z_i^\tau(t)=I(x_i(t)<\widehat{VaR}_i^\tau(t))}, for \eqn{1\leq i \leq N}, \eqn{t\in[0,1]}, and \eqn{\tau \in [0,1]},\cr
#' where \eqn{I(\cdot)} is an indicator function.
#'
#' @seealso \code{\link{var.forecast}} \code{\link{var.backtest}}
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Forecasting Value at Risk via intra-day return curves. International Journal of Forecasting. <doi:10.1016/j.ijforecast.2019.10.006>.\cr
var.vio <- function(yd,var_curve){
  if( ncol(yd)!=ncol(var_curve) ) warning('the number of observations must match for calculating the violation process!')
  point_grid=nrow(yd)
  n=ncol(yd)
  y_vio=matrix(0,point_grid,n)
  for (i in 1:n){
    ind=as.matrix(yd)[,i]<as.matrix(var_curve)[,i]
    y_vio[ind,i]=1
  }
  return(y_vio)
}



#' Backtest Intra-day VaR forecasts
#'
#' @description var.backtest function backtests the unbiasedness and the independence hypotheses for the intra-day VaR curve forecasts.
#'
#' @param vio A (grid_point) x (number of observations) matrix drawn from the violation process curves.
#' @param tau The nominal/true quantile of the VaR curves.
#' @param K The maximal lagged autocorrelation considered for the independence test. If it is missing, a default value "K=20" is used.
#'
#' @return
#' @return List of generated processes:
#' @return unbias: the p-value of the unbiasedness test;
#' @return independent: the p-value of the independence test.
#'
#' @export
#'
#' @details
#' Given the violation process \eqn{Z_i^\tau(t)} at the quantile \eqn{\tau}, the function computes P-values of two hypothesis tests:
#' (1) unbiasedness \eqn{H_0}: \eqn{E(Z_i^\tau(t)-\tau)=0}, for all \eqn{t \in[0,1]}, \eqn{1\leq i\leq N}. The test statistics,\cr
#' \eqn{T_N=N|| \bar{Z}(t)-\tau||^2} is employed, where \eqn{||\cdot ||} is the \eqn{L^2} norm, and \eqn{\bar{Z}(t)=1/N\sum_{i=1}^N Z_i(t)}.
#' (2) independence \eqn{H_0}: \eqn{Z_i^\tau(t)} is independent, a portmanteau test statistic in used (Kokoszka et al., 2017),\cr
#' \eqn{V_{N,K}=N\sum_{h=1}^K||\hat{\gamma}_{h,Z}||^2}, \cr
#' where \eqn{K} is a pre-set maximum lag length, and the autocovariance function \eqn{\hat{\gamma}_{h,Z}(t,s)=\frac{1}{N}\sum_{i=1}^{N-h}[Z_i(t)-\bar{Z}_i(t)][Z_{i+h}(s)-\bar{Z}(s)]}, for \eqn{||\cdot ||} is the \eqn{L^2} norm, and \eqn{\bar{Z}(t)=1/N\sum_{i=1}^N Z_i(t)}.
#'
#' @seealso \code{\link{var.forecast}} \code{\link{var.vio}}
#'
#' @examples
#' \dontrun{
#' # generate discrete evaluations of the FGARCH(1,1) process.
#' grid_point = 50; N = 200
#' yd = dgp.fgarch(grid_point, N, "garch")
#' yd = yd$garch_mat
#'
#' # extract data-driven basis functions through the truncated FPCA method.
#' basis_est = basis.est(yd, M=2, "tfpca")$basis
#'
#' # fit the curve data and the conditional volatility by using a FGARCH(1,1) model with M=1.
#' fd = fda::Data2fd(argvals=seq(0,1,len=grid_point),y=yd,fda::create.bspline.basis(nbasis=32))
#' garch11_est = est.fGarch(fd, basis_est[,1])
#' diag_garch = diagnostic.fGarch(garch11_est, basis_est[,1], yd)
#'
#' # get the in-sample fitting of conditional variance.
#' sigma_fit = diag_garch$sigma2[,1:N]
#' error_fit = diag_garch$eps
#'
#' # get in-sample intra-day VaR curve by assuming a point-wisely Gaussian distributed error term.
#' var_obj = var.forecast(yd, sigma_fit, error_fit, quantile_v=0.01, Method="normal")
#' intra_var = var_obj$intraday_VaR
#'
#' # get the violation curves.
#' intra_vio = var.vio(yd,intra_var)
#'
#' # backtesting the Unbiasedness Hypothesis for the violation curve.
#' pvalues = var.backtest(vio=intra_vio, tau=0.01, K=10)
#' pvalues$unbias
#' pvalues$independent
#' }
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Forecasting Value at Risk via intra-day return curves. International Journal of Forecasting. <doi:10.1016/j.ijforecast.2019.10.006>.\cr
#'
var.backtest <- function(vio, tau, K=NULL){
  #########################
  ### unbiasedness test ###
  #########################
  # the test T statistics
  H0_Tn<-function(z,tau){
    point_grid=nrow(z)
    n=ncol(z)
    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}
    Tn=n*int_approx(colMeans(z)-tau)^2
    return(Tn)
  }
  # get the critical values
  Tn_cv<-function(z,cv_N){
    N=ncol(z)
    point_grid=nrow(z)

    # autocovariance operator from -N to N
    cov_est<-function(x,h){
      N=ncol(x)
      c_est=x[,1:(N-h)]%*%t(x[,(h+1):N])/(N-h)
      return(c_est)
    }

    Gam=0
    Gam=Gam+cov_est(z,0)#long_run_cov2(z, C0 = 3, H = 3)#c

    eig=eigen(Gam)
    eigvals=eig$values
    eigvals=as.vector(eigvals)
    vect=eigvals/point_grid

    int_approx=function(x){
      temp_n=NROW(x)
      return((1/temp_n)*sum(x))}

    lim_sum=matrix(0,cv_N,1)

    lim_sum=0
    for (j in 1:length(vect)){
      lim_sum=lim_sum+vect[j]*rnorm(cv_N,mean=0,sd=1)}

    cv=quantile(lim_sum,probs=c(0.90,0.95,0.99))
    return(list(cv,lim_sum))
  }

  cv_N=5000 # sample size for approximating critical values
  Tn_stat=H0_Tn(vio,tau)
  limit=Tn_cv(vio,cv_N)
  emplim=limit[[2]]
  unbias_p = 1-ecdf(emplim)(Tn_stat)

  #########################
  ### independent test ###
  #########################
  if(is.null(K) == TRUE) {
    K=20
  }

  T_statistic <- function(fdata, lag_val){
    T = nrow(fdata)
    p = ncol(fdata)
    # calculate autocovariance function
    gamma_l <- function(fdata, lag){
      T = nrow(fdata)
      center_dat = t(scale(fdata, center = TRUE, scale = FALSE))

      gamma_lag_sum = 0
      for(ij in 1:(T - lag)){
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*% t(as.matrix(center_dat[,(ij + lag)])))
      }
      return(gamma_lag_sum/T)
    }
    gamma_hat = gamma_l(fdata = fdata, lag = lag_val)
    T_h = T * ((1/p)^2) * sum(gamma_hat^2)
    return(T_h)
  }

  # Portmanteau test statistics
  gaPort.stat <- function(H, datmat){
    vals = rep(0,H)
    for(j in 1:H){
      vals[j] = T_statistic(t(datmat), j)
    }
    return(sum(vals))
  }

  # useful inner functions
  eta <- function(i, j, datmat){
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = max(c(i,j))
    sum1 = array(0, c(p, p, p, p))
    for(k in (l+1):T){
      sum1 = sum1 + datmat[,k-i] %o% datmat[,k] %o% datmat[,k-j] %o% datmat[,k]
    }
    return(sum1/T)
  }

  etaM <- function(i, datmat){
    T = dim(datmat)[2]
    p = dim(datmat)[1]
    l = i
    sum1 = array(0, c(p,p))
    for(k in (l+1):T){
      sum1 = sum1 + (datmat[,k-i])^2 %o% (datmat[,k])^2
    }
    return(sum1/T)
  }

  etapopM <- function(H, datmat){
    etal = list()
    for(j in 1:H){
      etal[[j]] = list()
    }
    for(k in 1:H){
      etal[[k]] = etaM(k, datmat)
    }
    return(etal)
  }

  mean.W.2 <- function(x.e, H, datmat){
    sum1 = 0
    p = dim(datmat)[1]
    for(j in 1:H){
      sum1 = sum1 + sum(x.e[[j]])
    }
    return(sum1/p^2)
  }

  etapopNvarMC2 <- function(H, datmat, len1, len2){
    sum1 = 0

    rref = runif(len1, 0, 1)
    rref = sort(rref)
    rrefind = round(rref * dim(datmat)[1])
    rrefind[which(rrefind == 0)] = 1

    rrefG = c(0, rref)
    xd = diff(rrefG)
    gmat = xd %o% xd %o% xd %o% xd

    for(m in 1:H){
      x = eta(m, m, datmat[rrefind,])
      sum1 = sum1 + sum((x*x)*gmat)
    }

    rref2 = runif(len2, 0, 1)
    rref2 = sort(rref2)
    rrefind2 = round(rref2 * dim(datmat)[1])
    rrefind2[which(rrefind2 == 0)] = 1
    rrefG2 = c(0, rref2)
    xd2 = diff(rrefG2)
    gmat2 = xd2 %o% xd2 %o% xd2 %o% xd2

    if(H > 1){
      for(j in 1:(H-1)){
        for(k in (j+1):H){
          if((abs(k-j)) < 3){
            x = eta(j, k, datmat[rrefind,])
            #sum1=sum1+sum( (x.e[[j]][[k]])^2)
            sum1 = sum1 + 2 * sum((x*x)*gmat)
          }#end if

          if((abs(k-j)) >= 3){
            x = eta(j, k, datmat[rrefind2,])
            #sum1=sum1+sum( (x.e[[j]][[k]])^2)
            sum1 = sum1 + 2 * sum((x*x)*gmat2)
          }
        }
      }
    }
    return(2 * sum1)
  }

  vio = vio - rowMeans(vio)
  len1 = 20;len2 = 15; # parameters for Monte Carlo integral, can be altered.
  res = gaPort.stat(K, vio)
  x.e = etapopM(K, vio)
  res2 = mean.W.2(x.e, K, vio)
  res3 = etapopNvarMC2(K, vio, len1, len2)
  beta = res3/(2 * res2)
  nu = (2 * res2^2)/res3
  independent_p = 1 - pchisq(res/beta, df = nu)

  return(list(unbias = unbias_p, independent = independent_p))
}
