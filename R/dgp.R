#' Data Generating Process - Independent Process
#'
#' @description dgp.fiid function generates iid functional curve data following the Ornstein–Uhlenbeck process.
#'
#' @param grid_point The number of grid point in each curve observation.
#' @param N The sample size.
#'
#' @return A (grid_point) x (number of observations) matrix for iid sequences, where the finite realization of curves are stored in columns.
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @details
#' \eqn{x_i(t)=e^{-t/2}W_i(e^t)}, \eqn{t \in [0,1]},\cr
#' where \eqn{W_i(t)} is a standard Brownian Motion.
#'
#' @seealso \code{\link{dgp.farch}} \code{\link{dgp.fgarch}}
#' @examples
#' # generate discrete evaluations of 100 iid curves that each curve is realized on 50 grid points.
#' yd_iid = dgp.fiid(50,100)
#'
#' # smooth discrete data into functional curves.
#' fd = fda::Data2fd(argvals=seq(0,1,len = 50),y=yd_iid,fda::create.bspline.basis(nbasis = 32))
dgp.fiid <- function(grid_point, N){
  times=1:grid_point/grid_point
  # control the covariance structure
  comat=matrix(NA,grid_point,grid_point)
  for (i in 1:grid_point){
    comat[i,]=exp(1)^(-times[i]/2-times/2)*pmin(exp(1)^(times[i]),exp(1)^(times))
  }
  fiid=mvrnorm(n = N, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)

  return(t(fiid))
}




#' Data Generating Process - Functional ARCH Process
#'
#' @description dgp.farch function generates functional curve data following the functional ARCH(1) process.
#'
#' @param grid_point The number of grid point in each curve observation.
#' @param N The sample size.
#' @param beta_par The kernel coefficient function in the conditional volatility equation. If it is missing, the default function "\eqn{16 * t * (1-t) * s * (1-s)}" is used, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#'
#' @return List of generated processes:
#' @return arch_mat: FARCH sequences, where the finite realization of curves are stored in columns,
#' @return sigma_mat: Conditional volatility sequences,  where the finite realization of curves are stored in columns.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @details
#' \eqn{x_i(t)} is generated to follow an FARCH(1) process, such that,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds},\cr
#' where the innovation \eqn{\varepsilon_i(t)} follows an Ornstein–Uhlenbeck process \code{\link{dgp.fiid}}, and the kernel coefficients \eqn{\omega(t)=0.1t(1-t)}, \eqn{\alpha(t,s)=12t(1-t)s(1-s)}.
#'
#' @seealso \code{\link{dgp.fiid}} \code{\link{dgp.fgarch}}
#' @examples
#' # generate discrete evaluations of 100 FARCH curves that each curve is realized on 50 grid points.
#' yd = dgp.farch(50,100)
#' yd_arch = yd$arch_mat
#'
#' # smooth discrete data into functional curves.
#' fd = fda::Data2fd(argvals=seq(0,1,len = 50),y=yd_arch,fda::create.bspline.basis(nbasis = 32))
#' @references
#' Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the ARCH model. Econometric Theory, 29(2), 267-288.
#'
dgp.farch <- function(grid_point, N, beta_par=NULL){
  times = 1:grid_point/grid_point
  error_sim <- function(grid_point, N){
    ti=1:grid_point/grid_point
    comat=matrix(NA,grid_point,grid_point)
    for (i in 1:grid_point){
      comat[i,]=exp(1)^(-ti[i]/2-ti/2)*pmin(exp(1)^(ti[i]),exp(1)^(ti))
    }
    epsilon=mvrnorm(n = N, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)
    return(t(epsilon))
  }

  int_approx <- function(x){
    temp_n = nrow(x)
    return((1/(temp_n)) * sum(x))
  }

  no_proc = 1000 + N
  sim_sigma2_matrix = sim_arch_matrix = matrix(NA, grid_point, no_proc)
  error_matrix = error_sim(grid_point = grid_point, N = no_proc)

  if(is.null(beta_par) == TRUE) {
    beta_par = function(t,s){
      return(16 * t * (1-t) * s * (1-s))
    }
  }

  delta_par = function(t){
    return(0.1* t * (1-t))
  }
  #fill the initial values
  sim_sigma2_matrix[,1] = delta_par(times)
  sim_arch_matrix[,1] = delta_par(times) * error_matrix[,1]

  #fill in the rest of sigma2 and functional garch process matrices:
  for(j in 2:no_proc){
    #first fill in sigma2 column:
    for(i in 1:grid_point){
      vector_op = beta_par(times[i], times) * ((sim_arch_matrix[,(j-1)])^2)
      sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(as.matrix(vector_op))
    }
    #fill in garch process values for column:
    sim_arch_matrix[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
  }
  return(list(arch_mat = sim_arch_matrix[,1001:(1000 + N)],
              sigma_mat = sim_sigma2_matrix[,1001:(1000 + N)]))
}




#' Data Generating Process - Functional GARCH Process
#'
#' @description dgp.fgarch function generates an functional GARCH(1,1) process.
#'
#' @param grid_point The number of grid point in each curve observation.
#' @param N The sample size.
#' @param alpha_par The ARCH kernel coefficient function in the conditional volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}" is used, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#' @param beta_par The GARCH kernel coefficient function in the conditional volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}" is used, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#'
#'
#' @return List of generated processes:
#' @return garch_mat: FARCH sequences, where the finite realization of curves are stored in columns,
#' @return sigma_mat: Conditional volatility sequences,  where the finite realization of curves are stored in columns.
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#'
#' @details
#' \eqn{x_i(t)} is generated to follow an FGARCH(1,1) process, such that,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds+\int \beta(t,s) \sigma^2_{i-1}(s)ds},\cr
#' where the innovation \eqn{\varepsilon_i(t)} follows an Ornstein–Uhlenbeck process \code{\link{dgp.fiid}}, and the kernel coefficients \eqn{\omega(t)=0.1t(1-t)}, \eqn{\alpha(t,s)=\beta(t,s) = 12t(1-t)s(1-s)}.
#'
#' @seealso \code{\link{dgp.fiid}} \code{\link{dgp.farch}}
#'
#' @examples
#' # generate discrete evaluations of 100 FGARCH curves that each curve is realized on 50 grid points.
#' yd = dgp.fgarch(50,100)
#' yd_garch = yd$garch_mat
#'
#' # smooth discrete data into functional curves.
#' fd = fda::Data2fd(argvals=seq(0,1,len = 50),y=yd_garch,fda::create.bspline.basis(nbasis = 32))
#'
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis, 38(1), 3-21.
#'
dgp.fgarch <- function(grid_point, N, alpha_par=NULL, beta_par=NULL){
  times = 1:grid_point/grid_point
  error_sim <- function(grid_point, N){
    ti=1:grid_point/grid_point
    comat=matrix(NA,grid_point,grid_point)
    for (i in 1:grid_point){
      comat[i,]=exp(1)^(-ti[i]/2-ti/2)*pmin(exp(1)^(ti[i]),exp(1)^(ti))
    }
    epsilon=mvrnorm(n = N, mu = c(rep(0,grid_point)), Sigma = comat, empirical = FALSE)
    return(t(epsilon))
  }

  int_approx <- function(x){
    temp_n = nrow(x)
    return((1/(temp_n)) * sum(x))
  }

  no_proc = 1000 + N
  sim_sigma2_matrix = sim_garch_matrix = matrix(NA, grid_point, no_proc)
  error_matrix = error_sim(grid_point = grid_point, N = no_proc)


  if( is.null(alpha_par) == TRUE) {
    alpha_par = function(t,s){
      return(12 * t * (1-t) * s * (1-s))
    }
  }

  if( is.null(beta_par) == TRUE) {
    beta_par = function(t,s){
      return(12 * t * (1-t) * s * (1-s))
    }
  }

  delta_par = function(t){
    return(0.5* t * (1-t))
  }
  #fill the initial values
  sim_sigma2_matrix[,1] = delta_par(times)
  sim_garch_matrix[,1] = delta_par(times) * error_matrix[,1]

  #fill in the rest of sigma2 and functional garch process matrices:
  for(j in 2:no_proc){
    #first fill in sigma2 column:
    for(i in 1:grid_point){
      vector_for_alpha_op = alpha_par(times[i], times) * ((sim_garch_matrix[,(j-1)])^2)
      vector_for_beta_op = beta_par(times[i], times) * ((sim_sigma2_matrix[,j-1]))
      sim_sigma2_matrix[i,j] = delta_par(times[i]) + int_approx(as.matrix(vector_for_alpha_op)) + int_approx(as.matrix(vector_for_beta_op))
    }
    #fill in garch process values for column:
    sim_garch_matrix[,j] = sqrt(sim_sigma2_matrix[,j]) * error_matrix[,j]
  }
  return(list(garch_mat = sim_garch_matrix[,1001:(1000 + N)],
              sigma_mat = sim_sigma2_matrix[,1001:(1000 + N)]))
}







#' @title Sample Data for Analysis
#'
#' @description A dataset of S&P 500 intra-day price at one-minute frequence from 9:30 a.m. to 4:00 p.m.. The sample ranges from 02/January/2015 to 29/March/2018, and trading days with missing intra-day observations have been eliminated from the sample.
#'
#' @format a 390 by 788 matrix that each column containing the raw intra-day price discretely observed at daily basis.
#'
#' @seealso \code{\link{intra.return}}
#' @examples
#' yd = as.matrix(sample_data)
#' grid_point = dim(yd)[1]
#' N = dim(yd)[2]
#'
#' # smooth data into functional objects by using B-spline bases.
#' fd = fda::Data2fd(y=yd, argvals=seq(0,1,len = grid_point), fda::create.bspline.basis(nbasis = 32))
"sample_data"
