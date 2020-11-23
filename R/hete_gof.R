#' Test for Conditional Heteroscedasticity of Functional Data
#'
#' @description fun_hetero applies a test of the null hypothesis that the objective functional curve data is not conditionally heteroscedastic, small p-values suggest that the curves exhibit conditional heteroscedasticity.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N intra-day return curves.
#' @param K The lag autocorrelation coefficients served for the conditional heteroscedasticity test. If it is missing, a default value "K=20" is used.
#' @param stat_Method A string to indicate which test will be implemented: "norm" - \eqn{V_{N,K}}; "functional" - \eqn{M_{N,K}}. Or if it is missing, the "functional"-type method is used.
#' @param pplot an optional argument to compute and plot the P-values (either "norm" or "functional") as a function of K, K=1,2,...,20. If pplot=1, the p-values will be computed and figured; if it is missing, this step will be skipped.
#'
#' @return List of objects:
#' @return stats: the test statistic.
#' @return p_value: the P_value.
#'
#'
#' @export
#' @importFrom graphics abline plot
#'
#' @import stats
#'
#' @details
#' Given the objective curve data \eqn{x_i(t)}, for \eqn{1\leq i \leq N}, \eqn{t\in[0,1]}, the test aims at distinguishing the hypotheses:\cr
#' \eqn{H_0}: the sequence \eqn{x_i(t)} is IID; \cr
#' \eqn{H_1}: the sequence \eqn{x_i(t)} is conditionally heteroscedastic. Two portmanteau type statistics are applied: \cr
#' the norm-based statistic: \eqn{V_{N,K}=N\sum_{h=1}^K\hat{\rho}_h^2}, where \eqn{\hat{\rho}_h} is the sample autocorrelation of the time series \eqn{||x_1||^2,\dots,||x_N||^2}, and \eqn{K} is a pre-set maximum lag length.\cr
#' the fully functional statistic \eqn{M_{N,K}=N\sum_{h=1}^K||\hat{\gamma}_h||^2}, where the autocovariance kernel \eqn{\hat{\gamma}_{h}(t,s)=\frac{1}{N}\sum_{i=1}^{N-h}[x_i^2(t)-\bar{x}^2(t)][x^2_{i+h}(s)-\bar{x}(s)]}, for \eqn{||\cdot ||} is the \eqn{L^2} norm, and \eqn{\bar{x}(t)=1/N\sum_{i=1}^N x^2_i(t)}.
#'
#' @seealso \code{\link{sample_data}}
#'
#' @examples
#' # generate discrete evaluations of the iid curves under the null hypothesis.
#' yd = dgp.fiid(50, 100)
#'
#' # test the conditional heteroscedasticity.
#' fun_hetero(yd, K=5, "functional")
#'
#' @references
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Tests for conditional heteroscedasticity of functional data. Journal of Time Series Analysis.
fun_hetero <- function (yd, K=NULL, stat_Method, pplot=NULL){

  if(is.null(K) == TRUE) {
    K=20
  }
  ## functions to compute the norm statistics
  normx2 <- function(x){ # compute the L2 norm
    N=ncol(x)
    grid_point=nrow(x)
    y1=matrix(NA,1,N)
    x=x-rowMeans(x)
    for(i in 1:N){
      y1[,i]=(1/grid_point)*sum(x[,i]^2)
    }
    y1=y1-mean(y1)
    return(y1)
  }

  gamma1 <- function(x,lag_val){ # compute the covariance of normx2 objective series
    N=ncol(x)
    x=x-mean(x)
    gamma_lag_sum=0
    for(j in 1:(N-lag_val)){
      gamma_lag_sum=gamma_lag_sum+x[j]*x[(j+lag_val)]
    }
    gamma1_h=gamma_lag_sum/N
    return(gamma1_h)
  }

  Port1.stat <- function(dy,Hf){ # the norm statistics
    vals=rep(0,Hf)
    N=ncol(dy)
    y2=normx2(dy)
    sigma2=gamma1(y2,0)
    for(h in 1:Hf){
      gam1hat=gamma1(y2,h)
      vals[h]=N*(gam1hat/sigma2)^2 # the square of autocorrelation function
    }
    stat=sum(vals)
    pv=(1-pchisq(stat,df=Hf))

    return(list(stats = stat, p_value = pv))
  }

  ## functions to compute the fully functional statistics
  func1 <- function(x){
    N=ncol(x)
    grid_point=nrow(x)
    x=x-rowMeans(x)
    y2=matrix(NA,grid_point,N)
    for(i in 1:N){
      y2[,i]=x[,i]^2
    }
    y2=y2-rowMeans(y2)
    return(y2)
  }

  gamma2 <- function(x, lag_val){ # compute the covariance matrix
    N=ncol(x)
    gam2_sum=0
    for (j in 1:(N-lag_val)){
      gam2_sum=gam2_sum+x[,j]%*%t(x[,(j+lag_val)])
    }
    gam2=gam2_sum/N
    return(gam2)
  }

  Port2.stat <- function(dy,Hf){ # the fully functional statistics
    N=ncol(dy)
    grid_point=nrow(dy)
    vals=rep(0,Hf)
    y2=func1(dy)
    for(h in 1:Hf){
      gam2hat=gamma2(y2,h)
      vals[h]=N*((1/grid_point)^2)*sum(gam2hat^2)}
    stat=sum(vals)
    return(stat)
  }

  # functions to compute the p_value of the fully functional statistics

  etaM <- function(dy){
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    dy=func1(dy)
    sum1=0
    for(k in 1:N){
      sum1=sum1+(dy[,k])^2
    }
    return(sum1/N)
  }

  mean.W.2 <- function(dy, Hf){ # mean function
    sum1=0
    grid_point=dim(dy)[1]
    etal=etaM(dy)
    sum1=(sum(etal)/grid_point)^2
    sum1=sum1*Hf
    return(sum1)
  }

  etaM_cov <- function(dy){
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    dy=func1(dy)
    sum1=0
    for(k in 1:N){
      sum1=sum1+(dy[,k])%*%t(dy[,k])
    }
    return(sum1/N)
  }

  etapopNvarMC2 <- function(dy, Hf){ # covariance function
    N=dim(dy)[2]
    grid_point=dim(dy)[1]
    x=etaM_cov(dy)^2
    sum1=sum(x)/(grid_point^2)
    sum1=2*Hf*(sum1^2)
    return(sum1)
  }

  Port.test.2 <- function(dy, Hf){ # Welch–SaNerthwaite approximation
    stat=Port2.stat(dy, Hf)
    res2=mean.W.2(dy, Hf)
    res3=etapopNvarMC2(dy, Hf)
    beta=res3/(2*res2)
    nu=(2*res2^2)/res3
    pv = 1-pchisq(stat/beta, df=nu)
    return(list(stats = stat, p_value = pv))
  }

  if(is.null(stat_Method) == TRUE) {
    stat_Method="functional"
  }

  switch(stat_Method,
         norm = {
           stat_p=Port1.stat(yd,K)
         },
         functional = {
           stat_p=Port.test.2(yd,K)
         },
         stop("Enter something to switch me!"))

  if(is.null(pplot) == TRUE) {
    stat_p=stat_p
  }
  else if(pplot==1){
    kseq=1:20
    pmat=c(rep(1,20))
    for (i in 1:20){
      switch(stat_Method,
             norm = {
               pmat[i]=Port1.stat(yd,kseq[i])[2]
             },
             functional = {
               pmat[i]=Port.test.2(yd,kseq[i])[2]
             },
             stop("Enter something to switch me!"))
    }
    x<-1:20
    plot(x,pmat, col="black",type="l",ylim=c(0,1.1),xlab="Lag",ylab="P-Values",main="P-values of conditional heteroscedasticity test")
    abline(h=0.1, col="red")
  }

  return(stat_p)
}





#' Goodness-of-fit Test for Functional ARCH/GARCH Model
#'
#' @description gof.fgarch function approximates the P-value of the \eqn{M_{N,K}} statistics accounting for the effect of functional GARCH parameter estimation.
#'
#' @param yd A (grid_point) x (number of observations) matrix drawn from N discrete evaluation curves.
#' @param basis The M-dimensional basis functions.
#' @param model A string to indicate which model will be estimated: "arch" - FARCH(1); "garch" - FGARCH(1,1).
#' @param K The statistic will be based on H lag autocorrelation coefficients. If it is missing, a default value "K=20" is used.
#' @param pplot An optional argument to compute and plot the P-values as a function of K, K=1,2,...,20. If pplot=1, the p-values will be computed and figured; if it is missing, this step will be skipped.
#' @param max_eval The maximum number of evaluations of the optimization function.
#'
#' @return The P_value of the \eqn{M_{N,K}} statistic.
#'
#' @export
#' @importFrom graphics abline plot
#' @importFrom nloptr cobyla
#'
#' @details
#' The test statistic used is as the same as the \eqn{M_{N,K}} statistic in \code{\link{fun_hetero}}. However, the asymptotic distribution is adjusted to account for the estimation effect, because the model residual depends on the joint asymptotics of the innovation process and the estimated parameters. We assume that the kernel parameters are consistently estimated by the Least Squares method proposed in Aue et al. (2017). Then, the asymptotic distribution of the statistic \eqn{M_{N,K}} is given in Theorem 3.1 in Rice et al. (2020).
#'
#' @seealso \code{\link{basis.tfpca}} \code{\link{basis.score}} \code{\link{est.fGarch}} \code{\link{diagnostic.fGarch}}
#'
#' @examples
#' # generate discrete evaluations of the FARCH process.
#' grid_point=50; N=200
#' yd = dgp.fgarch(grid_point, N, "garch")
#' yd = yd$garch_mat
#' ba = basis.tfpca(yd, M=2)
#' basis_est = ba$basis
#'
#' # test the adequacy of the FGARCH(1,1) model.
#' gof.fgarch(yd, basis_est[,1], "arch", K=5)
#'
#' @references
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis, 38(1), 3-21.\cr
#' Rice, G., Wirjanto, T., Zhao, Y. (2020). Tests for conditional heteroscedasticity of functional data. Journal of Time Series Analysis.
#'
gof.fgarch <- function (yd, basis, model, K=NULL, pplot=NULL, max_eval=10000){

  basis = as.matrix(basis)
  M = ncol(basis)
  sample_size = ncol(yd)
  grid_point = nrow(yd)
  int_approx=function(x){
    temp_n=NROW(x)
    return((1/temp_n)*sum(x))}
  y2m=yd*yd
  y_vec=matrix(0,nrow=M,ncol=sample_size)
  for(i in 1:sample_size){
    for(j in 1:M){y_vec[j,i]=int_approx(y2m[,i]*basis[,j])}}

  switch(model,
         arch = {
           get_theta=function(theta,M){
             #first M entries are d
             d=theta[1:M]
             #Next p*M^2 entries are A, A_inds are indices of A matrix params
             num_mat_params=(M^2-M*(M-1)/2)
             A_inds=M+(1:(num_mat_params))
             theta_A=theta[A_inds]

             curr_A_vals=theta_A[1:num_mat_params]
             A=matrix(0,ncol=M,nrow=M)
             diag(A)=curr_A_vals[1:M]
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             A=t(A)
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             return(list(ds=d,As=A))
           }

           function_to_minimize2=function(x,data=y_vec){
             sample_size=ncol(data)
             M=nrow(data)
             pam = get_theta(x,M)
             d=pam$d
             A=pam$As
             sigma_sq_proj_coefs=matrix(0,M,sample_size)
             sigma_sq_proj_coefs[,1]=d

             for(i in 2:(sample_size)){
               first_part=d
               second_part=A%*%data[,i-1]
               sigma_sq_proj_coefs[,i]=first_part+second_part
             }

             s_theta=0
             for (i in 2:sample_size) {s_theta=s_theta+t(data[,i]-sigma_sq_proj_coefs[,i])%*%(data[,i]-sigma_sq_proj_coefs[,i])}
             r=sum(s_theta)
             return(r)
           }

           #nonlinear inequality constraints on parameters
           eval_g0<-function(x,Mn=M){
             conpam=get_theta(x,Mn)
             A = conpam$As
             upb = 1-(10^-20)
             h <- -(sum(A)-upb)
             return(h)
           }
           num_params=M^2-M*(M-1)/2+M
           stav=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))
           ress=nloptr::cobyla(x0 = stav,fn = function_to_minimize2, lower=c(rep(10^-20,num_params)),
                               upper=c(rep(1,num_params)),
                               hin = eval_g0, nl.info = FALSE,
                               control = list(maxeval=max_eval))

           para=as.numeric(ress$par)
           pam_hat = get_theta(para,M)
           d_hat=pam_hat$d
           A_hat=pam_hat$As

           d_M=0
           for (i in 1:M){
             d_M=d_M+d_hat[i]*basis[,i]}

           alpha_M=matrix(0,grid_point,grid_point)
           for(i in 1:grid_point){
             for(j in 1:grid_point){
               temp_alpha_m=0
               for(t in 1:M){
                 for(s in 1:M){temp_alpha_m=temp_alpha_m+A_hat[t,s]*basis[i,t]*basis[j,s]                  }
               }
               alpha_M[i,j]=temp_alpha_m
             }}

           sigma_fit=matrix(1,grid_point,(sample_size+1))
           sigma_fit[,1]=d_M
           for(j in 2:(sample_size+1)){
             #first fill in sigma2 column:
             for(i in 1:grid_point){
               fit_alpha_op = alpha_M[i,] * ((yd[,(j-1)])^2)
               sigma_fit[i,j] = d_M[i] + int_approx(fit_alpha_op)
             }      }

           error_fit=yd/sqrt(abs(sigma_fit[,1:sample_size]))

           sigma_2_proj_coefs=matrix(0,M,sample_size)
           sigma_2_proj_coefs[,1]=d_hat
           for(i in 2:(sample_size)){
             first_part=d_hat
             second_part=A_hat%*%y_vec[,i-1]
             sigma_2_proj_coefs[,i]=first_part+second_part
           }
         },
         garch = {
           get_theta=function(theta,M){
             #first M entries are d
             d=theta[1:M]

             #Next p*M^2 entries are A, A_inds are indices of A matrix params
             num_mat_params=(M^2-M*(M-1)/2)
             A_inds=M+(1:(num_mat_params))
             B_inds=sfsmisc::last(A_inds)+(1:(num_mat_params))

             theta_A=theta[A_inds]
             theta_B=theta[B_inds]

             curr_A_vals=theta_A[1:num_mat_params]
             A=matrix(0,ncol=M,nrow=M)
             diag(A)=curr_A_vals[1:M]
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]
             A=t(A)
             A[lower.tri(A)]=curr_A_vals[(M+1):length(curr_A_vals)]

             #see above comment
             curr_B_vals=theta_B[1:num_mat_params]
             B=matrix(0,ncol=M,nrow=M)
             diag(B)=curr_B_vals[1:M]
             B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]
             B=t(B)
             B[lower.tri(B)]=curr_B_vals[(M+1):length(curr_B_vals)]

             return(list(ds=d,As=A,Bs=B))
           }

           function_to_minimize2=function(x,data=y_vec){
             sample_size=ncol(data)
             M=nrow(data)
             pam = get_theta(x,M)
             d=pam$d
             A=pam$As
             B=pam$Bs
             sigma_sq_proj_coefs=matrix(0,M,sample_size)
             sigma_sq_proj_coefs[,1]=d

             for(i in 2:(sample_size)){
               first_part=d
               second_part=0
               for (l in 1:(i-1)){
                 first_part=first_part+B^(l-1)%*%d
                 second_part=second_part+B^(l-1)%*%(A%*%data[,i-l])}
               sigma_sq_proj_coefs[,i]=first_part+second_part
             }
             s_theta=0
             for (i in 2:sample_size) {s_theta=s_theta+t(data[,i]-sigma_sq_proj_coefs[,i])%*%(data[,i]-sigma_sq_proj_coefs[,i])}
             r=sum(s_theta)
             return(r)
           }

           #nonlinear inequality constraints on parameters
           eval_g0<-function(x,Mn=M){
             h <- numeric(2)
             upb = 1-(10^-20)
             conpam=get_theta(x,M)
             A=conpam$As
             B=conpam$Bs
             h[1] <- -(sum(A)-upb)
             h[2] <- -(sum(B)-upb)
             return(h)
           }
           num_params=2*(M^2-M*(M-1)/2)+M
           stav=c(runif(M,10^-10,(1-10^-10)),runif(num_params-M,0,1))
           ress=nloptr::cobyla(x0 = stav,fn = function_to_minimize2, lower=c(rep(10^-20,num_params)),
                               upper=c(rep(1,num_params)),
                               hin = eval_g0, nl.info = FALSE,
                               control = list(maxeval=max_eval))

           para=as.numeric(ress$par)
           pam_hat = get_theta(para,M)
           d_hat=pam_hat$d
           A_hat=pam_hat$As
           B_hat=pam_hat$Bs

           d_M=0
           for (i in 1:M){
             d_M=d_M+d_hat[i]*basis[,i]}

           alpha_M=matrix(0,grid_point,grid_point)
           beta_M=matrix(0,grid_point,grid_point)
           for(i in 1:grid_point){
             for(j in 1:grid_point){
               temp_alpha_m=0
               temp_beta_m=0
               for(t in 1:M){
                 for(s in 1:M){
                   temp_alpha_m=temp_alpha_m+A_hat[t,s]*basis[i,t]*basis[j,s]
                   temp_beta_m=temp_beta_m+B_hat[t,s]*basis[i,t]*basis[j,s]
                 }
               }
               alpha_M[i,j]=temp_alpha_m
               beta_M[i,j]=temp_beta_m
             }}

           sigma_fit=matrix(1,grid_point,(sample_size+1))
           sigma_fit[,1]=d_M
           for(j in 2:(sample_size+1)){
             #first fill in sigma2 column:
             for(i in 1:grid_point){
               fit_alpha_op = alpha_M[i,] * ((yd[,(j-1)])^2)
               fit_beta_op = beta_M[i,] * ((sigma_fit[,j-1]))
               sigma_fit[i,j] = d_M[i] + int_approx(fit_alpha_op) + int_approx(fit_beta_op)
             }     }

           error_fit=yd/sqrt(abs(sigma_fit[,1:sample_size]))

           sigma_2_proj_coefs=matrix(0,M,sample_size)
           sigma_2_proj_coefs[,1]=d_hat
           for(i in 2:(sample_size)){
             first_part=d_hat
             second_part=A_hat%*%y_vec[,i-1]+B_hat%*%sigma_2_proj_coefs[,i-1]
             sigma_2_proj_coefs[,i]=first_part+second_part
           }
         },
         stop("Enter something to switch me!"))

  #### goodness-of-fit test

  gof_pvalue<-function(error_fit,sigma_fit,ortho_basis_matrix,y_2_proj_coefs,sigma_2_proj_coefs,para,K){

    ## functions to calculate the fully functional statistics.
    int_approx_tensor<-function(x){# x is a 4-dimensional tensor
      dt=length(dim(x))
      temp_n=nrow(x)
      return(sum(x)/(temp_n^dt))}

    func1 <- function(x){
      N=ncol(x)
      grid_point=nrow(x)
      x=x-rowMeans(x)
      y2=matrix(NA,grid_point,N)
      for(i in 1:N){
        y2[,i]=x[,i]^2
      }
      y2=y2-rowMeans(y2)
      return(y2)
    }

    gamma2 <- function(x, lag_val){ # compute the covariance matrix
      N=ncol(x)
      gam2_sum=0
      for (j in 1:(N-lag_val)){
        gam2_sum=gam2_sum+x[,j]%*%t(x[,(j+lag_val)])
      }
      gam2=gam2_sum/N
      return(gam2)
    }

    Port2.stat <- function(dy,Hf){ # the fully functional statistics
      N=ncol(dy)
      grid_point=nrow(dy)
      vals=rep(0,Hf)
      y2=func1(dy)
      for(h in 1:Hf){
        gam2hat=gamma2(y2,h)
        vals[h]=N*((1/grid_point)^2)*sum(gam2hat^2)}
      stat=sum(vals)
      return(stat)
    }

    ## functions to calculate covariance operator

    func1 <- function(x){
      N=ncol(x)
      grid_point=nrow(x)
      x=x-rowMeans(x)
      y2=matrix(NA,grid_point,N)
      for(i in 1:N){
        y2[,i]=x[,i]^2
      }
      y2=y2-rowMeans(y2)
      return(y2)
    }

    H0_f<-function(y_2_proj_coefs1){
      # H0 is a M x (M+2M^2) matrix
      M=nrow(y_2_proj_coefs1)
      n=ncol(y_2_proj_coefs1)

      ded=diag(M)
      dey=matrix(0, M, M^2)
      h0_sum=0
      #des=matrix(0, M, M^2)########  for garch
      for (i in 2:n){
        if (M>1){
          for (j in 1:M){
            dey[j,(((j-1)*M+1):(j*M))]=y_2_proj_coefs1[,(i-1)]}}
        else {dey[1,1]=y_2_proj_coefs1[,(i-1)]}
        h0_c=cbind(ded,dey)
        h0_sum=h0_sum+h0_c}

      return(h0_sum/(n-1))
      #return(rbind(ded,dey,des))  #  for garch
    }


    Q0_f<-function(y_2_proj_coefs1){
      # Q0 is a (M+2M^2) x (M+2M^2) matrix
      M=nrow(y_2_proj_coefs1)
      n=ncol(y_2_proj_coefs1)

      dey=matrix(0, M, M^2)
      #des=matrix(0, M, M^2)########  for garch
      q0_sum=0
      for (i in 2:n){
        if (M>1){
          for (j in 1:M){
            dey[j,(((j-1)*M+1):(j*M))]=y_2_proj_coefs1[,(i-1)]}}
        else {dey[1,1]=y_2_proj_coefs1[,(i-1)]}
        ded=diag(M)
        q0_c=cbind(ded,dey)
        q0_e=t(q0_c)%*%q0_c
        q0_sum=q0_sum+q0_e
      }
      return(q0_sum/(n-1))
    }


    J0_f<-function(y_2_proj_coefs1,sigma_2_proj_coefs){
      # J0 is a M x M matrix
      sumdys=0
      for (i in 1: ncol(y_2_proj_coefs1)){
        dys=y_2_proj_coefs1[,i]-sigma_2_proj_coefs[,i]
        sumdys=sumdys+dys%*%t(dys)
      }
      return(sumdys/ncol(y_2_proj_coefs1))
    }

    covariance_op_new<-function(error_fit,sigma_fit,ortho_basis_matrix,y_2_proj_coefs1,sigma_2_proj_coefs,para,h,g){
      ortho_basis_matrix=ortho_basis_matrix/(sum(ortho_basis_matrix)/nrow(error_fit))
      TT = dim(error_fit)[2]
      p = dim(error_fit)[1]

      M = dim(y_2_proj_coefs1)[1]
      sum1 = 0
      sum1_m = 0
      error_fit = func1(error_fit)
      for (k in 1:TT){
        sum1 = sum1 + (error_fit[,k]) %*% t(error_fit[,k])
        sum1_m = sum1_m + (error_fit[,k])^2
      }

      A = (sum1/TT)%o%(sum1/TT)
      cov_meanA = (sum1_m/TT)

      H0_est = H0_f(y_2_proj_coefs1)
      Q0_est = Q0_f(y_2_proj_coefs1)
      J0_est = J0_f(y_2_proj_coefs1,sigma_2_proj_coefs)
      innermatrix = solve(Q0_est)%*%t(H0_est)%*%J0_est%*%H0_est%*%solve(Q0_est)

      # kronecker %x%
      #sum2 = array(0, c(p, p, p, p))
      rep_permu<-function(M){
        permu=matrix(NA,M^2,2)
        permu[,1]=c(rep(c(1:M),each=M))
        permu[,2]=c(rep(c(1:M),M))
        return(permu)
      }
      permu = rep_permu(M)

      d_hat = para[1:M]
      A_hat = matrix(para[(M+1):(M+M^2)],nrow=M,ncol=M)

      G_g = array(0, c((M+M^2),p,p))

      for (k in 2:(TT-g)){
        p2_g=matrix(0,(M+M^2),p)
        parti_g=matrix(0,(M+M^2),1)
        for (m in 1:(M+M^2)){
          if ( m <= M ){
            vecm1=c(replicate(M,0))
            vecm1[m]=1
            partd = sum(vecm1)
            parti_g[m,1]=partd
          }
          else{
            vecm2=matrix(0,M,M)
            vecm2[permu[m-M,1],permu[m-M,2]]=1
            parta = sum(vecm2%*%y_2_proj_coefs1[,k+g-1])
            parti_g[m,1]= parta
          }
        }
        for (j in 1:(M+M^2)){
          if (j<=M){p2_g[j,]=ortho_basis_matrix[,j]*parti_g[j,1]}#(j-(ceiling(j/M)-1)*M)
          else{
            p2_g[j,]=(rowSums(ortho_basis_matrix[,(permu[j-M,1])]%o%ortho_basis_matrix[,(permu[j-M,2])])/p)*parti_g[j,1]}
        }
        G_g = G_g + ((t(replicate((M+M^2),(1/sigma_fit[,k+g])))*p2_g)%o%(error_fit[,k]))
      }
      G_g = G_g/(TT-g-1)

      G_h = array(0, c((M+M^2),p,p))

      for (k in 2:(TT-h)){
        p2_h=matrix(0,(M+M^2),p)
        parti_h=matrix(0,(M+M^2),1)
        for (m in 1:(M+M^2)){
          if ( m <= M ){
            vecm1=c(replicate(M,0))
            vecm1[m]=1
            partd = sum(as.matrix(vecm1))
            parti_h[m,1]=partd
          }
          else{
            vecm2=matrix(0,M,M)
            vecm2[permu[m-M,1],permu[m-M,2]]=1
            parta = sum(vecm2%*%y_2_proj_coefs1[,k+h-1])
            parti_h[m,1]= parta
          }

        }
        for (j in 1:(M+M^2)){
          if (j<=M){p2_h[j,]=ortho_basis_matrix[,j]*parti_h[j,1]}#(j-(ceiling(j/M)-1)*M)
          else{
            p2_h[j,]=(rowSums(ortho_basis_matrix[,(permu[j-M,1])]%o%ortho_basis_matrix[,(permu[j-M,2])])/p)*parti_h[j,1]}
        }
        G_h=G_h + ((t(replicate((M+M^2),(1/sigma_fit[,k+h])))*p2_h)%o%(error_fit[,k]))
      }
      G_h = G_h/(TT-h-1)

      arr_3_prod<-function(mat1,Lmat,mat2){
        p=dim(mat1)[2]
        L=dim(Lmat)[1]
        mat1t=aperm(mat1)
        x=array(0,c(p,p,p,p))
        for (i in 1:p){
          for (j in 1:p){
            x[i,j,,] = x[i,j,,] + mat1t[i,,]%*%Lmat%*%mat2[,,j]
          }
        }
        return(x)
      }
      B = arr_3_prod(G_h,innermatrix,G_g)
      rm(parti_h)
      rm(parti_g)
      rm(p2_h)
      rm(p2_g)

      sum3 = array(0, c(p, p, p, p))

      arr_prod<-function(mat1,Lmat){
        p=dim(mat1)[2]
        L=dim(Lmat)[1]
        mat1t=aperm(mat1)
        x=matrix(0,p,p)
        for (i in 1:p){
          x[i,] = x[i,] + mat1t[i,,]%*%Lmat
        }
        return(x)
      }

      for (k in 2:(TT-g)){
        smtx=matrix(0,M,M^2)
        for (m in 1:M^2){
          vecm2=matrix(0,M,M)
          vecm2[permu[m,1],permu[m,2]]=1
          parta = sum(vecm2%*%y_2_proj_coefs1[,k-1])
          smtx[permu[m,1],m] = parta}#t(y_2_proj_coefs1[,(k-1)])}
        partii=cbind(diag(M),smtx)
        innercov2 = as.matrix((solve(Q0_est) %*% t(partii) %*% (y_2_proj_coefs1[,k]-sigma_2_proj_coefs[,k])),(M+M^2),1)
        EG_g=arr_prod(G_g,innercov2)
        sum3 = sum3 + (error_fit[,k]) %o% (error_fit[,k+h]) %o% EG_g
      }
      C=sum3/(TT-g-1)
      D=array(0, c(p, p, p, p))
      for (i in 1:p){
        for (j in 1:p){
          D[i,j,,]=C[,,i,j]
        }
      }

      cov_est=B+C+D
      cov_mean=matrix(NA,p,p)
      for (i in 1:p){
        for (j in 1:p){
          cov_mean[i,j]=cov_est[i,j,i,j]
        }
      }
      return(list(A,cov_est,cov_mean,cov_meanA))
    }


    ## compute the fully functional statistics and p_values

    stat = Port2.stat(error_fit,K)
    mu_cov=0
    sigma2_cov=0

    len=20 # parameters in Monte Carlo integration
    len2=15
    grid_point=dim(error_fit)[1]

    rrefind <- floor(grid_point * sort(runif(len, 0, 1)))
    rrefind[which(rrefind == 0)] = 1
    r_samp <- c(0,rrefind) # c(0, rrefind) we define v_0 = 0 for computation of the product in rTrap
    xd <- diff(r_samp / grid_point)
    gmat_A <- xd %o% xd
    gmat <- xd %o% xd %o% xd %o% xd


    rrefind2 <- floor(grid_point * sort(runif(len2, 0, 1)))
    rrefind2[which(rrefind2 == 0)] = 1
    r_samp2 <- c(0,rrefind2) # we define v_0 = 0 for computation of the product in rTrap
    xd2 <- diff(r_samp2 / grid_point)
    gmat2_A <- xd2 %o% xd2
    gmat2 <- xd2 %o% xd2 %o% xd2 %o% xd2

    mu_cov=0
    sigma2_cov=0

    if (K>1){ # for K is greater than 1
      for (hh in 1:K){
        mean_A=array(0, c(len))
        mean_mat=array(0, c(len, len))
        cov_mat_A=array(0, c(len, len, len, len))
        cov_mat=array(0, c(len, len, len, len))

        covop_trace=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,hh,hh)
        mean_A[2:len] = covop_trace[[4]][-len]
        mean_mat[2:len,2:len] = covop_trace[[3]][-len,-len]
        mu_cov=mu_cov+sum((1/2) * (covop_trace[[4]]+mean_A)*xd)^2
        mu_cov=mu_cov+sum((1/2) * (covop_trace[[3]]+mean_mat)*gmat_A)^2

        cov_mat_A[2:len,2:len,2:len,2:len] = covop_trace[[1]][-len,-len,-len,-len]
        cov_mat[2:len,2:len,2:len,2:len] = covop_trace[[2]][-len,-len,-len,-len]
        sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_trace[[1]]^2+cov_mat_A^2)*gmat)
        sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_trace[[2]]^2+cov_mat^2)*gmat)
      }

      for (h1 in 1:K){
        for (h2 in (h1+1):K){
          if((abs(h2-h1)) < 3){
            cov_mat=array(0, c(len, len, len, len))

            covop_est=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,h1,h2)
            cov_mat[2:len,2:len,2:len,2:len] = covop_est[[2]][-len,-len,-len,-len]
            sigma2_cov=sigma2_cov+2*2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat)
          }

          if((abs(h2-h1)) >= 3){
            cov_mat=array(0, c(len2, len2, len2, len2))
            covop_est=covariance_op_new(error_fit[rrefind2,],sigma_fit[rrefind2,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind2,],length(rrefind2),1),y_2_proj_coefs,sigma_2_proj_coefs,para,h1,h2)
            cov_mat[2:len2,2:len2,2:len2,2:len2] = covop_est[[2]][-len2,-len2,-len2,-len2]
            sigma2_cov=sigma2_cov+2*2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat2)
          }
        }
      }
    }
    else{ # for K=1
      mean_A=array(0, c(len))
      mean_mat=array(0, c(len, len))
      cov_mat_A=array(0, c(len, len, len, len))
      cov_mat=array(0, c(len, len, len, len))

      covop_est=covariance_op_new(error_fit[rrefind,],sigma_fit[rrefind,],as.matrix(as.matrix(ortho_basis_matrix)[rrefind,],length(rrefind),1),y_2_proj_coefs,sigma_2_proj_coefs,para,1,1)

      mean_A[2:len] = covop_est[[4]][-len]
      mean_mat[2:len,2:len] = covop_est[[3]][-len,-len]

      mu_cov=mu_cov+sum((1/2) * (covop_est[[4]]+mean_A)*xd)^2
      mu_cov=mu_cov+sum((1/2) * (covop_est[[3]]+mean_mat)*gmat_A)^2

      cov_mat_A[2:len,2:len,2:len,2:len] = covop_est[[1]][-len,-len,-len,-len]
      cov_mat[2:len,2:len,2:len,2:len] = covop_est[[2]][-len,-len,-len,-len]

      sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_est[[1]]^2+cov_mat_A^2)*gmat)
      sigma2_cov=sigma2_cov+2*sum((1/2)*(covop_est[[2]]^2+cov_mat^2)*gmat)
    }
    sigma2_cov=2*sigma2_cov

    beta = sigma2_cov/(2 * mu_cov)
    nu = (2 * mu_cov^2)/sigma2_cov

    return(1 - pchisq(stat/beta, df = nu))
  }

  if( is.null(K) == TRUE) {
    K=20}
  stat_pvalue=gof_pvalue(error_fit,sigma_fit,basis,y_vec,sigma_2_proj_coefs,para,K)

  if( is.null(pplot) == TRUE) {
    stat_pvalue=stat_pvalue
  }
  else if(pplot==1){
    kseq=1:20
    pmat=c(rep(1,20))
    for (i in 1:20){
      pmat[i]=gof_pvalue(error_fit,sigma_fit,basis,y_vec,sigma_2_proj_coefs,para,kseq[i])
    }
    x<-1:20
    plot(x,pmat, col="black",type="l",ylim=c(0,1.2),xlab="Lag",ylab="P-Values",main="P-values of Goodness of Fit test")
    abline(h=0.1, col="red")
  }

  return(stat_pvalue)
}
