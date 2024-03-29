% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{est.fGarch}
\alias{est.fGarch}
\title{Estimate Functional GARCH Model}
\usage{
est.fGarch(fdata, basis, p = 1, q = 1)
}
\arguments{
\item{fdata}{The functional data object with N paths.}

\item{basis}{The M-dimensional basis functions.}

\item{p}{order of the depedence on past volatilities.}

\item{q}{order of the depedence on past squared observations.}
}
\value{
List of model paramters:

d: d Parameter vector, for intercept function \eqn{\delta}.

As: A Matrices, for \eqn{\alpha} operators.

Bs: B Matrices, for \eqn{\beta} operators.
}
\description{
est.fGarch function estimates the Functional GARCH(p,q) model by using the Quasi-Maximum Likelihood Estimation method.
}
\examples{
\dontrun{
# generate discrete evaluations of the FGARCH process and smooth them into a functional data object.
yd = dgp.fgarch(grid_point=50, N=200, "garch")
yd = yd$garch_mat
fd = fda::Data2fd(argvals=seq(0,1,len=50),y=yd,fda::create.bspline.basis(nbasis=32))

# extract data-driven basis functions through the truncated FPCA method.
basis_est = basis.est(yd, M=2, "tfpca")$basis

# estimate an FGARCH(1,1) model with basis when M=1.
garch11_est = est.fGarch(fd, basis_est[,1])
}

}
\references{
Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019). Functional GARCH models: The quasi-likelihood approach and its applications. Journal of Econometrics. 209(2), 353-375. <doi:10.1016/j.jeconom.2019.01.006>.\cr
}
\seealso{
\code{\link{est.fArch}} \code{\link{est.fGarchx}} \code{\link{diagnostic.fGarch}}
}
