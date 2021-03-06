% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgp.R
\name{dgp.fiid}
\alias{dgp.fiid}
\title{Data Generating Process - Independent Process}
\usage{
dgp.fiid(grid_point, N)
}
\arguments{
\item{grid_point}{The number of grid point in each curve observation.}

\item{N}{The sample size.}
}
\value{
A (grid_point) x (number of observations) matrix for iid sequences, where the finite realization of curves are stored in columns.
}
\description{
dgp.fiid function generates iid functional curve data following the Ornstein–Uhlenbeck process.
}
\details{
\eqn{x_i(t)=e^{-t/2}W_i(e^t)}, \eqn{t \in [0,1]},\cr
where \eqn{W_i(t)} is a standard Brownian Motion.
}
\examples{
# generate discrete evaluations of 100 iid curves that each curve is realized on 50 grid points.
yd_iid = dgp.fiid(50, 100)

# smooth discrete data into functional curves.
fd = fda::Data2fd(argvals=seq(0,1,len = 50),y=yd_iid,fda::create.bspline.basis(nbasis = 32))

}
\seealso{
\code{\link{dgp.fgarch}}
}
