
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CurVol

<!-- badges: start -->

<!-- badges: end -->

The goal of CurVol is to document useful functions to analyze the
volatility of functional time series data. Methods and tools in this
package replicate hypothesis testing, model estimation, and backtesting
in a series of papers:

Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the
ARCH model. Econometric Theory, 29(2), 267-288.

Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized
autoregressive conditional heteroskedasticity. Journal of Time Series
Analysis, 38(1), 3-21.

Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019).
Functional GARCH models: The quasi-likelihood approach and its
applications. Journal of Econometrics, 209(2), 353-375.

Rice, Wirjanto, and Zhao (2020) Tests for conditional heteroscedasticity
of functional data, Journal of Time Series Analysis. 41(6), 733-758.

Rice, Wirjanto, and Zhao (2020) Forecasting Value at Risk via intra-day
return curves, International Journal of Forecasting.

Rice, Wirjanto, and Zhao (2020) Functional GARCH-X models with an
application to forecasting crude oil return curves, Working paper.

## Installation

You can install the released version of CurVol from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CurVol")
```

## R topics documented

1.  backtest.var - backtest the intra-day VaR forecasts.
2.  basis.est - estimate non-negative basis functions.
3.  dgp.fgarch - generate functional data following the functional
    ARCH(1) or GARCH(1,1) process.
4.  dgp.fiid - generate iid functional data following Ornstein–Uhlenbeck
    process.
5.  diagnostic.fGarch - estimation parameters as the inputs for
    diagnostic purposes.
6.  est.fArch - estimate functional ARCH (q) model.
7.  est.fGarch - estimate Functional GARCH (p,q) model.
8.  est.fGarchx - estimate Functional GARCH-X model.
9.  fun\_hetero - test conditional heteroscedasticity for functional
    data.
10. gof.fgarch - goodness-of-fit test for functional ARCH/GARCH model.
11. intra.return - intra-day return curves: intra-day return (IDR),
    cumulative intra-day return (CIDR), overnight cumulative intra-day
    return (OCIDR).
12. sample\_data - a sample data containing S\&P 500 intra-day price.
13. var.forecast - forecast daily/intra-day Value-at-Risk.
14. var.vio - a violation process for the intra-day VaR curves.
