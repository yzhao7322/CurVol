
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CurVol

<!-- badges: start -->

<!-- badges: end -->

The goal of CurVol is to document useful functions to analyse the
volatility of functional curve data. Methods and tools in this package
replicate the hypothesis testing, model estimation, and backtesting
Value-at-risk in a series of papers including:

Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the
ARCH model. Econometric Theory, 29(2), 267-288.

Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized
autoregressive conditional heteroskedasticity. Journal of Time Series
Analysis, 38(1), 3-21.

Cerovecki, C., Francq, C., Hormann, S., Zakoian, J. M. (2019).
Functional GARCH models: The quasi-likelihood approach and its
applications. Journal of Econometrics, 209(2), 353-375.

Rice, Wirjanto, and Zhao (2020) Tests for conditional heteroscedasticity
of functional data, Journal of Time Series Analysis.

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

1.  backtest.inde - backtest the independence of intra-day VaR.
2.  backtest.unbias - backtest the unbiasedness of intra-day VaR.
3.  basis.fsnn - sparse and non-negative functional principal
    components.
4.  basis.pf - non-negative predictive factors.
5.  basis.pp - exponential and Bernstein basis functions.
6.  basis.tfpca - non-negative truncated functional principal
    components.
7.  basis.score - functional scores by projecting the squared process
    onto given basis functions.
8.  dgp.farch - generate functional data following the functional
    ARCH(1) process.
9.  dgp.fgarch - generate functional data following the functional
    GARCH(1,1) process.
10. dgp.fiid - generate iid functional data following Ornstein–Uhlenbeck
    process.
11. diagnostic.fGarch - estimation parameters as the inputs for
    diagnostic purposes.
12. est.fArch - estimate functional ARCH (q) model.
13. est.fGarch - estimate Functional GARCH (p,q) model.
14. est.fGarchx - estimate Functional GARCH-X model.
15. fun\_hetero - test conditional heteroscedasticity for functional
    data.
16. gof.fgarch - goodness-of-fit test for functional ARCH/GARCH model.
17. intra.return - intra-day return curves: intra-day return (IDR),
    cumulative intra-day return (CIDR), overnight cumulative intra-day
    return (OCIDR).
18. sample\_data - a sample data containing S\&P 500 intra-day price.
19. var.forecast - forecast daily/intra-day Value-at-Risk.
20. var.vio - calculate the violation process based on Value-at-Risk.
