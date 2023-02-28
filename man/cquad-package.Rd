\name{cquad-package}

\alias{cquad-package}

\docType{package}

\title{Conditional Maximum Likelihood for Quadratic Exponential Models for Binary Panel Data}

\description{
Estimation, based on conditional maximum likelihood, of the quadratic exponential model proposed by Bartolucci & Nigro (2010) and of a simplified and a modified version of this model. The quadratic exponential model is suitable for the analysis of binary longitudinal data when state dependence (further to the effect of the covariates and a time-fixed individual intercept) has to be taken into account. Therefore, this is an alternative to the dynamic logit model having the advantage of easily allowing conditional inference in order to eliminate the individual intercepts and then getting consistent estimates of the parameters of main interest (for the covariates and the lagged response). The simplified version of this model does not distinguish, as the original model does, between the last time occasion and the previous occasions. The modified version formulates in a different way the interaction terms and it may be used to test in a easy way state dependence as shown in Bartolucci, Nigro & Pigini (2018). The package also includes estimation of the dynamic logit model by a pseudo conditional estimator based on the quadratic exponential model, as proposed by Bartolucci & Nigro (2012).}

\details{
\tabular{ll}{
Package: \tab cquad\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2019-03-26\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Francesco Bartolucci (University of Perugia, IT), Claudia Pigini
(University of Perugia, IT), Francesco Valentini (University of Ancona
"Politecnica delle Marche")

Maintainer: Francesco Bartolucci <francesco.bartolucci@unipg.it>
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F. and Nigro, V. (2012). Pseudo conditional maximum likelihood estimation of the dynamic logit model for binary panel data, \emph{Journal of Econometrics}, \bold{170}, 102-116.

Bartolucci, F. and Pigini, C. (2017). cquad: An R and Stata package for conditional maximum likelihood estimation of dynamic binary panel data models, \emph{Journal of Statistical Software}, \bold{78}, 1-26, doi:10.18637/jss.v078.i07.

Bartolucci, F., Nigro, V., & Pigini, C. (2018). Testing for state dependence in binary panel data with individual covariates by a modified quadratic exponential model. \emph{Econometric Reviews}, \bold{37(1)}, 61-88.

Bartolucci, F., Valentini. F., & Pigini, C. (2021), Recursive Computation of the Conditional Probability Function of the Quadratic Exponential Model for Binary Panel Data, \emph{Computational Economics}, https://doi.org/10.1007/s10614-021-10218-2.

Cox, D. R. (1972), The Analysis of multivariate binary data, \emph{Applied Statistics}, \bold{21}, 113-120.
}

\keyword{panel data}

\examples{
# example based on simulated data
data(data_sim)
data_sim = data_sim[1:500,]   # to speed up the example, remove otherwise
# static model
out1 = cquad(y~X1+X2,data_sim)
# dynamic model
out2 = cquad(y~X1+X2,data_sim,dyn=TRUE)
}
