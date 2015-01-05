\name{cquad-package}

\alias{cquad-package}

\alias{cquad}

\docType{package}

\title{Conditional ML for Quadratic Exponential Models for Binary Panel Data}

\description{
The package estimates, by conditional maximum likelihood, the quadratic exponential model proposed by Bartolucci & Nigro (2010) and a simplified and a modified version of this model. The quadratic exponential model is suitable for the analysis of binary longitudinal data when state dependence (further to the effect of the covariates and a time-fixed individual intercept) has to be taken into account. Therefore, this is an alternative to the dynamic logit model having the advantage of easily allowing conditional inference in order to eliminate the individual intercepts and then getting consistent estimates of the parameters of main interest (for the covariates and the lagged response). The simplified version of this model does not distinguish, as the original model does, between the last time occasion and the previous occasions. The modified version formulates in a different way the interaction terms and it may be used to test in a easy way state dependence as shown in Bartolucci, Nigro & Pigini (2013). The package also includes estimation of the dynamic logit model by a pseudo conditional estimator based on the quadratic exponential model, as proposed by Bartolucci & Nigro (2012).}

\details{
\tabular{ll}{
Package: \tab cquad\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-01-05\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Francesco Bartolucci (University of Perugia, IT), Claudia Pigini (University of Perugia, IT)

Maintainer: Francesco Bartolucci <francesco.bartolucci@unipg.it>
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F. and Nigro, V. (2012). Pseudo conditional maximum likelihood estimation of the dynamic logit model for binary panel data, \emph{Journal of Econometrics}, \bold{170}, 102-116.

Bartolucci, F., Nigro, V. and Pigini, C. (2013). Testing for state dependence in binary panel data with individual covariates, MPRA Paper 48233, University Library of Munich, Germany.

Cox, D. R. (1972), The Analysis of multivariate binary data, \emph{Applied Statistics}, \bold{21}, 113-120.
}

\keyword{panel data}

\examples{
# example about PSID data
data(PSID)
PSID = PSID[1:2000,]   # only to speed up the example, to be removed in general
# static model
id = PSID[,1]; yv = PSID[,12]; X = PSID[,6:11]
out1 = cquad_basic(id,yv,X)
# dynamic model
out2 = cquad_basic(id,yv,X,dyn=TRUE)
}
