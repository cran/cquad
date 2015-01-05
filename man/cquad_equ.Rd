\name{cquad_equ}

\alias{cquad_equ}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Conditional maximum likelihood estimation for the modified version of the quadratic exponential model (to test for state dependence)}

\description{Fit by conditional maximum likelihood a modified version of the model for binary logitudinal data proposed by Bartolucci & Nigro (2010), in which the interaction terms have an extended form. This modified version is used to test for state dependence as described in Bartolucci et al. (2013).}

\usage{cquad_equ(id, yv, X = NULL, be = NULL, w = rep(1, n))}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{id}{list of the reference unit of each observation}
  \item{yv}{corresponding vector of response variables}
  \item{X}{corresponding matrix of covariates (optional)}
  \item{be}{intial vector of parameters (optional)}
  \item{w}{vector of weights (optional)}
}

\value{
	\item{lk}{conditional log-likelihood value}
	\item{be}{estimate of the regression parameters (including for the lag-response)}
	\item{scv}{matrix of individual scores}
	\item{J}{Hessian of the log-likelihood function}
	\item{se}{standard errors}
	\item{Tv}{number of time occasions for each unit}
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F., Nigro, V. and Pigini, C. (2013), Testing for state dependence in binary panel data with individual covariates, MPRA Paper 48233, University Library of Munich, Germany.
}
\author{
Francesco Bartolucci (University of Perugia), Claudia Pigini (University of Perugia)
}

\examples{
# example about PSID data
data(PSID)
PSID = PSID[1:2000,]   # only to speed up the example, to be removed in general
# static model
id = PSID[,1]; yv = PSID[,12]; X = PSID[,6:11]
out = cquad_equ(id,yv,X)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}
