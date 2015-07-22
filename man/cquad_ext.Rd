\name{cquad_ext}

\alias{cquad_ext}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Conditional maximum likelihood estimation of the quadratic exponential model for panel data}

\description{Fit by conditional maximum likelihood the model for binary logitudinal data proposed by Bartolucci & Nigro (2010).}

\usage{cquad_ext(id, yv, X = NULL, be = NULL, w = rep(1, n))}
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
	\item{ser}{robust standard errors}	
	\item{Tv}{number of time occasions for each unit}
}
	
\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator. \emph{Econometrica}, \bold{78}, pp. 719-733.
}

\author{
Francesco Bartolucci (University of Perugia), Claudia Pigini (University of Ancona "Politecnica delle Marche")
}

\examples{
# example based on simulated data
data(data_sim)
data_sim = data_sim[1:500,]   # to speed up the example, remove otherwise
id = data_sim$id; yv = data_sim$y; X = cbind(X1=data_sim$X1,X2=data_sim$X2)
# static model
out = cquad_ext(id,yv,X)
summary(out)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}