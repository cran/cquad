\name{cquad_pseudo}

\alias{cquad_pseudo}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Pseudo conditional maximum likelihood estimation of the dynamic logit model}

\description{Estimate the dynamic logit model for binary longitudinal data by the pseudo conditional maximum likelihood method proposed by Bartolucci & Nigro (2012).}

\usage{
cquad_pseudo(id, yv, X = NULL, be = NULL, w = rep(1,n), Ttol=10)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{id}{list of the reference unit of each observation}
  \item{yv}{corresponding vector of response variables}
  \item{X}{corresponding matrix of covariates (optional)}
  \item{be}{initial vector of parameters (optional)}
  \item{w}{vector of weights (optional)}
  \item{Ttol}{Threshold individual observations that activates the recursive algorithm (default=10)}
}


\value{
	\item{formula}{formula defining the model}
	\item{lk}{conditional log-likelihood value}
	\item{coefficients}{estimate of the regression parameters (including for the lag-response)}
	\item{vcov}{asymptotic variance-covariance matrix for the parameter estimates}	
	\item{scv}{matrix of individual scores}
	\item{J}{Hessian of the log-likelihood function}
	\item{se}{standard errors}
	\item{se2}{robust standard errors that also take into account the first step}
	\item{Tv}{number of time occasions for each unit}
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F. and Nigro, V. (2012), Pseudo conditional maximum likelihood estimation of the dynamic logit model for binary panel data, \emph{Journal of Econometrics}, \bold{170}, 102-116.
}

\author{
  Francesco Bartolucci (University of Perugia), Claudia Pigini (University of Ancona "Politecnica delle Marche"), Francesco Valentini (University
  of Ancona "Politecnica delle Marche")
}

\examples{
\dontrun{
# example based on simulated data
data(data_sim)
data_sim = data_sim[1:500,]   # to speed up the example, remove otherwise
id = data_sim$id; yv = data_sim$y; X = cbind(X1=data_sim$X1,X2=data_sim$X2)
# estimate dynmic logit model
out = cquad_pseudo(id,yv,X, Ttol=10)
summary(out)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}