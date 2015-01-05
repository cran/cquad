\name{cquad_pseudo}

\alias{cquad_pseudo}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Pseudo conditional maximum likelihood estimation of the dynamic logit model}

\description{Estimate the dynamic logit model for binary logitudinal data by the pseudo conditional maximum likelihood method proposed by Bartolucci & Nigro (2012).}

\usage{
cquad_pseudo(id, yv, X = NULL, be = NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{id}{list of the reference unit of each observation}
  \item{yv}{corresponding vector of response variables}
  \item{X}{corresponding matrix of covariates (optional)}
  \item{be}{intial vector of parameters (optional)}
}


\value{
	\item{lk}{conditional log-likelihood value}
	\item{be}{estimate of the regression parameters (including for the lag-response)}
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
Francesco Bartolucci (University of Perugia), Claudia Pigini (University of Perugia)
}

\examples{
\dontrun{
# example on PSID data with addition of time dummies
data(PSID)
# estimate dynmic logit model
id = PSID[,1]; yv = PSID[,12]; X = PSID[,6:11]
X[,5] = X[,5]/1000             # rescale income
TD = rep(1,1446)\%x\%rbind(0,0,diag(5))
colnames(TD) = c("y1989","y1990","y1991","y1992","y1993") 
X = cbind(X,TD)  
out = cquad_pseudo(id,yv,X)
summary(out)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}