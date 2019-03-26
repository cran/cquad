\name{sim_panel_logit}

\alias{sim_panel_logit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Simulate data from the dynamic logit model}

\description{Simulate data from the dynamic logit model given a set of covariates and a vector of parameters.}

\usage{sim_panel_logit(id, al, X = NULL, eta, dyn = FALSE)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{id}{list of the reference unit of each observation}
  \item{al}{list of individual specific effects}
  \item{X}{corresponding matrix of covariates (optional)}
  \item{eta}{vector of parameters}
  \item{dyn}{TRUE if in the dynamic version; FALSE for the static version (by default)}
}

\value{
	\item{yv}{simulated vector of binary response variables}
	\item{pv}{vector of probabilities of "success"}
}

\author{
Francesco Bartolucci (University of Perugia), Claudia Pigini (University
of Ancona "Politecnica delle Marche")
}

\examples{
# simulate data from the static logit model
n = 1000; TT = 5                 # sample size, number of time occasions
id = (1:n)\%x\%rep(1,TT)           # vector of indices
al = rnorm(n)                    # simulate alpha
X = matrix(rnorm(2*n*TT),n*TT,2) # simulate two covariates
eta1 = c(1,-1)                 # vector of parameters
out = sim_panel_logit(id,al,X,eta1)
y1 = out$yv

# simulate data from the dynamic logit model
eta2 = c(1,-1,2)            # vector of parameters including state dependence
out = sim_panel_logit(id,al,X,eta2,dyn=TRUE)
y2 = out$yv
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}
