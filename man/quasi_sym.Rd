\name{quasi_sym}

\alias{quasi_sym}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Recursive computation of the conditional likelihood for the
  Quadratic Exponential Model proposed in Bartolucci & Nigro (2010)}

\description{Recursively compute the denominator of the individual
  conditional likelihood function for the Quadratic Exponential Model,
  adapted from Krailo & Pike (1984).}

\usage{
quasi_sym(eta,s,dyn=FALSE,y0=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{eta}{individual vector of products between covariate and parameters}
  \item{s}{total score of the individual}
  \item{dyn}{TRUE if in the dynamic version; FALSE for the static
    version (by default)}
  \item{y0}{Individual initial observation for dynamic models}
  
}


\value{
  \item{f}{value of the denominator}
  \item{d1}{first derivative of the recursive function}
  \item{dl1}{a component of the score function}
  \item{D2}{second derivative of the recursive function}
  \item{Dl2}{a component of the Hessian matrix}	
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F., Valentini. F., & Pigini, C. (2021), Recursive Computation of the Conditional Probability Function of the Quadratic Exponential Model for Binary Panel Data, \emph{Computational Economics}, https://doi.org/10.1007/s10614-021-10218-2.

Krailo, M. D., & Pike, M. C. (1984). Algorithm AS 196: conditional multivariate logistic analysis of stratified case-control studies, \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, \bold{33(1)}, 95-103.
}

\author{
  Francesco Bartolucci (University of Perugia), Claudia Pigini (University of Ancona "Politecnica delle Marche"), Francesco Valentini (University
  of Ancona "Politecnica delle Marche")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}