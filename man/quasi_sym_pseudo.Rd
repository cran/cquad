\name{quasi_sym_pseudo}

\alias{quasi_sym_pseudo}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Recursive computation of pseudo conditional maximum likelihood
method proposed by Bartolucci & Nigro (2012).}

\description{Recursively compute the denominator of the individual
  conditional likelihood function for the pseudo conditional maximum
  likelihood method proposed by Bartolucci & Nigro (2012) recursively,
  adapted from Krailo & Pike (1984).}

\usage{
quasi_sym_pseudo(eta,qi,s,y0=NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{eta}{individual vector of products between covariate and parameters}
  \item{s}{total score of the individual}
  \item{qi}{Vector of quantities from first step estimation}
  \item{y0}{Individual initial observation for dynamic models}
  
}


\value{
  \item{f}{value of the denominator}
  \item{d1}{first derivative of the recursive function}
  \item{dl1}{a component of the score function}
  \item{D2}{second derivative of the recursive function}
  \item{Dl2}{a component for the Hessian matrix}	
}

\references{
Bartolucci, F. and Nigro, V. (2010), A dynamic model for binary panel data with unobserved heterogeneity admitting a root-n consistent conditional estimator, \emph{Econometrica}, \bold{78}, 719-733.

Bartolucci, F. and Nigro, V. (2012), Pseudo conditional maximum
likelihood estimation of the dynamic logit model for binary panel data,
\emph{Journal of Econometrics}, \bold{170}, 102-116.

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