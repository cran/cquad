\name{data_sim}
\alias{data_sim}
\docType{data}
\title{Simulated dataset}

\description{It contains a dataset simulated from the dynamic logit model}

\usage{data(data_sim)}
\format{
  The observations are for 1000 sample units at 5 five time occasions:
  \describe{
    \item{\code{id}}{list of the reference unit of each observation}
    \item{\code{time}}{number of the time occasion}
    \item{\code{X1}}{first covariate}
    \item{\code{X2}}{second covariate}
    \item{\code{y}}{response}
  }
}

% \references{
% Bartolucci, F. and Farcomeni, A. (2009), A multivariate extension of the dynamic logit model for longitudinal data based on a latent Markov heterogeneity structure, \emph{Journal of the American Statistical Association},  \bold{104}, pp. 816-831.	
	
% Bartolucci, F., Nigro, V. and Pigini, C. (2013), Testing for state dependence in binary panel data with individual covariates, MPRA Paper 48233, University Library of Munich, Germany.
% }

\examples{
data(data_sim)
head(data_sim)
}

\keyword{datasets}
\keyword{panel data}