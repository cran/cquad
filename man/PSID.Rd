\name{PSID}
\alias{PSID}
\docType{data}
\title{PSID dataset}

\description{It contains data about a sample of 1446 women followrd from 1987 to 1993 according to certain variables (source PSID:  http://psidonline.isr.umich.edu).}

\usage{data(PSID)}
\format{
  A data frame with 1510 observations on the following 12 items:
  \describe{
    \item{\code{id}}{list of the reference unit of each observation}
    \item{\code{year}}{yead of each observation}
    \item{\code{race}}{race equal to 0 for which and 1 for black}
    \item{\code{age86}}{age in 1986}
    \item{\code{edu86}}{number of years of schooling in 1986}
    \item{\code{child_1_2}}{number of children of age 1-2 (in previous year)}
    \item{\code{child_3_5}}{number of children of age 3-5 (in previous year)}
    \item{\code{child_6_13}}{number of children of age 6-13 (in previous year)}
    \item{\code{child_14}}{number of children of age 14- (in previous year)}
    \item{\code{income}}{household income (previous year)}
    \item{\code{fertility}}{binary variable equal to 1 if the woman has given birth to a child}
    \item{\code{occupation}}{binary variable equal to 1 if the woman is employed}
  }
}

\references{
Bartolucci, F. and Farcomeni, A. (2009), A multivariate extension of the dynamic logit model for longitudinal data based on a latent Markov heterogeneity structure, \emph{Journal of the American Statistical Association},  \bold{104}, pp. 816-831.	
	
Bartolucci, F., Nigro, V. and Pigini, C. (2013), Testing for state dependence in binary panel data with individual covariates, MPRA Paper 48233, University Library of Munich, Germany.
}

\examples{
data(PSID)
head(PSID)
}

\keyword{datasets}
\keyword{panel data}