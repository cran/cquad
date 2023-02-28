\name{cquad}

\alias{cquad}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Interface for functions fitting different versions of cquad}

\description{Fit by conditional maximum likelihood each of the models in cquad package.}

\usage{cquad(formula, data, index = NULL, model = c("basic","equal","extended","pseudo"),
             w = rep(1, n), dyn = FALSE, Ttol=10)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{formula}{formula with the same syntax as in plm package}
  \item{data}{data.frame or pdata.frame}
  \item{index}{to denote panel structure as in plm package}
  \item{model}{type of model = "basic", "equal", "extended", "pseudo"}
  \item{w}{vector of weights (optional)}
  \item{dyn}{TRUE if in the dynamic version; FALSE for the static version (by default)}
  \item{Ttol}{Threshold individual observations that activates the recursive algorithm (default=10)}
}

\value{
	\item{formula}{formula defining the model}
	\item{lk}{conditional log-likelihood value}
	\item{coefficients}{estimate of the regression parameters}
	\item{vcov}{asymptotic variance-covariance matrix for the parameter estimates}		
	\item{scv}{matrix of individual scores}
	\item{J}{Hessian of the log-likelihood function}
	\item{se}{standard errors}
	\item{ser}{robust standard errors}
	\item{Tv}{number of time occasions for each unit}
}

\author{
Francesco Bartolucci (University of Perugia), Claudia Pigini (University
of Ancona "Politecnica delle Marche"), Francesco Valentini (University
of Ancona "Politecnica delle Marche")
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# example based on simulated data
data(data_sim)
data_sim = data_sim[1:500,]   # to speed up the example, remove otherwise
# basic (static) model
out1 = cquad(y~X1+X2,data_sim)
summary(out1)
# basic (dynamic) model
out2 = cquad(y~X1+X2,data_sim,dyn=TRUE)
summary(out2)
# equal model
out3 = cquad(y~X1+X2,data_sim,model="equal")
summary(out3)
# extended model
out4 = cquad(y~X1+X2,data_sim,model="extended")
summary(out4)
# psuedo CML for dynamic model
out5 = cquad(y~X1+X2,data_sim,model="pseudo")
summary(out5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{panel data}
