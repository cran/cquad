\name{sq}

\alias{sq}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Generate binary sequences}

\description{Generate binary sequences of a certain length and with a certain sum.}

\usage{sq(J, s = NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{length of the binary sequences}
  \item{s}{sum of the binary sequences (optional)}
}

\value{
	\item{M}{Matrix of binary configurations}
}

\author{
Francesco Bartolucci (University of Perugia)
}

\examples{
# generage all sequence of 5 binary variables
sq(5)
# generage all sequence of 5 binary variables, with sum equal 2
sq(5,2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{conditional inference}