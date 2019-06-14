\name{bhpm.print.summary.stats}
\alias{bhpm.print.summary.stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print the Summary Statistics of Posterior Distributions}
\description{
The function prints the variable names, the mean, the HPI interval, the standard distribtion and the
MCMC standard error for the simulated sample.}
\usage{
	bhpm.print.summary.stats(summ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summ}{
The output from a call to bhpm.summary.stats
}
}
\value{
Nothing
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
summ = bhpm.summary.stats(raw)
bhpm.print.summary.stats(summ)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bhpm.print.summary.stats}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Cluster} % __ONLY ONE__ keyword per line
\keyword{Adverse event} % __ONLY ONE__ keyword per line
\keyword{Adverse outcome} % __ONLY ONE__ keyword per line
