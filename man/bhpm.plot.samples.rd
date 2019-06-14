\name{bhpm.plot.samples}
\alias{bhpm.plot.samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Posterior Distribution}
\description{
This function plots a graph of the sampled posterior distribution.
}
\details{
Two graphs are displayed on the same panel. The left graph is the traceplot of the chains. The right graph is
a plot of the distribution.
}
\usage{
	bhpm.plot.samples(samples, title)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samples}{
An array of samples indexed by \emph{chain}.
}
  \item{title}{
The graph title.
}
}
\value{
Nothing is returned.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
data(bhpm.cluster.data1)
raw = bhpm.cluster.BB.hier3(bhpm.cluster.data1)
sample = raw$theta[,1,1,1,1,]
bhpm.plot.samples(sample, "Sample plot")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bhpm.plot.samples}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Cluster} % __ONLY ONE__ keyword per line
\keyword{Adverse event} % __ONLY ONE__ keyword per line
\keyword{Adverse outcome} % __ONLY ONE__ keyword per line
