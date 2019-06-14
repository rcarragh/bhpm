\name{bhpm.ptheta}
\alias{bhpm.ptheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reports the posterior probability that theta (the increase in the log-odds) is greater than zero, zero, and less than zero for each outcome}
\description{
This function reports the posterior probability that theta is positive negative or zero, i.e. that there is an increase, decrease, or no difference in the log
odds of an outcome being associated with treatment.
}
\usage{
	bhpm.ptheta(raw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{
The output from a call to one of bhpm.cluster.BB.hier3, bhpm.cluster.1a.hier3, bhpm.cluster.BB.hier2, bhpm.cluster.1a.hier2.
}
}
\value{
A data frame containing the columns: \emph{Trt.Grp}, \emph{Cluster}, \emph{Outcome.Grp}, \emph{Outcome}, \emph{ptheta}, \emph{ptheta.pos}, \emph{ptheta.zero}, \emph{ptheta.neg}.
}
\author{
R. Carragher
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

\dontrun{
data(bhpm.cluster.data1)
raw = bhpm.BB(bhpm.cluster.data1)
}

\dontrun{
p = bhpm.BB.ptheta(raw)

head(p, 2)
  Trt.Grp   Cluster Outcome.Grp   Outcome     ptheta ptheta.pos ptheta.zero ptheta.neg
1       2 0.0-180.0   Bdy-sys_1  Adv-Ev_1 0.83874074 0.83874074           0  0.1612593
2       2 0.0-180.0   Bdy-sys_1 Adv-Ev_10 0.06785185 0.06785185           0  0.9321481
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{bhpm.ptheta}
\keyword{Bayesian} % __ONLY ONE__ keyword per line
\keyword{Hierarchy} % __ONLY ONE__ keyword per line
\keyword{Cluster} % __ONLY ONE__ keyword per line
\keyword{Adverse event} % __ONLY ONE__ keyword per line
\keyword{Adverse outcome} % __ONLY ONE__ keyword per line
