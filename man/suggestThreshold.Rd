\name{suggestThreshold}
\alias{suggestThreshold}
\title{
Automatically suggest suitable threshold for marginal regulator activities.
}
\description{
The algorithm fits a mixture of a beta(1,beta) and beta(alpha,1) distribution to observed marginal regulator activities. Based on this mixture a cutoff is chosen such that the expected false positive rate is below a defined threshold.
}
\usage{
suggestThreshold(prob, fpr=0.001)
}
\arguments{
\item{prob}{marginal probability obtained from \code{\link{birteRun}}}
\item{fpr}{threshold for accepted false positive rate}
}
\value{
a cutoff for marginal activity probabilities
}
\references{
Froehlich, H. and Klau, G. (2013). Reconstructing Consensus Bayesian Network Structures with Application to Learning Molecular Interaction Networks. In: Beissbarth, T., Kollmar, M., Leha, A., Morgenstern, B., Schultz, A.-K., Waack, S., and Wingender, E., editors, Proc. German Conference on Bioinformatics, Open Access Series in Informatics, pages 46 - 55. Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik, Dagstuhl Publishing, Germany.
}
\examples{
freq = 0.2*rbeta(100, 1, 10) + 0.8*rbeta(100, 5, 1)
thresh = suggestThreshold(freq)
}
\author{
Holger Froehlich
}
