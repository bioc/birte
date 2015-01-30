\name{proposeInteractions}
\alias{proposeInteractions}
\title{
Propose possible regulator-regulator interactions that could be worthwhile to be tested into the biRte model.
}
\description{
The purpose is to identify interesting interactions terms for our model. Assuming the design matrix to be binary, the interaction term between two variables X1 and X2 is 1, if X1=1 AND X2=1 (zero, otherwise). The function computes the relative overlap between all pairs of variables X1 and X2 via the Tanimoto-Jaccard index. Variable pairs with cutoff.lower < Tanimoto-Jaccard < cutoff.upper are considered interesting candidates.
}
\usage{
proposeInteractions(affinities, cutoff.lower=0.1, cutoff.upper=0.8)
}
\arguments{
\item{affinities}{original regulator-target gene network}
\item{cutoff.lower}{lower range for Tanimoto-Jaccard index}
\item{cutoff.upper}{upper range for Tanimoto-Jaccard index}
}
\value{
A list of the same format as affinities$TF and affinites$miRNA: Each list entry corresponds to the intersection of two target gene sets.
}
\examples{
# artificial data
data(humanNetworkSimul)
affinities.int = proposeInteractions(affinities2)
}

\author{
Holger Froehlich
}

