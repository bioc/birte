\name{simplify}
\alias{simplify}
\title{
Simplify regulator-target gene network via clustering.
}
\description{
Several regulators (specifically TFs) might have extremely overlapping target genes. In order to identify clusters of highly similar regulators (mainly TFs) we implemented a network simplification algorithm in biRte: We construct the biadjacency matrix of the complete bipartite regulator target-gene graph and then calculate a single linkage clustering of regulators based on the Tanimoto-Jaccard similarity of their target genes. The dendrogram is cut at a defined height (default: 0.1) to idenfity resulting groups. The algorithm is meant to simplify the inference of active regulators, because the resulting regulator clusters have more dissimilar target gene profiles.
}
\usage{
simplify(affinities, cutoff=0.9)
}
\arguments{
\item{affinities}{original regulator-target gene network}
\item{cutoff}{cut dendrogram at height 1 - cutoff (i.e. similarity cutoff)}
}
\value{
clustered / simplified network
}
\examples{
# artificial data
data(humanNetworkSimul)
affinities2 = simplify(affinities2)
}

\author{
Holger Froehlich
}

