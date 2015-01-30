\name{getPotentialSwaps}
\alias{getPotentialSwaps}
\title{
Calculate swap partner for TF-/miRNA-target graph. This function is usually only called internally, but may be used to speed up repetitive calls to \code{\link{birteRun}}.
}
\description{
Calculates for TF-/miRNA-target graph all potential swap partner.
}
\usage{
getPotentialSwaps(genesets, perc.overlap.cutoff=0.8, integer.id=TRUE) 
}
\arguments{
  \item{genesets}{Each entry corresponds to a regulator (miRNA, TF, other) and contains its target genes.}  
  \item{perc.overlap.cutoff}{Percentage cutoff of minimal overlap between two miRNAs or TFs to be possible swap partner.}
  \item{integer.id}{If TRUE, the swap partner are not output as characters, but as integer indices.}  
}
\value{
A list object, where each element corresponds to one regulator.  
}
\examples{
# artificial data
data(humanNetworkSimul)
genesets = c(sapply(affinities2$TF, names), 
sapply(affinities2$miRNA, names), sapply(affinities2$other, names))
swaps = getPotentialSwaps(genesets)
}
\author{
Holger Froehlich
}

