\name{estimateNetwork}
\alias{estimateNetwork}
\title{
Estimate network between active regulators using Nested Effects Models (NEMs).
}
\description{
Given a biRte model, this function makes posterior inference about possible upstream-downstream relationships between active regulators. This is done based on observed differential expression of putative target genes. The idea is that regulator A acts upstream of regulator B, if differentially expressed targets of B are a subset of those of A.
}
\usage{
estimateNetwork(model, thresh=0.1, select=c("marginal", "MAP"), method="pairwise",
de.genes, bootstrap=0, typeII=0.1)
}
\arguments{
\item{model}{biRte model}
\item{thresh}{cutoff for marginal posterior probabilities}
\item{select}{"marginal": select regulators based on marginal posterior probabilities; "MAP": select regulators based on MAP configuration}
\item{method}{algorithm used for NEM based network inference, see \code{\link{nem}}}
\item{de.genes}{set of differentially expressed genes}

\item{bootstrap}{optional: number of bootstrap replicates to draw (non-parameteric bootstrap)}
\item{typeII}{assumed type-II error rate}
}
\value{
nem-model
}
\examples{
# see vignette
}
\author{
Holger Froehlich
}

