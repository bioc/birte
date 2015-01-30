\name{simulateData}
\alias{simulateData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate expression data.
}
\description{
The function draws expression data from a multivariate normal distribution with block structured co-variance matrix. First, data is drawn for a control condition (no active regulators). Then data is generated for the situation that a certain fraction of regulators is turned 'on' (treatment condition). Regulator activity states are sampled from a Bernoulli distribution.
}
\usage{
simulateData(affinities, nrep = 5, miRNAExpressions = TRUE, fn.targets = 0.1, 
fp.targets = 0.2, exp.nTF = 5, exp.nmiR = 5, exp.interact = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{affinities}{
regulator-target gene network (see \code{\link{humanNetworkSimul}})
}
  \item{nrep}{
number of replicates per condition
}
  \item{miRNAExpressions}{
Should miRNA expression data be simulated?
}
  \item{fn.targets}{
fraction of false negative target predictions (i.e. missing edges per regulator in the bipartite regulator-gene graph)
}
  \item{fp.targets}{
fraction of false positive target predictions
}
 
  \item{exp.nTF}{
expected number of active TFs
}
  \item{exp.nmiR}{
expected number of active miRNAs
}

  \item{exp.interact}{
expected number of active interaction terms
}
}
\details{
If active interaction terms should be simulated, a set of possible interaction terms has to be defined in \code{affinities$other}. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{dat.mRNA }{mRNA data -- active regulators are expected to induce a log FC of 1}
\item{dat.miRNA }{miRNA data -- active miRNAs are expected to show a log FC of 1}
\item{dat.TF}{TF expression data -- active miRNAs are expected to show a log FC of 0.5}
\item{miRNAstates}{simulated miRNA activities in treatment condition}
\item{TFstates}{simulated TF activities in treatment condition}
\item{inter.states}{simulated regulator interaction activities in treatment condition}
%% ...
}
\author{
Holger Froehlich
}
\examples{
data(humanNetworkSimul)
sim = simulateData(affinities2)
}
