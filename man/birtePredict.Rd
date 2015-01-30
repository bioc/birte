\name{birtePredict}
\alias{birtePredict}
\title{
Prediction of gene expression via biRte.
}
\description{
Given a biRte model, this function makes posterior inference about gene expression data. 
}
\usage{
birtePredict(model, test.genes, method=c("Bayes", "MAP"))
}
\arguments{
\item{model}{output of \code{\link{birteRun}}}
\item{test.genes}{Set of gene IDs. Gene IDs should be contained into the defined regulator-target gene network. Note that  expression data is generally not required to be available for these genes.}
\item{method}{Bayes: estimate expectation of posterior predictive distribution. MAP: Use previously fitted ridge regression model (\code{\link{birteFitRidge}}).}
}
\value{
#conditions x #replicates matrix containing data.frame objects with
\item{gene}{gene, for which predictions are made}
\item{mean}{expected expression}
\item{sd}{SD of predictions (only for method="Bayes")}
}
\examples{
# artificial data
data(humanNetworkSimul)
sim = simulateData(affinities2)
limmamRNA = limmaAnalysis(sim$dat.mRNA, design=NULL, "treated - control")

# burnin and sampling size is much too small in reality
result = birteLimma(dat.mRNA=sim$dat.mRNA, data.regulators=NULL, 
limmamRNA=limmamRNA, 
affinities=affinities2, niter=100, nburnin=100, thin=2)

est = birtePredict(result, rownames(sim$dat.mRNA))
}
\author{
Holger Froehlich
}

