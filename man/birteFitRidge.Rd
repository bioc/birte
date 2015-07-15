\name{birteFitRidge}
\alias{birteFitRidge}
\title{
Fit ridge regression model given a defined set of active regulators.
}
\description{
Given a most likely configuration of active regulators identified by biRte, this method fits a conventional ridge regression model to explain gene expression. This function is required, if one would like to use MAP based prediction of gene expression instead of Bayesian predictions (see \code{\link{birtePredict}}). To fit the ridge regression model the R-package \code{ridge} is employed, which provides an efficient tuning of the regularization hyperparameter. 
}
\usage{
birteFitRidge(model, mRNA.train, ref.cond=1)
}
\arguments{
\item{model}{output of \code{\link{birteRun}}}
\item{mRNA.train}{vector of gene expression values}
\item{ref.cond}{condition to consider}
}
\value{
an object of class "ridgeLinear" (see \code{\link{linearRidge}})
}
\details{
In order to make predictions with the fitted ridge regression model (\code{\link{birtePredict}}) store it into a slot "fit.ridge" of the object returned by \code{\link{birteRun}} and \code{\link{birteLimma}}, respectively.
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

fit.ridge = birteFitRidge(result, sim$dat.mRNA[,1])
}
\author{
Holger Froehlich
}

