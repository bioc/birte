\name{birteRun}
\alias{birteRun}
\alias{birteLimma}
\title{
Main interface for Bayesian Inference of Regulatory Influence on Expression (biRte).
}
\description{
The function birteRun estimates regulator activities from gene expression data plus a given regulator-target gene network via MCMC sampling. The function assumes experimental data to be of one of the following formats:
i) expression matrix (#genes x #samples) 
ii) mRNA log fold changes

In the first case it is now also possible to estimate regulator activity for each individual sample. In addition one can estimate condition specific regulator activity, if there are exactly two conditions.

In addition to miRNA, mRNA and TF expression data, biRte also allows for integrating other data types (e.g. CNV data - data type 'other'). Note that in such a case also the corresponding regulator-target gene relationships have to be defined. 

birteLimma is a convenience function, which allows to directly pass results from a previous limma analysis (see \code{\link{limmaAnalysis}}) to birteRun. When working with relative expression levels (log fold changes) birteLimma allows to deal with arbitrary complex statistical designs, including e.g. time as a covariate.
}
\usage{
birteRun(dat.mRNA, mRNA.Sigma=NULL, nrep.mRNA=c(5, 5), 
    df.mRNA=sum(nrep.mRNA)-2,
    data.regulators=NULL, sigma.regulators=NULL, 
    nrep.regulators=NULL, diff.regulators=NULL, 
  	init.regulators=NULL, theta.regulators=NULL, 
    reg.interactions=FALSE, affinities, use.affinities=FALSE,    
    niter=100000, nburnin=100000, thin=50, 
    potential_swaps=NULL, only_switches=FALSE, only.diff.TFs=TRUE, 
    explain.LFC=TRUE, 
    single.sample=FALSE, single.sample.estimator=c("mpost", "MAP"),
    model=c("no-plug-in", "all-plug-in"))
    
birteLimma(dat.mRNA=NULL, limmamRNA,
  	data.regulators=NULL, limma.regulators=NULL, 
    fdr.regulators=NULL, lfc.regulators=NULL,
		init.regulators=NULL, theta.regulators=NULL, 
    reg.interactions=FALSE, affinities, use.affinities=FALSE,
	  niter=100000, nburnin=100000, thin=50, 
    potential_swaps=NULL, only_switches=FALSE, only.diff.TFs=TRUE, 
    explain.LFC=TRUE, 
    single.sample=FALSE, single.sample.estimator=c("mpost", "MAP"),
    model=c("no-plug-in", "all-plug-in"))
}
\arguments{
  \item{dat.mRNA}{mRNA expression data matrix. IMPORTANT: Replicates must be ordered according to nrep.mRNA}
  \item{mRNA.Sigma}{gene expression variances (array data). IMPORTANT: Names have to match the row names in dat.mRNA. If mRNA.Sigma = NULL, birte.run tries to deduce variances from a limma analysis (see \code{\link{limmaAnalysis}})}
  \item{nrep.mRNA}{number of replicates per condition.}
   \item{df.mRNA}{residual degrees of freedom of linear model for mRNA data.}
  
  \item{data.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains one data matrix. IMPORTANT: Samples in data matrices have to be grouped according to conditions. That means first there are all samples from the first condition, then those from the second condition, etc.}
  \item{sigma.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains one named vector of expression variances (array data) or dispersion parameters (RNAseq data). IMPORTANT: Names have to fit to the row names of the data matrix.}  
  \item{nrep.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains the number of replicates per condition} 
  \item{diff.regulators}{list with at most 3 components (miRNA, TF, other). Each component is a character vector with differentially expressed regulators. Has to be subset of row names of the data matrix}
  \item{init.regulators}{list with at most 3 components (miRNA, TF, other). Each component is matrix of #conditions x length(affinities[[regulator type]]): initial states for regulators. In case this matrix is not provided (i.e. NULL) initial states are assumed to be 0. IMPORTANT: column names have to match names(affinities[[regulator type]])}
  \item{theta.regulators}{list with at most 3 components (miRNA, TF, other). If single numbers are provided, each component contains the expected fraction of active regulators. If vectors are provided, each vector entry corresponds to the individual probability of a specific regulator to be active. Accordingly, vectors should be named in agreement with the regulator-target gene network. If affinities$other corresponds to interaction terms between regulators, theta.regulators can also be provided as a #regulators x #regulators matrix.}
  \item{reg.interactions}{If TRUE, entries of affinities$other are interpreted as interaction terms between regulators.}
  
  \item{affinities}{Regulator-target gene interactions. This is a list with at most three components (TF, miRNA, other). Each of these lists again contains a weighted adjacency list representation. See \code{\link{humanNetworkSimul}} for an example. IMPORTANT: gene names used in this network have to match with row names of dat.mRNA. Moreover, regulator names have to fit to row names of the corresponding data matrices. It is assumed that in presence of regulators of type "other", there exist also miRNA regulators. If you want to circumvent this behavior (e.g. you want to use CNVs, but *no* miRNAs), store these regulators in the slot "miRNA".} 
  \item{use.affinities}{Should weights given in the bipartite regulator-target gene graph given a specific meaning? If yes, it is assumed that weights correspond to quantitative influences of regulators on their targets.}
  
  \item{niter}{Number of MCMC iterations (AFTER burnin).}
  \item{nburnin}{Number of MCMC iterations UNTIL burnin is assumed to be finished.}
  \item{thin}{Thinning of Markov chain: only use every thin's sample for posterior computation.}
  \item{potential_swaps}{Pre-computed potential swaps (OPTIONAL, see get_potential_swaps).}
  \item{only_switches}{Should only switches be performed?}    
  \item{only.diff.TFs}{Should, in case of TF expression data, only the information for differentially expressed TFs be considered? Note that this makes fewer assumption about the relation of mRNA and protein expression data, but typically leads to less conservative results (i.e. more TFs predicted to be active).}
  \item{limmamRNA}{results of limma analysis for mRNA data according to \code{\link{limmaAnalysis}}}
  \item{limma.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains the results of a limma analysis for regulator data according to \code{\link{limmaAnalysis}}}
  \item{lfc.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains the log fold change cutoff for differential expression. It is assumed to be 0, if not provided.}
  \item{fdr.regulators}{list with at most 3 components (miRNA, TF, other). Each component contains the FDR cutoff for differential expression (DEFAULT: 0.05).}
  \item{explain.LFC}{If yes, biRte tries to explain mRNA log fold changes rather than expression levels itself.}
  \item{single.sample}{If yes, biRte tries to explain mRNA data for each individual sample. The output is a #samples x #regulators matrix.}
  \item{single.sample.estimator}{Which type of estimate for regulator activity should be provided? mpost = marginal posterior activity probability for each regulator; MAP = most likely regulator configuration found during MCMC sampling}
  \item{model}{If "no-plug-in", for marginal log likelihoods are considered for regulator specific expression data. Otherwise, (posterior) variance estimates are used directly.}
}
\value{
If single.sample is FALSE, the function returns a list containing the following entries:
  \item{post}{#regulators x #conditions matrix containing the marginal probability for each regulator to influence mRNA expression.}
  \item{map}{#regulators x #conditions matrix containing the regulator configuration with highest joint probability.}       
  \item{coef}{matrix of #conditions x #replicates: Each entry is itself a matrix (embedded into a list) of #coefficients x #effective samples. The data contains the expected regression coefficients.}
  \item{log_lik_trace}{(Marginal) log-likelihood trace of MCMC sampling.}
  \item{eff_sample_size}{effective sample size after burnin and thinning}  
  
  \item{contains.interactions}{TRUE, if affinities$other corresponds to interaction terms between regulators, FALSE otherwise}
  \item{explain.LFC}{TRUE, if model explains mRNA log fold change, FALSE otherwise}
  \item{nburnin}{number of burnin iterations - as provided as an argument}
  \item{affinities}{original regulator-target gene network}
  \item{C_cnt}{number of conditions}  
  \item{design}{design matrix effectively used for model training}
  \item{param}{estimated parameters for mRNA precision (i.e. inverse variance) distribution}
  
If single.sample is TRUE the output is a #samples x #regulator matrix.
}
\author{
Holger Froehlich
}
\examples{
# artificial data
data(humanNetworkSimul)
sim = simulateData(affinities2)
limmamRNA = limmaAnalysis(sim$dat.mRNA, design=NULL, "treated - control")
limmamiRNA = limmaAnalysis(sim$dat.miRNA, design=NULL, "treated - control")
limmaTF = limmaAnalysis(sim$dat.TF, design=NULL, "treated - control")

# burnin and sampling size is much too small in reality
result = birteLimma(dat.mRNA=sim$dat.mRNA, 
data.regulators=list(miRNA=sim$dat.miRNA, TF=sim$dat.TF), 
limmamRNA=limmamRNA, limma.regulators=list(miRNA=limmamiRNA, TF=limmaTF), 
affinities=affinities2, niter=100, nburnin=100, thin=2) 
plotConvergence(result)
pred = birtePredict(result, rownames(sim$dat.mRNA))
MSE.Bayes = mean((pred[[1]][[1]]$mean - limmamRNA$pvalue.tab[rownames(sim$dat.mRNA),"logFC"])^2)
MSE.Bayes


# real data
library(Biobase)
data(EColiOxygen)
# prepare network
affinities = list(TF=sapply(names(EColiNetwork$TF), function(tf){
w = rep(1, length(EColiNetwork$TF[[tf]])); 
names(w)= EColiNetwork$TF[[tf]]; w}))
# prepare data
mydat = exprs(EColiOxygen)
colnames(mydat) = make.names(paste(pData(EColiOxygen)$GenotypeVariation, 
pData(EColiOxygen)$GrowthProtocol, sep="."))
mydat = cbind(mydat[,colnames(mydat) =="wild.type.aerobic"], 
exprs(EColiOxygen)[,colnames(mydat) == "wild.type.anaerobic"])

# more realistic sampling
\dontrun{
result = birteRun(dat.mRNA=mydat, 
nrep.mRNA=c(3,4), affinities=affinities, niter=10000, nburnin=10000)
plotConvergence(result)
}
}

