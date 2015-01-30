\name{humanNetworkSimul}
\docType{data}
\alias{humanNetworkSimul}
\alias{affinities2}
\title{Subset of regulator-target gene network for human}
\description{
The human regulatory network was constructed as follows: For miRNA target gene prediction we used MiRmap (Vejnar and Zdobnov, 2012) and converted the reported scores into z-scores. P-values were then calculated based on the observed approximate normality of z-scores under the null hypothesis and corrected for multiple testing via the FDR method under dependencies of tests (Benjamini and Yekutieli, 2001). In consequence we arrived at 356 miRNAs regulating between 1 and 318 target genes (median: 5).

A TF-target gene network was compiled by computing TF binding affinities to promoter sequences of all human genes according to the TRAP model (Roider et al., 2007) via the author's R implementation. Upstream sequences of genes were retrieved here from the ENSEMBL database via biomaRt (Durinck et al., 2009). We assumed that promoter sequences were located in the range 0 - 2Kbp upstream to the transcription start site of a gene. 556 TRANSFAC (Wingender, 2008; public version) as well as 130 JASPAR (Bryne et al., 2008; accessed Dec. 2011) TFBS matrices were used. As significant we considered those TFBS, for which a FDR < 5\% was reported (Benjamini-Yekutieli method). 

In a subsequent step from each significant TFBS we extracted the set of TFs, for which the corresponding TFBS matrix had been defined by parsing the original TRANSFAC and JASPAR files, respectively. Ambiguities were corrected via manual curation with the help of the commercial MetaCore software. In consequence we arrived at a TF-target gene network of 344 TFs regulating between 1 and 19,392 target genes (median: 155). 

For simulation purposes we randomly selected subset of 999 human genes together with the set of corresponding regulators. After network simplification (\code{\link{simplify}}) we arrived at 13 miRNA and 111 TF clusters.
}
\usage{humanNetworkSimul}
\format{\describe{
  \item{affinities2}{TF and miRNA target gene network. For miRNAs scores from z-scores derived from MiRmap are reported. For TFs -log10(p-values).}  
 }
}
\value{see above}
\details{
  \itemize{
		\item{affinities2}{
			\itemize{
			  \item{TF}{TF-target gene graph -- list object}
			  \item{miRNA}{miRNA-target gene graph -- list object}
			 }
			Both, the TF-target and the miRNA target gene graphs are given in form of a weighted adjacency list representation.
		}		
	}	
}
\note{
	The object 'affinities2' may be extended by one further component 'other' containing additional regulatory information (e.g. CNVs, predefined interaction terms between regulators).
}
\keyword{datasets}
\references{
Vejnar, C. E. and Zdobnov, E. M. (2012). Mirmap: comprehensive prediction of microrna target repression strength. Nucleic Acids Res, 40(22):11673-11683.

Roider, H. G., Kanhere, A., Manke, T., and Vingron, M. (2007). Predicting transcription factor affinities to dna from a biophysical model. Bioinformatics, 23(2):134-141.
}
