\name{limmaAnalysis}
\alias{limmaAnalysis}
\title{
Simple limma analysis on expression data with one contrast.
}
\description{
This is a convenience function, which can be used to perform a limma analysis with one contrast. Subsequently results may be passed to \code{\link{birteLimma}}
}
\usage{
limmaAnalysis(dat, design=NULL, contrast)
}
\arguments{
  \item{dat}{A matrix or ExpressionSet containing the expression values.}
  \item{design}{A design matrix. If no design matrix is provided, the method tries to infer a design matrix by taking column names as a factor. That means the design matrix is constructed according to  ~0 + factor(colnames(dat).}
  \item{contrast}{Contrast for the linear model. Currently only one contrast is supported. If more contrasts are desired separate biRte analyses should be run for each of them.}
}
\value{
Returns a list containing the following entries:
  \item{pvalue.tab}{Containing the result of the topTable function from limma.}
  \item{lm.fit}{Linear fit to the model, i.e. output of function 'eBayes'.}
  \item{design}{The design used.}
  \item{contrast}{The contrasts used.}
}
\references{
G. K. Smyth. Limma : Linear Models for Microarray Data. \emph{Bioinformatics}, (2005):397-420.
}
\examples{
# see birte.run
}
\author{
Holger Froehlich
}
