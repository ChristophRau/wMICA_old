\name{wMICA-package}
\alias{wMICA-package}
\alias{wMICA}
\docType{package}
\title{
wMICA- weighted Maximal Information Component Analysis
}
\description{
co-expression network analysis which combines maximal information (with MINErva) with a modified version of the Interaction Component Modeling for genes (ICMg) functions which have been altered to allow for weighted as opposed to binary gene/gene links.  Provides functions for computing weighted module memberships for each gene in a variety of ways.

This package is built off of the ICMg package and blatantly copies many of its functions.  Many thanks are owed to Juuso Parkkinen for developing the original package.
}
\details{
\tabular{ll}{
Package: \tab wMICA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-11\cr
License: \tab GNU GPL 3 or any later version (at your option)\cr
}
ICMg.links.sampler.SeedsandWeights can be used to sample components for weighted correlation data and/or networks with predefined 'seeds' for one or more modules.  There are two versions of the iteration functions available (R and C) and it is highly recommended to use the faster C versions.
}
\author{
Christoph Rau <chrau@ucla.edu>
Maintainer: Christoph Rau <chrau@ucla.edu>
}
\references{TBA}
\examples{
	## Load data and set parameters
	data(osmo)
	C.boost = 1 ## Use faster C funtions
	Seeds=c()
  Weight=0
  alpha = 10
	beta = 0.01
	B.num = 2
	B.size = 5
	S.num = 2  
	S.size = 5
 	C = 10
	pm0 = 0
	V0 = 1               
	V = 0.1

	## Run sampling
	res = ICMg.lins.sampler.SeedsAndWeights(osmo$ppi, C, Seeds,Weight, alpha, beta, pm0, V0, V, B.num, B.size, S.num, S.size, C.boost) 
	## Compute component membership probabilities for nodes
  	res$comp.memb <- ICMg.get.comp.memberships(osmo$ppi, res)
  	## Compute (hard) clustering for nodes
  	res$clustering <- apply(res$comp.memb, 2, which.max)
}
\keyword{ package }
