\name{wMICA.Run}
\alias{wMICA.Run}
\title{ wMICA: weighted Maximal Information Component Analysis }
\description{
  wMICA.Run is the core function of the wMICA package and runs the entire package to generate weighted module memberships using the MINE algorithm and a modified weighted version of the ICMg algorithm.
}
\usage{
wMICA.Run(input,output_prefix="MICA",n.cores=detectCores()-1,
threshold=.4,C=20,Seeds=c(),Weight=0,alpha=10,beta=.01,
B.num=10,B.size=10,S.num=10,S.size=10,MIC_ARRAY_FILE=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ A NxM matrix where N are samples and M are genes  }
  \item{output_prefix}{ prefix to affix to outputted files }
  \item{n.cores}{Number of cores of computer to use when calculating MINE score.  Defaults to one fewer than the total number of cores on the machine}
  \item{threshold}{lowest value edge to consider.  A lower threshold (<=.3) is more accurate, but takes significantly longer}
  \item{C}{The Number of Modules}
  \item{Seeds}{a Kx2 matrix (where K is the number of genes to be seeded into a module) with the first column containing the gene number (as in L) and the second column being the desired module to be seeded into}
  \item{Weight}{the number of iterations to weight each seeded gene by}
  \item{alpha}{ Hyperparameter describing the global distribution over components, larger alpha gives a more uniform distribution. }
  \item{beta}{ Hyperparameter describing the component-wise distributions over nodes, larger beta gives a more uniform distribution. }
  \item{B.num}{ Number of burnin rounds.* }
  \item{B.size}{ Size of one burnin round.* }
  \item{S.num}{ Number of sample rounds.* }
  \item{S.size}{ Size of one sample round.* }
  \item{MIC_ARRAY_FILE}{ Will accept an already created MIC array (if you want to pick a different threshold, for instance) }
}
\details{
  * One run involves three components.  First, the Maximal Information score between each pair of genes is calculated.  Then, Edges with a MIC score greater than the threshold are compiled and given to the module identification algorithm.  Modules are determined using ICMg, whose details follow:
  
  * One run consists of two parts, during burnin the sampler is expected to mix, after which the samples are taken. Information about convergence (convN and convL are estimates of convergence for link and node sampling, respectively) and component sizes are printed after each burnin/sample round. For example: B.num=8, B.size=100, S.num=20, S.size=10, runs 800 burnin iterations in 8 rounds and then takes 20 samples with an interval of 10 iterations.
  
  *Finally, Module Memberships are extracted and all files (MIC array, Links, Memberships) are saved.
}
\value{
  Returns 3 output files to the working directory with prefix X:
  \item{X_MICARRAY.csv}{ N*N matrix of Maximal Information Scores }
  \item{X_Links.csv}{ All included links used in wMICA } 
  \item{X_MODULEMEMBERSHIPS.csv}{C*N matrix of results, where the rows are the modules and the columns are the genes } 
}
\author{ Christoph Rau }
\seealso{}
\examples{
	## Load data and set parameters
	data(MICASample)
	Seeds=c()
  Weight=0
  C=20
  wMICA.Run(indata,Seeds=Seeds,Weight=Weight,C=C)
}
\keyword{ methods }
