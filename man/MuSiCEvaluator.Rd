\name{MuSiCEvaluator}
\alias{MuSiCEvaluator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate Cell Type Proportions of Bulk RNAseq Dataset Estimated via the MuSiC Method
}
\description{
Evaluate cell type proportions of bulk RNAseq dataset estimated via the MuSiC method.
}
\usage{
MuSiCEvaluator(
  esBulk, 
  esRef, 
  sid, 
  cellType, 
  nBoot = 20, 
  seed = 1234567, 
  nPerm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{esBulk}{
SingleCellExperiment object of bulk RNAseq with unknown cell type info.
}
  \item{esRef}{
SingleCellExperiment object of reference scRNAseq.
}
  \item{sid}{
character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating subject id.
}
  \item{cellType}{
character. column variable name for phenotype data in 'esRef' indicating cell type.
}
  \item{nBoot}{
integer. number of bootstrapping
}
  \item{seed}{
integer. random seed to generate exactly the same results for each call of 'MuSiCEvaluator'.
}
  \item{nPerm}{
integer. number of permutations to calculate null distributions of mean and sd of pseudo R2.  
  }
}
\details{
Pseudo \eqn{R^2} and its stability are used to evaluate cell type proportions estimated via MuSiC.
\eqn{R^2} and its stability are used to evaluate cell type proportions estimated via NNLS.
Stability is based on the bootstrapping reference scRNAseq dataset.

Esemble estimated cell type proportions are used to improve the accuracy of the estimation.

}
\value{
A list of 11 elements.
  \item{esBulk }{updated SingleCellExperiment object including estimated cell type proportions in phenotype data slot so that we can adjust cell type proportions in down-stream analysis}
  
  \item{prop.real.Ref }{nSubj.Ref x nCellType matrix of true cell type proportions for reference scRNAseq}

  \item{prop.real.Ref.LF}{[nSubj.Ref x nCellType] x 3 data frame of long format of prop.real.Ref.
  1st column is real cell type proportions for all subjects (subject by subject);
  2nd column is subject id;
  3rd column is cell type.}

  \item{prop.est}{a list of 2 elements. 
  1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on original reference scRNAseq;
  2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on original reference scRNAseq.}
  
  \item{prop.est.LF}{a list of 2 elements. Long format of the 2 elements in prop.est.
  1st element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on original reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is subject id;
  3rd column is cell type.
  
  2nd element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on original reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is subject id;
  3rd column is cell type.  }

  \item{prop.est.boot}{a list of 2 elements. 
  1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on ensembled bootstrappings of reference scRNAseq;
  2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on ensembled bootstrappings of reference scRNAseq}
  
  \item{prop.est.boot.LF}{a list of 2 elements. Long format of the 2 elements in prop.est.boot.
  1st element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on ensembled bootstrappings of reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is subject id;
  3rd column is cell type.
  
  2nd element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on ensembled bootstrappings of reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is subject id;
  3rd column is cell type.}
  
  \item{R2.MuSiCFrame.boot}{long format pseudo R2 via MuSiC used to draw trajectory of R2 across bootstrapping with 3 columns.
  1st column is R2.MuSiC; 2nd column indicates if results are based on original reference scRNAseq or
  based on which bootstrapping reference scRNAseq.}
  
  \item{R2.NNLSFrame.boot}{long format pseudo R2 via MuSiC used to draw trajectory of R2 across bootstrapping with 3 columns.
  1st column is R2.MuSiC; 2nd column indicates if results are based on original reference scRNAseq or
  based on which bootstrapping reference scRNAseq.}
  
  
  \item{R2meanSDMat}{\code{nSubj x 4} data matrix. Rows are subjects. 
  1st column is the sample sd of pseudo R2 (MuSiC) for each subject based on bootstrapping. 
  2nd column is the sample mean of pseudo R2 (MuSiC) for each subject based on bootstrapping. 
  3rd column is the sample sd of R2 (NNLS) for each subject based on bootstrapping. 
  4th column is the sample mean of R2 (NNLS) for each subject based on bootstrapping. 
  }
  
  \item{R2meanSD.frame}{a data frame with 6 columns. 
  1st column indicates permutation.
  2nd column indicates subject id.
  3rd column is the sample sd of pseudo R2 (MuSiC) for each subject based on bootstrapping for each permutation. 
  4th column is the sample mean of pseudo R2 (MuSiC) for each subject based on bootstrapping for each permutation. 
  5th column is the sample sd of R2 (NNLS) for each subject based on bootstrapping for each permutation.
  6th column is the sample mean of R2 (NNLS) for each subject based on bootstrapping for each permutation.
  }

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Weiliang Qiu <Weiliang.Qiu@Sanofi.com> \cr 
Mindy Zhang <Mindy.Zhang@sanofi.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(SummarizedExperiment)
library(SingleCellExperiment)

# load bulk RNAseq data
data(esBulk.XinT2D.sce)
print(esBulk.XinT2D.sce)

# load reference scRNAseq data
data(esRef.EMTAB.sce)
print(esRef.EMTAB.sce)

# First do cell type deconvolution, then evaluate the estimates.
# For testing purpose, we set number of bootstrapping 'nBoot=3'.
# In real applications, please set nBoot >= 20.
# Since it is computationally intensive, it would take a long running time 
#   when 'nBoot' is too large.
res = MuSiCEvaluator(esBulk = esBulk.XinT2D.sce, esRef = esRef.EMTAB.sce, 
                              sid = "SubjectName", 
                              cellType = "cellType", 
                              nBoot = 3, 
                              seed = 1234567,
                              nPerm = NULL)
print(names(res))

##############
# trajectories of pseudo R2 across different reference scRNAseq datasets;
# cell type proportions estimated via MuSiC 
# note that pseudo R2 can be negative
library(iPlotFun)
R2.MuSiC.boot.LF3 = res$R2.MuSiCFrame.boot

g = trajPlot1Grp(dat = R2.MuSiC.boot.LF3, 
         time = "grp", 
         y = "R2.MuSiC", 
         sid = "sid", 
         numericTime = FALSE,
         order.time = NULL, 
         xlab = "", 
         ylab = "pseudo R2", 
         title = "MuSiC: esBulk.XinT2D and esRef.EMTAB", 
         xAngle = 45
)
print(g)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ method }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
