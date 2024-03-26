\name{MuSiCEvaluator2scRNAseq}
\alias{MuSiCEvaluator2scRNAseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate MuSiC Using Two scRNAseq Datasets: One as Reference; the Other as Bulk RNAseq
}
\description{
Evaluate MuSiC using two scRNAseq datasets: one as reference; the other as bulk RNAseq.
}
\usage{
MuSiCEvaluator2scRNAseq(
  es4Bulk, 
  esRef, 
  sid, 
  cellType, 
  nBoot = 20, 
  seed = 1234567, 
  nPerm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es4Bulk}{
SingleCellExperiment object of scRNAseq used to construct bulk RNAseq.
}
  \item{esRef}{
SingleCellExperiment object of reference scRNAseq.
}
  \item{sid}{
character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating subject id.
}
  \item{cellType}{
character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating cell type.
}
  \item{nBoot}{
integer. number of bootstrapping.
}
  \item{seed}{
integer. random seed to generate exactly the same results for each call of 'MuSiCEvaluator2scRNAseq'.
}
  \item{nPerm}{
integer. number of permutations to calculate null distributions of mean and sd of pseudo R2.  
  }
}
\details{
The correlation between estimated cell type proportions and real cell type proportions might be confounded by cell type information. We can use partial correlation, instead of Pearson correlation. To calculate partial correlation, we first build a linear regression \eqn{cellType.est ~ cellType + \epsilon_{est}} and use residuals \eqn{\hat{\epsilon}_{est}}  to replace \eqn{cellType.est} (i.e., the estimated cell type proportions for all subjects). We next build another linear regression \eqn{cellType.real ~ cellType + \epsilon_{real}} and use residuals \eqn{\hat{\epsilon}_{real}}  to replace \eqn{cellType.real} (i.e., the real cell type proportions for all subjects). The partial correlation is equal to the Pearson correlation between \eqn{cellType.est} and \eqn{cellType.real}.

Similarly, \code{RMSE.resi} and \code{mAD.resi} are the root mean squared error and mean absolute deviation between \eqn{cellType.est} and \eqn{cellType.real}, respectively.
}
\value{
A list with 16 elements.
  \item{esBulk }{SingleCellExperiment object constructed based on the scRNAseq dataset \code{es4Bulk}. It includes estimated cell type proportions in phenotype data slot so that we can adjust cell type proportions in down-stream analysis.}
  
  \item{prop.real.Ref }{nSubj.Ref x nCellType matrix of true cell type proportions for reference scRNAseq.}
  
  \item{prop.real.Bulk}{nSubj.Bulk x nCellType matrix of true cell type proportion for constructed Bulk RNAseq.}
  
  \item{prop.est}{a list of 2 elements. 
  1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on original reference scRNAseq;
  2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on original reference scRNAseq.}
  
  \item{prop.est.LF}{a list of 2 elements. Long format of the 2 elements in prop.est.
  1st element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on original reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is real cell type proportions;
  3rd column is subject id;
  4th column is cell type.
  
  2nd element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on original reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is real cell type proportions;
  3rd column is subject id;
  4th column is cell type.}
  
  \item{prop.est.LF.resi}{
  
same structure as prop.est.LF. 'resi' indicates the estimated and real cell type proportions
are residuals after regressing out cell types.
  }

  \item{performMat}{a 2x8 matrix of performance indices. 
  rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE (root mean square error), mAD (mean absolute deviation), R (Pearson correlation), R.rank (Spearman rank correlation) before and after regressing out cell types.}
  
  \item{prop.est.boot}{a list of 2 elements. 
  1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on ensembled bootstrappings of reference scRNAseq;
  2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on ensembled bootstrappings of reference scRNAseq.}
  
  \item{prop.est.boot.LF}{a list of 2 elements. Long format of the 2 elements in prop.est.boot.
  1st element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  based on ensembled bootstrappings of reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is real cell type proportions;
  3rd column is subject id;
  4th column is cell type.
  
  2nd element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  based on ensembled bootstrappings of reference scRNAseq.
  1st column is estimated cell tye proportions;
  2nd column is real cell type proportions;
  3rd column is subject id;
  4th column is cell type.
  }
  
  \item{prop.est.LF.resi}{
same structure as prop.est.boot.LF. 'resi' indicates the estimated and real cell type proportions
are residuals after regressing out cell types.
  }
  
  \item{performMat.boot}{a 2x8 matrix of performance indices based on ensembled bootstrappings of reference scRNAseq. 
  rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE (root mean squared error), mAD (mean absolute deviation), R (Pearson correlation), and R.rank (Spearman rank correlation) before and after regressing out cell types.}
  
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
Jinglan Qiu <jinglanqiu1@uvic.ca>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
