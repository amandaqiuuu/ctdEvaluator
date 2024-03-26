\name{esBulk.XinT2D.sce}
\alias{esBulk.XinT2D.sce}
\docType{data}
\title{
Bulk RNAseq Dataset Constructed based on XinT2D
}
\description{
SingleCellExperiment object stores gene expression levels for bulk RNAseq dataset constructed from scRNAseq dataset XinT2D.
}
\usage{data("esBulk.XinT2D.sce")}
\details{
The dataset \code{esBulk.XinT2D.sce} is an SingleCellExperiment object
storing a Bulk RNAseq dataset that was constructed by the scRNAseq dataset \code{es4Bulk.XinT2D.sce} by averaging gene expression across all cells for each subject. The cell type proportions are hence known for \code{esBulk.XinT2D.sce}. It contains 39849 genes and 18 samples (i.e., subjects).
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat Commun 10, 380 (2019). https://doi.org/10.1038/s41467-018-08023-x
}
\examples{
library(SummarizedExperiment)
library(SingleCellExperiment)

data(esBulk.XinT2D.sce)

print(esBulk.XinT2D.sce)

# unique subject ids
print(unique(esBulk.XinT2D.sce$SubjectName))

# true cell type proportions
colDat = SummarizedExperiment::colData(esBulk.XinT2D.sce)
print(dim(colDat))
print(colDat)
}
\keyword{datasets}
