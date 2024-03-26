\name{es4Bulk.XinT2D.sce}
\alias{es4Bulk.XinT2D.sce}
\docType{data}
\title{
XinT2D scRNAseq Dataset to be Constructed as Bulk RNAseq Dataset
}
\description{
ExpressionSet object stores gene expression levels for scRNAseq dataset XinT2D.
}
\usage{data("es4Bulk.XinT2D.sce")}
\details{
The dataset \code{es4Bulk.XinT2D.sce} is an ExpressionSet object
storing a scRNAseq dataset to be used to construct Bulk RNAseq with known cell type proportions. It contains 39849 genes and 1492 samples (i.e., cells).
The 1492 cells belong to 6 T2D subjects and 12 non-T2D subjects.
Among the 1492 cells, there are 4 cell types: alpha (\eqn{n=886} cells), beta (\eqn{n=472} cells), delta (\eqn{n=49} cells), and gamma (\eqn{n=85} cells). 
}
\source{
This dataset is listed in Wang et al.'s tutorial website \url{https://xuranw.github.io/MuSiC/articles/MuSiC.html}.
}
\references{
Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat Commun 10, 380 (2019). https://doi.org/10.1038/s41467-018-08023-x
}
\examples{
data(es4Bulk.XinT2D.sce)

print(es4Bulk.XinT2D.sce)

# unique subject ids
print(unique(es4Bulk.XinT2D.sce$SubjectName))

# frequencies of cell types
print(table(es4Bulk.XinT2D.sce$cellType))

}
\keyword{datasets}
