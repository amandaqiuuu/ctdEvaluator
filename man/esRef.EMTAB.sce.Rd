\name{esRef.EMTAB.sce}
\alias{esRef.EMTAB.sce}
\docType{data}
\title{
EMTAB scRNAseq Dataset Can be Used as Reference scRNAseq in Cell Type Deconvolution
}
\description{
ExpressionSet object stores gene expression levels for scRNAseq dataset EMTAB.
}
\usage{data("esRef.EMTAB.sce")}
\details{
The dataset \code{esRef.EMTAB.sce} is an ExpressionSet object
storing a scRNAseq dataset to be used as a reference scRNAseq in MuSiC algorithm. It contains 25453 genes and 748 samples (i.e., cells).
The 748 cells belong to 6 non-T2D subjects.
Among the 748 cells, there are 4 cell types: alpha (\eqn{n=443} cells), beta (\eqn{n=171} cells), delta (\eqn{n=59} cells), and gamma (\eqn{n=75} cells). 
}
\source{
This dataset is a subset of EMTAB scRNAseq dataset listed in Wang et al.'s tutorial website \url{https://xuranw.github.io/MuSiC/articles/MuSiC.html}.
}
\references{
Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat Commun 10, 380 (2019). https://doi.org/10.1038/s41467-018-08023-x
}
\examples{
data(esRef.EMTAB.sce)

print(esRef.EMTAB.sce)

# unique subject ids
print(unique(esRef.EMTAB.sce$SubjectName))

# frequencies of cell types
print(table(esRef.EMTAB.sce$cellType))
}
\keyword{datasets}
