\name{esRef.EMTAB.fakeCT.sce}
\alias{esRef.EMTAB.fakeCT.sce}
\docType{data}
\title{
Revised EMTAB scRNAseq Dataset Can be Used as Reference scRNAseq in Cell Type Deconvolution wit Fake Cell Types
}
\description{
ExpressionSet object stores gene expression levels for revised scRNAseq dataset EMTAB.
}
\usage{data("esRef.EMTAB.fakeCT.sce")}
\details{
The cell types
acinar, ductal, PSC, unclassified endocrine were relabelled as gamma, delta, alpha, and beta, respectively.
Poor results would be resulted in when we use \code{esRef.EMTAB.fakeCT.sce} as reference to deconvolute cell types alpha, beta, delta, and gamma.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(esRef.EMTAB.fakeCT.sce)

print(esRef.EMTAB.fakeCT.sce)

# unique subject ids
print(unique(esRef.EMTAB.fakeCT.sce$SubjectName))

# frequencies of true cell types
print(table(as.character(esRef.EMTAB.fakeCT.sce$cellType.true)))

# frequencies of relabeled cell types
print(table(esRef.EMTAB.fakeCT.sce$cellType))

# cross table

print(table(as.character(esRef.EMTAB.fakeCT.sce$cellType.true), esRef.EMTAB.fakeCT.sce$cellType))
}
\keyword{datasets}
