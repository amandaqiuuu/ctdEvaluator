# ctdEvaluator

**ctdEvaluator**: Functions to Evaluate Results of Cell Type Deconvolution Algorithms

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](LICENSE)

---

## About

**ctdEvaluator** is an R package providing a collection of functions to evaluate the results of cell type deconvolution, particularly for bulk RNAseq datasets using single-cell RNAseq (scRNAseq) as reference. The package implements novel evaluation methods for cell type proportion estimation, supporting robust benchmarking and application of deconvolution algorithms such as MuSiC.

This package is associated with a manuscript to be published soon. Please cite the manuscript when using this package in your research (citation details will be updated upon publication).

---

## Installation

### Prerequisites
- R (>= 3.6.0)
- [Bioconductor](https://bioconductor.org/)

### Required Packages
- Biobase
- MuSiC ([GitHub: xuranw/MuSiC](https://github.com/xuranw/MuSiC))
- xbioc
- nnls
- knitr, rmarkdown
- ggplot2
- iPlotFun (contact Weiliang Qiu for access)
- SingleCellExperiment, SummarizedExperiment

### Install ctdEvaluator
```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biobase", "SingleCellExperiment", "SummarizedExperiment"))

# Install MuSiC from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('xuranw/MuSiC')

# Install ctdEvaluator (local source)
# setwd('path/to/ctdEvaluator')
# install.packages(".", repos = NULL, type = "source")
```

---

## Included Datasets

The package provides four datasets for demonstration and benchmarking:

- **esRef.EMTAB.sce**: Reference scRNAseq (SingleCellExperiment) from EMTAB, 748 cells, 4 cell types (alpha, beta, delta, gamma).
- **esRef.EMTAB.fakeCT.sce**: Reference scRNAseq with relabeled cell types, 299 cells, 4 relabeled types.
- **es4Bulk.XinT2D.sce**: scRNAseq for constructing bulk RNAseq, 1492 cells, 4 cell types.
- **esBulk.XinT2D.sce**: Bulk RNAseq constructed from es4Bulk.XinT2D, 18 samples, known cell type proportions.

See the [vignette](vignettes/ctdEvaluator.Rmd) for detailed descriptions and usage examples.

---

## Main Functions

### `MuSiCEvaluator`
Estimate cell type proportions in bulk RNAseq using a reference scRNAseq dataset. Suitable for real applications where bulk cell type proportions are unknown.

**Usage:**
```r
MuSiCEvaluator(esBulk, esRef, sid, cellType, nBoot = 20, seed = 1234567)
```
- `esBulk`: ExpressionSet/SummarizedExperiment of bulk RNAseq
- `esRef`: ExpressionSet/SummarizedExperiment of reference scRNAseq
- `sid`: Column name for subject ID
- `cellType`: Column name for cell type
- `nBoot`: Number of bootstraps (default 20)
- `seed`: Random seed

### `MuSiCEvaluator2scRNAseq`
Evaluate deconvolution performance using two scRNAseq datasets (one as reference, one to construct bulk RNAseq with known proportions).

**Usage:**
```r
MuSiCEvaluator2scRNAseq(es4Bulk, esRef, sid, cellType, nBoot = 20, seed = 1234567)
```
- `es4Bulk`: scRNAseq used to construct bulk RNAseq
- `esRef`: Reference scRNAseq
- `sid`, `cellType`, `nBoot`, `seed`: as above

---

## Example Analysis

See the [vignette](vignettes/ctdEvaluator.Rmd) for a full workflow, including:
- Loading example datasets
- Running deconvolution and evaluation
- Interpreting results and visualizations

**Quick start:**
```r
library(ctdEvaluator)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Load example data
data(esBulk.XinT2D.sce)
data(esRef.EMTAB.sce)

# Estimate cell type proportions
res <- MuSiCEvaluator(esBulk.XinT2D.sce, esRef.EMTAB.sce, sid = "SubjectName", cellType = "cellType")
```

---

## License

This package is distributed under the [GNU General Public License v2.0](LICENSE).

---

## Contact & Acknowledgments

- Maintainer: Jinglan Qiu (<jinglanqiu1@uvic.ca>)
- Contributors: Weiliang Qiu, Xuekui Zhang
- For iPlotFun package access: Weiliang Qiu (<Weiliang.Qiu@sanofi.com>)

We thank the authors of MuSiC and the Bioconductor community for foundational tools and datasets.

---

## Changelog

See [NEWS](NEWS) for version history and recent changes.

---

## Citation

If you use **ctdEvaluator** in your research, please cite the associated manuscript (details to be updated upon publication) and this package.
