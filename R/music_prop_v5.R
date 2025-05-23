# modified on March 30, 2021
#  (1) move definitions of functions genRefscRNAseqData* to file 'genRefscRNAseqData.R'

# v5 created on Jan. 8, 2021
#  (1) return R2.NNLS 

# v4 created on Jan. 4, 2021
#  (1) add function 'genRefscRNAseqData.default'
#  (2) add function 'genRefscRNAseqData.default2'
#  (3) add function 'genRefscRNAseqData2'
#
# v3 created on Jan. 3, 2021
#  (1) sub-sampling to generate 'nSim' scRNAseq datasets that similar to the original 
#      reference scRNAseq datasets
#  (2) obtain MLE of the cell type proportions based on 'nSim' estimated cell type proportions
#  (3) remove permFlag
#
# v2 created on Jan. 2, 2021
#  (1) output residual sum of square (i.e., deviance)
#  (2) add 'permFlag' option

###################


#
#bulk.eset	#ExpressionSet for bulk data
#sc.eset	#ExpressionSet for single cell data
#markers	#vector or list of gene names, default as NULL. If NULL, use all genes that provided by both bulk and single cell dataset.
#clusters	#character, the phenoData of single cell dataset used as clusters;
#samples	#character,the phenoData of single cell dataset used as samples;
#select.ct	#vector of cell types, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#cell_size	#data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#ct.cov	#logical. If TRUE, use the covariance across cell types;
#verbose	#logical, default as TRUE.
#iter.max	#numeric, maximum iteration number
#nu	#regulation parameter, take care of weight when taking recipical
#eps	#Thredshold of convergence
#centered	#logic, substract avg of Y and D
#normalize	#logic, divide Y and D by their standard deviation
#
music_prop3 = function (bulk.eset, sc.eset, markers = NULL, clusters, samples, 
  select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE, 
  iter.max = 1000, nu = 1e-04, eps = 0.01, centered = FALSE, 
  normalize = FALSE, ...) 
{
  bulk.gene = rownames(bulk.eset)[base::rowMeans(SingleCellExperiment::counts(bulk.eset)) != 
    0]
  bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
  if (is.null(markers)) {
    sc.markers = bulk.gene
  }
  else {
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  
  sc.basis = MuSiC::music_basis(sc.eset, non.zero = TRUE, markers = sc.markers, 
    clusters = clusters, samples = samples, select.ct = select.ct, 
    cell_size = cell_size, ct.cov = ct.cov, verbose = verbose)
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  if (is.null(markers)) {
    if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.eset))) 
      stop("Too few common genes!")
  }
  else {
    if (length(cm.gene) < 0.2 * length(unlist(markers))) 
      stop("Too few common genes!")
  }
  if (verbose) {
    message(paste("Used", length(cm.gene), "common genes..."))
  }
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
  m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]
  M.S = colMeans(sc.basis$S, na.rm = T)
  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }
    else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    }
    else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 
      1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]), 
      ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  Yjg = MuSiC::relative.ab(SingleCellExperiment::counts(bulk.eset)[m.bulk, ])
  N.bulk = ncol(bulk.eset)
  if (ct.cov) {
    Sigma.ct = sc.basis$Sigma.ct[, m.sc]
    if (sum(Yjg[, i] == 0) > 0) {
      D1.temp = D1[Yjg[, i] != 0, ]
      Yjg.temp = Yjg[Yjg[, i] != 0, i]
      Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
      if (verbose) 
        message(paste(colnames(Yjg)[i], "has common genes", 
          sum(Yjg[, i] != 0), "..."))
    }
    else {
      D1.temp = D1
      Yjg.temp = Yjg[, i]
      Sigma.ct.temp = Sigma.ct
      if (verbose) 
        message(paste(colnames(Yjg)[i], "has common genes", 
          sum(Yjg[, i] != 0), "..."))
    }
    lm.D1.weighted = MuSiC::music.iter.ct(Yjg.temp, D1.temp, M.S, 
      Sigma.ct.temp, iter.max = iter.max, nu = nu, eps = eps, 
      centered = centered, normalize = normalize)
    Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
    Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
    weight.gene.temp = rep(NA, nrow(Yjg))
    weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
    Weight.gene = cbind(Weight.gene, weight.gene.temp)
    r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
    Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
  }
  else {
    Sigma = sc.basis$Sigma[m.sc, ]
    valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) == 
      0) & (!is.na(M.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    if (verbose) {
      message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    }
    D1 = D1[, valid.ct]
    M.S = M.S[valid.ct]
    Sigma = Sigma[, valid.ct]
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    R2.NNLS = NULL

    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp = D1[Yjg[, i] != 0, ]
        Yjg.temp = Yjg[Yjg[, i] != 0, i]
        Sigma.temp = Sigma[Yjg[, i] != 0, ]
        if (verbose) 
          message(paste(colnames(Yjg)[i], "has common genes", 
            sum(Yjg[, i] != 0), "..."))
      }
      else {
        D1.temp = D1
        Yjg.temp = Yjg[, i]
        Sigma.temp = Sigma
        if (verbose) 
          message(paste(colnames(Yjg)[i], "has common genes", 
            sum(Yjg[, i] != 0), "..."))
      }
      lm.D1.weighted = music.iter2(Yjg.temp, D1.temp, M.S, 
        Sigma.temp, iter.max = iter.max, nu = nu, eps = eps, 
        centered = centered, normalize = normalize)

      R2.NNLS = c(R2.NNLS, lm.D1.weighted$R2.NNLS)

      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)

  names(R2.NNLS)=colnames(Yjg)

  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene, 
    Weight.gene = Weight.gene, r.squared.full = r.squared.full, 
    Var.prop = Var.prop, 
    R2.NNLS = R2.NNLS
    ))
}

