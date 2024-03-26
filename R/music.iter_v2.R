# v2 created on Jan. 8, 2021
#  (1) call music.basic2

music.iter2 = function (Y, D, S, Sigma, iter.max = 1000, nu = 1e-04, eps = 0.01, 
  centered = FALSE, normalize = FALSE) 
{
  if (length(S) != ncol(D)) {
    common.cell.type = intersect(colnames(D), names(S))
    if (length(common.cell.type) <= 1) {
      stop("Not enough cell types!")
    }
    D = D[, match(common.cell.type, colnames(D))]
    S = S[match(common.cell.type, names(S))]
  }
  if (ncol(Sigma) != ncol(D)) {
    common.cell.type = intersect(colnames(D), colnames(Sigma))
    if (length(common.cell.type) <= 1) {
      stop("Not enough cell type!")
    }
    D = D[, match(common.cell.type, colnames(D))]
    Sigma = Sigma[, match(common.cell.type, colnames(Sigma))]
    S = S[match(common.cell.type, names(S))]
  }
  k = ncol(D)
  common.gene = intersect(names(Y), rownames(D))
  if (length(common.gene) < 0.1 * min(length(Y), nrow(D))) {
    stop("Not enough common genes!")
  }
  Y = Y[match(common.gene, names(Y))]
  D = D[match(common.gene, rownames(D)), ]
  Sigma = Sigma[match(common.gene, rownames(Sigma)), ]
  X = D
  if (centered) {
    X = X - mean(X)
    Y = Y - mean(Y)
  }
  if (normalize) {
    X = X/sd(as.vector(X))
    S = S * sd(as.vector(X))
    Y = Y/sd(Y)
  }
  else {
    Y = Y * 100
  }
  lm.D = music.basic2(Y, X, S, Sigma, iter.max = iter.max, 
    nu = nu, eps = eps)
  return(lm.D)
}

