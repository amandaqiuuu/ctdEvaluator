# v3 created on Jan. 9, 2021
#  (1) only output R2.NNLS; R2.MuSiC = R.squared
#
# v2 created on Jan. 8, 2021
#  (1) output weighted adjusted R2 'adjR2.weight'
#  (2) output mean deviance (sum of square of weighted residuals)

music.basic2 = function (Y, X, S, Sigma, iter.max, nu, eps) 
{
  k = ncol(X)
  nGenes = length(Y)

  lm.D = nnls(X, Y)
  r = resid(lm.D)

  # adjusted R square
  var.res.NNLS=var(r, na.rm = TRUE)
  var.total.NNLS=var(Y, na.rm = TRUE)
  R2.NNLS = 1 - var.res.NNLS/var.total.NNLS

  weight.gene = 1/(nu + r^2 + colSums((lm.D$x * S)^2 * t(Sigma)))
  Y.weight = Y * sqrt(weight.gene)
  D.weight = sweep(X, 1, sqrt(weight.gene), "*")
  lm.D.weight = nnls(D.weight, Y.weight)
  p.weight = lm.D.weight$x/sum(lm.D.weight$x)
  p.weight.iter = p.weight
  r = resid(lm.D.weight)
  for (iter in 1:iter.max) {
    weight.gene = 1/(nu + r^2 + colSums((lm.D.weight$x * 
      S)^2 * t(Sigma)))
    Y.weight = Y * sqrt(weight.gene)
    D.weight = X * as.matrix(sqrt(weight.gene))[, rep(1, 
      k)]
    lm.D.weight = nnls(D.weight, Y.weight)
    p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
    r.new = resid(lm.D.weight)
    if (sum(abs(p.weight.new - p.weight)) < eps) {
      p.weight = p.weight.new
      r = r.new
      R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)

      fitted = X %*% as.matrix(lm.D.weight$x)
      var.p = diag(solve(t(D.weight) %*% D.weight)) * 
        mean(r^2)/sum(lm.D.weight$x)^2
      return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, 
        fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), 
        p.weight = p.weight, q.weight = lm.D.weight$x, 
        fit.weight = fitted, resid.weight = Y - X %*% 
          as.matrix(lm.D.weight$x), weight.gene = weight.gene, 
        converge = paste0("Converge at ", iter), rsd = r, 
        R.squared = R.squared, 
	R2.NNLS = R2.NNLS, 
	var.p = var.p))
    }
    p.weight = p.weight.new
    r = r.new
  }
  fitted = X %*% as.matrix(lm.D.weight$x)
  R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)

  var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2

  return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x, 
    fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), p.weight = p.weight, 
    q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight = Y - 
      X %*% as.matrix(lm.D.weight$x), weight.gene = weight.gene, 
    converge = "Reach Maxiter", rsd = r, R.squared = R.squared, 
    R2.NNLS = R2.NNLS, 
    var.p = var.p))
}

