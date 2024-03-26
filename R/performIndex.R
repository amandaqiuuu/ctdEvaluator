# calculate performance indices, including RMSE, mAD, Pearson correlation,
#   and partial correlation between estimated and real cell type proportions;
#   also calculate root mean square error and mean absolute deviation after regressing out cell types

# prop.LF - a data frame object having 4 elements
#  1st element: prop.est (a vectorized estimated cell type proportions subject by subject);
#  2nd element: prop.real(a vectorized real cell type proportions subject by subject)
#  3rd element: sid (subject ids repeated nCellType times)
#  4th element: cellType (cell type names, each repeat nSubj times)

performIndexFunc.default = function(prop.LF)
{
  ######
  # calculate RMSE, mAD, and R based on cell type proportions
  ######
  # original data
  myx = prop.LF$prop.est
  myy = prop.LF$prop.real
  
  # RMSE
  RMSE = sqrt(mean( (myx - myy)^2  ))
  
  # mAD
  mAD = mean( abs(myx - myy)  ) 
  
  # partial correlation
  R = cor(x=myx, y=myy)
  
  R.rank = cor(x=myx, y=myy, method="spearman")

  ################################
  # regress out effect of cell type
  ################################
  prop.LF.resi = prop.LF
  tt.est = lm(prop.est ~ factor(cellType), data=prop.LF.resi)
  tt.real = lm(prop.real ~ factor(cellType), data=prop.LF.resi)
  
  # replace prop.est and prop.real by residuals
  prop.LF.resi$prop.est = tt.est$residuals
  prop.LF.resi$prop.real = tt.real$residuals
  
  ######
  # calculate RMSE, mAD, and R based on cell type proportions after regression out cell types
  ######
  # MuSiC
  # Pearson correlation
  myx = prop.LF.resi$prop.est
  myy = prop.LF.resi$prop.real
  
  # RMSE
  RMSE.resi = sqrt(mean( (myx - myy)^2  ))
  
  # mAD
  mAD.resi = mean( abs(myx - myy)  ) 
  
  # partial correlation
  R.resi = cor(x=myx, y=myy)
  R.rank.resi = cor(x=myx, y=myy, method="spearman")

  performIndex = c(RMSE, mAD, R, R.rank, RMSE.resi, mAD.resi, R.resi, R.rank.resi)  
  names(performIndex) = c("RMSE", "mAD", "R", "R.rank", "RMSE.resi", "mAD.resi", "R.resi", "R.rank.resi")
  
  res = list(performIndex = performIndex,
             prop.LF.resi = prop.LF.resi)
  
  invisible(res)
}


performIndexFunc = function(prop.est, prop.real.Bulk)
{
  ####
  # get long format estimated cell type proportions
  ####
  prop.est.MuSiC.LF = data.frame(prop.est=as.vector(prop.est$MuSiC), 
                                 prop.real=as.vector(prop.real.Bulk),
                                 sid = rep(rownames(prop.est$MuSiC), 
                                           ncol(prop.est$MuSiC)),
                                 cellType = rep(colnames(prop.est$MuSiC),
                                                each = nrow(prop.est$MuSiC)))
  
  prop.est.NNLS.LF = data.frame(prop.est=as.vector(prop.est$NNLS), 
                                prop.real=as.vector(prop.real.Bulk),
                                sid = rep(rownames(prop.est$NNLS), 
                                          ncol(prop.est$NNLS)),
                                cellType = rep(colnames(prop.est$NNLS),
                                               each = nrow(prop.est$NNLS)))
  
  prop.est.LF = list(MuSiC = prop.est.MuSiC.LF,
                     NNLS = prop.est.NNLS.LF)
  
  pIndex.MuSiC = performIndexFunc.default(prop.LF = prop.est.LF$MuSiC)
  pIndex.NNLS = performIndexFunc.default(prop.LF = prop.est.LF$NNLS)
  
  ################################
  # calculate RMSE, mAD, and R
  # RMSE: root mean squared error; mAD: mean absolute deviation; R: Pearson correlation
  # performMat is a 2x8 matrix. rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE, mAD, R, R.rank,
  #  RMSE.resi, mAD.resi, R.resi, R.rank.resi
  performMat = rbind(pIndex.MuSiC$performIndex, pIndex.NNLS$performIndex)
  rownames(performMat) = c("MuSiC", "NNLS")
  
  # a data frame containing residuals after regressing out cell types
  prop.MuSiC.LF.resi = pIndex.MuSiC$prop.LF.resi
  prop.NNLS.LF.resi = pIndex.NNLS$prop.LF.resi
  
  prop.est.LF.resi = list(MuSiC = prop.MuSiC.LF.resi,
                          NNLS = prop.NNLS.LF.resi)
  
  res = list(performMat = performMat, 
             prop.est.LF = prop.est.LF,
             prop.est.LF.resi = prop.est.LF.resi
             )
  
  invisible(res)
}