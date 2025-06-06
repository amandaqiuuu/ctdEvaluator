# created on June 11, 2021

# calculate sample mean and sd of pseudo R2 after permuting cell types in reference scRNAseq
MuSiCEvaluator2scRNAseq.perm = function(es4Bulk, esRef, sid, cellType, nBoot= 20, seed = 1234567, nPerm = 20)
{
  # cell type in reference scRNAseq
  pDat = colData(esRef)
  cellType.u = unique(pDat[, c(cellType)])
  nSamples = nrow(pDat)
  
  sidVec = NULL
  permVec = NULL
  
  set.seed(seed)  
  # for each permutation
  for(iPerm in 1:nPerm)
  {
    cat(" iPerm=", iPerm)
    # permute cell types
    cellTypeVeci = sample(cellType.u, size = nSamples, replace = TRUE)
    esRefi = esRef
    # upate cell types in reference scRNAseq dataset
    colData(esRefi)[, c(cellType)] = cellTypeVeci
    
    resi = MuSiCEvaluator2scRNAseq.default(es4Bulk = es4Bulk, 
                                           esRef = esRefi, 
                                           sid= sid, 
                                           cellType = cellType, 
                                           nBoot= nBoot, 
                                           seed = seed + iPerm)
    
    R2meanSDMat = calR2MeanSD(obj.Evaluator = resi)
    sidVec = c(sidVec, rownames(R2meanSDMat))
    permVec = c(permVec, rep(iPerm, nrow(R2meanSDMat)))
    if(iPerm > 1)
    {
      R2meanSDMat.LF = rbind(R2meanSDMat.LF, R2meanSDMat)
    } else {
      R2meanSDMat.LF = R2meanSDMat
    }
  }
  
  R2meanSD.frame = data.frame(perm=permVec, 
                              sid = sidVec, 
                              R2.MuSiC.sd = R2meanSDMat.LF[,1],
                              R2.MuSiC.m = R2meanSDMat.LF[,2],
                              R2.NNLS.sd = R2meanSDMat.LF[,3],
                              R2.NNLS.m = R2meanSDMat.LF[,4]
  )
  
  invisible(R2meanSD.frame)
  
  
}

# calculate sample mean and sd of pseudo R2 after permuting cell types in reference scRNAseq
MuSiCEvaluator.perm = function(esBulk, esRef, sid, cellType, nBoot= 20, seed = 1234567, nPerm = 20)
{
  # cell type in reference scRNAseq
  pDat = colData(esRef)
  cellType.u = unique(pDat[, c(cellType)])
  nSamples = nrow(pDat)
  
  sidVec = NULL
  permVec = NULL
  
  set.seed(seed)  
  # for each permutation
  for(iPerm in 1:nPerm)
  {
    cat(" iPerm=", iPerm)
    # permute cell types
    cellTypeVeci = sample(cellType.u, size = nSamples, replace = TRUE)
    esRefi = esRef
    # upate cell types in reference scRNAseq dataset
    colData(esRefi)[, c(cellType)] = cellTypeVeci
    
    resi = MuSiCEvaluator.default(esBulk = esBulk, 
                                           esRef = esRefi, 
                                           sid= sid, 
                                           cellType = cellType, 
                                           nBoot= nBoot, 
                                           seed = seed + iPerm)
    
    R2meanSDMat = calR2MeanSD(obj.Evaluator = resi)
    sidVec = c(sidVec, rownames(R2meanSDMat))
    permVec = c(permVec, rep(iPerm, nrow(R2meanSDMat)))
    if(iPerm > 1)
    {
      R2meanSDMat.LF = rbind(R2meanSDMat.LF, R2meanSDMat)
    } else {
      R2meanSDMat.LF = R2meanSDMat
    }
  }
  
  R2meanSD.frame = data.frame(perm=permVec, 
                              sid = sidVec, 
                              R2.MuSiC.sd = R2meanSDMat.LF[,1],
                              R2.MuSiC.m = R2meanSDMat.LF[,2],
                              R2.NNLS.sd = R2meanSDMat.LF[,3],
                              R2.NNLS.m = R2meanSDMat.LF[,4]
  )
  
  invisible(R2meanSD.frame)
  
  
}

