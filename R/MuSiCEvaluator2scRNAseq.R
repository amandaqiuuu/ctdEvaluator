# created on June 11, 2021
#  evaluate MuSiC using 2 scRNAseq datasets: one as reference; the other as bulk RNAseq

# es4Bulk - ExpressionSet object of scRNAseq used to construct bulk RNAseq
# esRef - ExpressionSet object of reference scRNAseq
# sid - character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating subject id
# cellType - character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating cell type
# nBoot - integer. number of bootstrapping
# seed - integer. random seed to generate exactly the same results for each call of 'MuSiCEvaluator2scRNAseq'
MuSiCEvaluator2scRNAseq.default = function(es4Bulk, esRef, sid, cellType, nBoot= 20, seed = 1234567)
{

  pDatBulk = colData(es4Bulk)
  sidBulk = pDatBulk[, c(sid)]
  if(is.null(sidBulk))
  {
    stop(paste("Error! 'es4Bulk' does not contain phenotype variable", sid, sep=""))
  }
  cellTypeBulk = pDatBulk[, c(cellType)]
  if(is.null(cellTypeBulk))
  {
    stop(paste("Error! 'es4Bulk' does not contain phenotype variable", cellType, sep=""))
  }
  
  pDatRef = colData(esRef)
  sidRef = pDatRef[, c(sid)]
  if(is.null(sidRef))
  {
    stop(paste("Error! 'esRef' does not contain phenotype variable", sid, sep=""))
  }
  cellTypeRef = pDatRef[, c(cellType)]
  if(is.null(cellTypeRef))
  {
    stop(paste("Error! 'esRef' does not contain phenotype variable", cellType, sep=""))
  }
  
  
  # reference scRNAseq
  esRef$subj=colData(esRef)[, c(sid)]
  esRef$cellType = factor(colData(esRef)[, c(cellType)])
  esRef$cellTypeID = esRef$cellType

  # obtain cell type proportions for esRef
  esRef.construct.full = MuSiC::bulk_construct(esRef, clusters = 'cellType', samples = 'subj')

  # calculate cell type proportions
  # rows are subjects; columns are cell types
  prop.real.Ref = data.matrix(MuSiC::relative.ab(esRef.construct.full$num.real, by.col = FALSE))

  # long format of prop.real.Ref
  prop.real.Ref.LF = data.frame(prop.real = as.vector(prop.real.Ref),
                                sid = rep(rownames(prop.real.Ref), ncol(prop.real.Ref)),
                                cellType = rep(colnames(prop.real.Ref), each = nrow(prop.real.Ref))
  )
  
  ##############

  # scRNAseq for constructing bulk RNAseq
  es4Bulk$subj=colData(es4Bulk)[, c(sid)]
  es4Bulk$cellType = factor(colData(es4Bulk)[, c(cellType)])
  es4Bulk$cellTypeID = es4Bulk$cellType
  
  # obtain cell type proportions for esRef
  es4Bulk.construct.full = MuSiC::bulk_construct(es4Bulk, clusters = 'cellType', samples = 'subj')

  # calculate cell type proportions
  # rows are subjects; columns are cell types
  prop.real.Bulk = data.matrix(MuSiC::relative.ab(es4Bulk.construct.full$num.real, by.col = FALSE))

  # long format of prop.real.Bulk
  prop.real.Bulk.LF = data.frame(prop.real = as.vector(prop.real.Bulk),
                                sid = rep(rownames(prop.real.Bulk), ncol(prop.real.Bulk)),
                                cellType = rep(colnames(prop.real.Bulk), each = nrow(prop.real.Bulk))
  )
  
  
  # constructed Bulk RNAseq
  esBulk = es4Bulk.construct.full$bulk.counts
  
  cellTypes.Ref = unique(esRef$cellTypeID)
  nCellTypes.Ref = length(unique(cellTypes.Ref))

  # Estimate cell type proportions for bulk RNAseq 
  Est.prop.esBulk = music_prop3(bulk.eset = esBulk, 
                               sc.eset = esRef,
                               clusters = 'cellTypeID', samples = 'subj', 
                               select.ct = cellTypes.Ref, verbose = TRUE)

  # pseudo R2 for MuSiC
  # nSubj.Bulk x 1 vector
  R2.MuSiC = Est.prop.esBulk$r.squared.full

  # R2 for NNLS
  # nSubj.Ref x 1 vector
  R2.NNLS = Est.prop.esBulk$R2.NNLS

  ####################################
  # estimated cell type proportions via MuSiC
  # rows are subjects; columns are cell types
  est.prop.MuSiC= data.matrix(Est.prop.esBulk$Est.prop.weighted)
  # make sure columns are matched
  pos.order.MuSiC=match(colnames(prop.real.Bulk), colnames(est.prop.MuSiC))
  est.prop.MuSiC2=est.prop.MuSiC[, pos.order.MuSiC]

  # estimated cell type deconvolution via NNLS
  # rows are subjects; columns are cell types
  est.prop.NNLS = data.matrix(Est.prop.esBulk$Est.prop.allgene)
  pos.order.NNLS=match(colnames(prop.real.Bulk), colnames(est.prop.NNLS))
  est.prop.NNLS2=est.prop.NNLS[, pos.order.NNLS]

  # list of estimated cell type proportions
  prop.est = list(MuSiC=est.prop.MuSiC2, 
                  NNLS = est.prop.NNLS2)
  

  #######################
  # calculate performance indices and indices after regressing out cell type
  ######################
  pIndex = performIndexFunc(prop.est = prop.est, prop.real.Bulk = prop.real.Bulk)
  performMat = pIndex$performMat
  prop.est.LF = pIndex$prop.est.LF
  prop.est.LF.resi = pIndex$prop.est.LF.resi
  

  ########################################################
  # bootstrapping reference scRNAseq data based on algorithm used by 
  # Menden et al., Sci. Adv. 2020; 6 : eaba2619 22 July 2020
  ########################################################

  # for nBoot bootstrap data

  ##
  R2.MuSiCMat = matrix(0, nrow=ncol(esBulk), ncol=nBoot)
  colnames(R2.MuSiCMat)=paste("boot", 1:nBoot, sep="")

  R2.NNLSMat = matrix(0, nrow=ncol(esBulk), ncol=nBoot)
  colnames(R2.NNLSMat)=paste("boot", 1:nBoot, sep="")

  ###

  est.MuSiC = matrix(0, nrow=ncol(esBulk), ncol=nCellTypes.Ref)
  est.NNLS = matrix(0, nrow=ncol(esBulk), ncol=nCellTypes.Ref)
  
  # set random seed
  set.seed(seed)
  for(m in 1:nBoot)
  {
    cat(" m =", m)
    # get a bootstrapping reference scRNAseq dataset
    esRef.m = genRefscRNAseqData.default(sc.eset = esRef, 
                                       cellType.var = "cellType", 
                                       samples.var = "subj")
  
    esRef.m$cellType=factor(esRef.m$cellType)
    esRef.m$cellTypeID = esRef.m$cellType
  
    # estimate cell type proportions
    Est.prop.esBulk.lst = music_prop3(bulk.eset = esBulk, 
                                          sc.eset = esRef.m,
                                          clusters = 'cellTypeID', samples = 'subj', 
                                          select.ct = cellTypes.Ref, verbose = FALSE)
  
    # pseudo R2
    R2.MuSiCMat[,m] = Est.prop.esBulk.lst$r.squared.full
    # R2
    R2.NNLSMat[,m] = Est.prop.esBulk.lst$R2.NNLS
    
    # Estimation evaluation
    est.prop.MuSiC.lst = data.matrix(Est.prop.esBulk.lst$Est.prop.weighted)
    # make sure columns are matched
    pos.order.MuSiC=match(colnames(prop.real.Bulk), colnames(est.prop.MuSiC.lst))
    est.prop.MuSiC.lst2=est.prop.MuSiC.lst[, pos.order.MuSiC]
    
    est.MuSiC = est.MuSiC + est.prop.MuSiC.lst2

    ###  
    est.prop.NNLS.lst = data.matrix(Est.prop.esBulk.lst$Est.prop.allgene)
    pos.order.NNLS=match(colnames(prop.real.Bulk), colnames(est.prop.NNLS.lst))
    est.prop.NNLS.lst2=est.prop.NNLS.lst[, pos.order.NNLS]

    est.NNLS = est.NNLS + est.prop.NNLS.lst2
  }

  # get mean cell type proportion estimates across bootstrapping
  for(i in 1:nrow(est.MuSiC))
  {
    est.MuSiC[i,] = est.MuSiC[i,]/sum(est.MuSiC[i,], na.rm=TRUE)
    est.NNLS[i,] = est.NNLS[i,]/sum(est.NNLS[i,], na.rm=TRUE)
  }

  # match columns 
  pos.order.MuSiC2=match(colnames(prop.real.Bulk), colnames(est.MuSiC))
  est.MuSiC.boot=est.MuSiC[, pos.order.MuSiC2]

  pos.order.NNLS2=match(colnames(prop.real.Bulk), colnames(est.NNLS))
  est.NNLS.boot=est.NNLS[, pos.order.NNLS2]

  prop.est.boot = list(MuSiC=est.MuSiC.boot, 
                   NNLS = est.NNLS.boot) 

  #######################
  # calculate performance indices and indices after regressing out cell type
  # for bootstrapping results
  ######################
  pIndex.boot = performIndexFunc(prop.est = prop.est.boot, prop.real.Bulk = prop.real.Bulk)
  performMat.boot = pIndex.boot$performMat
  prop.est.boot.LF = pIndex.boot$prop.est.LF
  prop.est.boot.LF.resi = pIndex.boot$prop.est.LF.resi
  
  ####
  # r.squared.full, i.e., pseudo R2
  R2.MuSiCFrame=as.data.frame(R2.MuSiCMat)
  # added R2 based on original reference scRNAseq
  R2.MuSiCFrame$orig=R2.MuSiC

  # long format pseudo R2 used to draw trajectory of R2 across bootstrapping
  R2.MuSiCFrame.boot = data.frame(
    R2.MuSiC = c(as.matrix(R2.MuSiCFrame)), 
    grp = rep(colnames(R2.MuSiCFrame), each=nrow(R2.MuSiCFrame)),
    sid = rep(rownames(prop.real.Bulk), ncol(R2.MuSiCFrame)))

  ######
  R2.NNLSFrame=as.data.frame(R2.NNLSMat)
  R2.NNLSFrame$orig=R2.NNLS

  # long format R2 used to draw trajectory of R2 across bootstrapping
  R2.NNLSFrame.boot = data.frame(
    R2.NNLS = c(as.matrix(R2.NNLSFrame)), 
    grp = rep(colnames(R2.NNLSFrame), each=nrow(R2.NNLSFrame)),
    sid = rep(rownames(prop.real.Bulk), ncol(R2.NNLSFrame)))

  ###############################################
  # update esBulk's phenotype data slot
  ###############################################
  pDat = colData(esBulk)
  
  prop.est.MuSiC = as.data.frame(prop.est.boot$MuSiC)
  colnames(prop.est.MuSiC) = paste(colnames(prop.est.MuSiC), ".MuSiC", sep="")
  prop.est.MuSiC$subj = rownames(prop.est.MuSiC)
  pDat2 = merge(x=pDat, y = prop.est.MuSiC, by.x="subj", by.y="subj", sort = FALSE)
  
  prop.est.NNLS = as.data.frame(prop.est.boot$NNLS)
  colnames(prop.est.NNLS) = paste(colnames(prop.est.NNLS), ".NNLS", sep="")
  prop.est.NNLS$subj = rownames(prop.est.NNLS)
  pDat3 = merge(x=pDat2, y = prop.est.NNLS, by.x="subj", by.y="subj", sort = FALSE)
  
  pos = match(pDat$subj, pDat3$subj)
  pDat4 = pDat3[pos,]
  rownames(pDat4) = rownames(pDat)
  
  colData(esBulk) = pDat4  
  
  ###############################################
  # esBulk: ExpressionSet object including estimated cell type proportions in phenotype data slot
  #         so that we can adjust cell type proportions in down-stream analysis
  
  # prop.real.Ref: nSubj.Ref x nCellType matrix of true cell type proportions for reference scRNAseq

  # prop.real.Bulk: # nSubj.Bulk x nCellType matrix of true cell type proportion for constructed Bulk RNAseq
  
  # prop.est: a list of 2 elements. 
  #   1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on original reference scRNAseq
  #   2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on original reference scRNAseq
  

  # prop.est.LF: a list of 2 elements. Long format of the 2 elements in prop.est.
  #   1st element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on original reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is real cell type proportions;
  #     3rd column is subject id;
  #     4th column is cell type.
  #   2nd element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on original reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is real cell type proportions;
  #     3rd column is subject id;
  #     4th column is cell type.  
  
  # prop.est.LF.resi: same structure as prop.est.LF. 'resi' indicates the estimated and real cell type proportions
  # are residuals after regressing out cell types
  
  
  # performMat: a 2x8 matrix of performance indices. 
  #   rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE, mAD, R, R.rank, RMSE.resi, mAD.resi, R.resi, R.rank.resi
  
  # performMat.resi: a 2x8 matrix of performance indices after regressing out cell type effect. 
  #   rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE, mAD, R, R.rank, RMSE.resi, mAD.resi, R.resi, R.rank.resi

  # prop.est.boot: a list of 2 elements. 
  #   1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on ensembled bootstrappings of reference scRNAseq
  #   2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on ensembled bootstrappings of reference scRNAseq

  # prop.est.boot.LF: a list of 2 elements. Long format of the 2 elements in prop.est.boot.
  #   1st element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on ensembled bootstrappings of reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is real cell type proportions;
  #     3rd column is subject id;
  #     4th column is cell type.
  #   2nd element is a [nSubj.Bulk x nCellType] x 4 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on ensembled bootstrappings of reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is real cell type proportions;
  #     3rd column is subject id;
  #     4th column is cell type.
  
  # prop.est.boot.LF.resi: same structure as prop.est.boot.LF. 'resi' indicates the estimated and real cell type proportions
  # are residuals after regressing out cell types

  # performMat.boot: a 2x6 matrix of performance indices based on ensembled bootstrappings of reference scRNAseq. 
  #   rows are methods (MuSiC and NNLS); columns are evaluation criteria: RMSE, mAD, and R, RMSE.resi, mAD.resi, R.resi

  # R2.MuSiCFrame.boot: long format pseudo R2 via MuSiC used to draw trajectory of R2 across bootstrapping with 3 columns
  #   1st column is R2.MuSiC; 2nd column indicates if results are based on original reference scRNAseq or
  #     based on which bootstrapping reference scRNAseq

  # R2.NNLSFrame.boot: long format pseudo R2 via MuSiC used to draw trajectory of R2 across bootstrapping with 3 columns
  #   1st column is R2.MuSiC; 2nd column indicates if results are based on original reference scRNAseq or
  #     based on which bootstrapping reference scRNAseq

  res = list(
    esBulk = esBulk,
    prop.real.Ref = prop.real.Ref,
    prop.real.Ref.LF = prop.real.Ref.LF, 
    prop.real.Bulk = prop.real.Bulk,
    prop.real.Bulk.LF = prop.real.Bulk.LF, 
    prop.est = prop.est, 
    prop.est.LF = prop.est.LF,
    prop.est.LF.resi = prop.est.LF.resi, 
    performMat = performMat,
    prop.est.boot = prop.est.boot,
    prop.est.boot.LF = prop.est.boot.LF,
    prop.est.boot.LF.resi = prop.est.boot.LF.resi,
    performMat.boot = performMat.boot,
    R2.MuSiCFrame.boot = R2.MuSiCFrame.boot,
    R2.NNLSFrame.boot = R2.NNLSFrame.boot)

  invisible(res)
}



MuSiCEvaluator2scRNAseq = function(es4Bulk, esRef, sid, cellType, nBoot= 20, seed = 1234567, nPerm = NULL)
{
  
  res = MuSiCEvaluator2scRNAseq.default(es4Bulk = es4Bulk, 
                                        esRef = esRef, 
                                        sid = sid, 
                                        cellType = cellType, 
                                        nBoot= nBoot, 
                                        seed = seed)
  
  # calculate sample mean and sd of pseudo R2
  R2meanSDMat = calR2MeanSD(obj.Evaluator = res)
  res$R2meanSDMat = R2meanSDMat
  
  if(!is.null(nPerm) && is.integer(nPerm) && nPerm > 0)
  {
    R2meanSD.frame = MuSiCEvaluator2scRNAseq.perm(es4Bulk = es4Bulk, 
                                                  esRef = esRef, 
                                                  sid = sid, 
                                                  cellType = cellType, 
                                                  nBoot= nBoot, 
                                                  seed = seed, 
                                                  nPerm = nPerm)
  } else {
    R2meanSD.frame = NULL
  }

  res$R2meanSD.frame = R2meanSD.frame
  
  invisible(res)
}

