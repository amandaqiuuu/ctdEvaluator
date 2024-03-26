# created on June 11, 2021
#  (1) MuSiC cell type deconvolution evaluation via pseudo R2

# esBulk - ExpressionSet object of bulk RNAseq with unknown cell type info
# esRef - ExpressionSet object of reference scRNAseq
# sid - character. column variable name for phenotype data in 'esRef' and 'es4Bulk' indicating subject id
# cellType - character. column variable name for phenotype data in 'esRef' indicating cell type
# nBoot - integer. number of bootstrapping
# seed - integer. random seed to generate exactly the same results for each call of 'MuSiCEvaluator'
MuSiCEvaluator.default = function(esBulk, esRef, sid, cellType, nBoot= 20, seed = 1234567)
{

  pDatBulk = colData(esBulk)
  sidBulk = pDatBulk[, c(sid)]
  if(is.null(sidBulk))
  {
    stop(paste("Error! 'esBulk' does not contain phenotype variable", sid, sep=""))
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

  ####################################################################

  # bulk RNAseq
  esBulk$subj = esBulk$SubjectName

  # Estimate cell type proportions for bulk RNAseq 
  cellTypes.Ref = unique(esRef$cellType)
  nCellTypes.Ref = length(unique(cellTypes.Ref))

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
  ####################################
  # estimated cell type proportions via MuSiC
  # rows are subjects; columns are cell types
  est.prop.MuSiC= data.matrix(Est.prop.esBulk$Est.prop.weighted)
  # make sure columns are matched
  pos.order.MuSiC=match(colnames(prop.real.Ref), colnames(est.prop.MuSiC))
  est.prop.MuSiC2=est.prop.MuSiC[, pos.order.MuSiC]

  # estimated cell type deconvolution via NNLS
  # rows are subjects; columns are cell types
  est.prop.NNLS = data.matrix(Est.prop.esBulk$Est.prop.allgene)
  pos.order.NNLS=match(colnames(prop.real.Ref), colnames(est.prop.NNLS))
  est.prop.NNLS2=est.prop.NNLS[, pos.order.NNLS]

  # list of estimated cell type proportions
  prop.est = list(MuSiC=est.prop.MuSiC2, 
                NNLS = est.prop.NNLS2) 
  
  ####
  # get long format estimated cell type proportions
  ####
  prop.est.MuSiC.LF = data.frame(prop.est=as.vector(prop.est$MuSiC), 
                                 sid = rep(rownames(prop.est$MuSiC), 
                                           ncol(prop.est$MuSiC)),
                                 cellType = rep(colnames(prop.est$MuSiC),
                                                each = nrow(prop.est$MuSiC)))
  
  prop.est.NNLS.LF = data.frame(prop.est=as.vector(prop.est$NNLS), 
                                sid = rep(rownames(prop.est$NNLS), 
                                          ncol(prop.est$NNLS)),
                                cellType = rep(colnames(prop.est$NNLS),
                                               each = nrow(prop.est$NNLS)))

  prop.est.LF = list(MuSiC=prop.est.MuSiC.LF, 
                     NNLS = prop.est.NNLS.LF) 
  


  ########################################################
  # for nBoot simulated data
  ##
  R2.MuSiCMat = matrix(0, nrow=ncol(esBulk), ncol=nBoot)
  colnames(R2.MuSiCMat)=paste("boot", 1:nBoot, sep="")

  R2.NNLSMat = matrix(0, nrow=ncol(esBulk), ncol=nBoot)
  colnames(R2.NNLSMat)=paste("boot", 1:nBoot, sep="")

  ###

  est.MuSiC = matrix(0, nrow=ncol(esBulk), ncol=nCellTypes.Ref)
  est.NNLS = matrix(0, nrow=ncol(esBulk), ncol=nCellTypes.Ref)

  set.seed(seed)
  for(m in 1:nBoot)
  {
    cat(" m =", m)
    esRef.m = genRefscRNAseqData.default(sc.eset = esRef, 
                                         cellType.var = "cellType", 
                                         samples.var = "subj")
  
    esRef.m$cellType=factor(esRef.m$cellType)
    esRef.m$cellTypeID = esRef.m$cellType
  
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
    est.MuSiC = est.MuSiC + est.prop.MuSiC.lst
  
    est.prop.NNLS.lst = data.matrix(Est.prop.esBulk.lst$Est.prop.allgene)
    est.NNLS = est.NNLS + est.prop.NNLS.lst


  }

  # get mean cell type proportion estimates across bootstrapping
  for(i in 1:nrow(est.MuSiC))
  {
    est.MuSiC[i,] = est.MuSiC[i,]/sum(est.MuSiC[i,], na.rm = TRUE)
    est.NNLS[i,] = est.NNLS[i,]/sum(est.NNLS[i,], na.rm = TRUE)
  }

  # match columns 
  pos.order.MuSiC2=match(colnames(prop.real.Ref), colnames(est.MuSiC))
  est.MuSiC.boot=est.MuSiC[, pos.order.MuSiC2]

  pos.order.NNLS2=match(colnames(prop.real.Ref), colnames(est.NNLS))
  est.NNLS.boot=est.NNLS[, pos.order.NNLS2]

  prop.est.boot = list(MuSiC=est.MuSiC.boot, 
                   NNLS = est.NNLS.boot) 

  
  ####
  # get long format estimated cell type proportions
  ####
  prop.est.boot.MuSiC.LF = data.frame(prop.est=as.vector(prop.est.boot$MuSiC), 
                                      sid = rep(rownames(prop.est.boot$MuSiC), 
                                                ncol(prop.est.boot$MuSiC)),
                                      cellType = rep(colnames(prop.est.boot$MuSiC),
                                                     each = nrow(prop.est.boot$MuSiC)))
  
  prop.est.boot.NNLS.LF = data.frame(prop.est=as.vector(prop.est.boot$NNLS), 
                                     sid = rep(rownames(prop.est.boot$NNLS), 
                                               ncol(prop.est.boot$NNLS)),
                                     cellType = rep(colnames(prop.est.boot$NNLS),
                                                    each = nrow(prop.est.boot$NNLS)))
  
  
  prop.est.boot.LF = list(MuSiC=prop.est.boot.MuSiC.LF, 
                          NNLS = prop.est.boot.NNLS.LF) 
  
  ####
  # r.squared.full, i.e., pseudo R2
  R2.MuSiCFrame=as.data.frame(R2.MuSiCMat)
  # added R2 based on original reference scRNAseq
  R2.MuSiCFrame$orig=R2.MuSiC

  # long format pseudo R2 used to draw trajectory of R2 across bootstrapping
  R2.MuSiCFrame.boot = data.frame(
    R2.MuSiC = c(as.matrix(R2.MuSiCFrame)), 
    grp = rep(colnames(R2.MuSiCFrame), each=nrow(R2.MuSiCFrame)),
    sid = rep(rownames(prop.est$MuSiC), ncol(R2.MuSiCFrame)))

  ######
  R2.NNLSFrame=as.data.frame(R2.NNLSMat)
  R2.NNLSFrame$orig=R2.NNLS

  # long format R2 used to draw trajectory of R2 across bootstrapping
  R2.NNLSFrame.boot = data.frame(
    R2.NNLS = c(as.matrix(R2.NNLSFrame)), 
    grp = rep(colnames(R2.NNLSFrame), each=nrow(R2.NNLSFrame)),
    sid = rep(rownames(prop.est$NNLS), ncol(R2.NNLSFrame)))

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

  # esBulk: updated ExpressionSet object including estimated cell type proportions in phenotype data slot
  #         so that we can adjust cell type proportions in down-stream analysis
  
  # prop.real.Ref: nSubj.Ref x nCellType matrix of true cell type proportions for reference scRNAseq
  
  # prop.real.Ref.LF: [nSubj.Ref x nCellType] x 3 data frame of long format of prop.real.Ref.
  #  1st column is real cell type proportions for all subjects (subject by subject);
  #  2nd column is subject id;
  #  3rd column is cell type; 
  
  # prop.est: a list of 2 elements. 
  #   1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on original reference scRNAseq
  #   2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on original reference scRNAseq

  # prop.est.LF: a list of 2 elements. Long format of the 2 elements in prop.est.
  #   1st element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on original reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is subject id;
  #     3rd column is cell type.
  #   2nd element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on original reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is subject id;
  #     3rd column is cell type.  

  # prop.est.boot: a list of 2 elements. 
  #   1st element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on ensembled bootstrappings of reference scRNAseq
  #   2nd element is a nSubj.Bulk x nCellType matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on ensembled bootstrappings of reference scRNAseq
  
  # prop.est.boot.LF: a list of 2 elements. Long format of the 2 elements in prop.est.boot.
  #   1st element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via MuSiC
  #     based on ensembled bootstrappings of reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is subject id;
  #     3rd column is cell type.
  #   2nd element is a [nSubj.Bulk x nCellType] x 3 matrix of estimated cell type proportions for bulk RNAseq via NNLS
  #     based on ensembled bootstrappings of reference scRNAseq.
  #     1st column is estimated cell tye proportions;
  #     2nd column is subject id;
  #     3rd column is cell type.

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
    prop.est = prop.est, 
    prop.est.LF = prop.est.LF, 
    prop.est.boot = prop.est.boot,
    prop.est.boot.LF = prop.est.boot.LF, 
    R2.MuSiCFrame.boot = R2.MuSiCFrame.boot,
    R2.NNLSFrame.boot = R2.NNLSFrame.boot)

  invisible(res)
}


MuSiCEvaluator = function(esBulk, esRef, sid, cellType, nBoot= 20, seed = 1234567, nPerm = NULL)
{
  
  res = MuSiCEvaluator.default(esBulk = esBulk, 
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
    R2meanSD.frame = MuSiCEvaluator.perm(esBulk = esBulk, 
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

