# created on March 30, 2021
#  (1) move definitions of functions genRefscRNAseqData* from file 'music_prop_v5.R' to this file ('genRefscRNAseqData.R')

# generate data based on algorithm used by 
# Menden et al., Sci. Adv. 2020; 6 : eaba2619 22 July 2020
genRefscRNAseqData.default2 = function(sc.eset, cellType.var, samples.var,
  prop.keep.low=0.90, prop.keep.upp=0.95)
{
  pDat=colData(sc.eset)
  cellTypes = pDat[, c(cellType.var)]
  cellTypes.u= unique(cellTypes)
  # number of cell types
  nCellTypes = length(cellTypes.u)

  subjs=pDat[, c(samples.var)]
  # number of subjects
  subjs.u = unique(subjs)
  nSubjs = length(subjs.u)

  #rcVec = runif(nCellTypes)
  #fcVec = rcVec/sum(rcVec)
  prop.keep=runif(nCellTypes, min=prop.keep.low, max=prop.keep.upp)

  pos.resample = NULL 
  for(j in 1:nSubjs)
  {
    posj = which(subjs == subjs.u[j])
    pDatj = pDat[posj, ,drop=FALSE]
    cellTypes.j=pDatj[, c(cellType.var)]
    # number of cells per cell types
    Nc = table(cellTypes.j)
    # unique cell types
    cellTypes.j.u = names(Nc)
    # new number of cells per cell types
    Nc2 = ceiling(Nc*prop.keep)
    for(k in 1:nCellTypes)
    {
      posjk = which(cellTypes.j == cellTypes.j.u[k]) 
      if(length(posjk))
      {
        # resample
        posjk.r=sample(posjk, size=Nc2[k], replace=FALSE)
        pos.resample = c(pos.resample, posj[posjk.r])
      }
    }
  }
  
  sc.eset.r = sc.eset[, c(pos.resample)]
  sampleNames(sc.eset.r)= paste("s", 1:ncol(sc.eset.r), sep="")

  invisible(sc.eset.r)
}

# generate data based on algorithm used by 
# Menden et al., Sci. Adv. 2020; 6 : eaba2619 22 July 2020
genRefscRNAseqData2 = function(sc.eset, cellType.var, samples.var, nSim=5,
  prop.keep.low=0.90, prop.keep.upp=0.95)
{
  sc.eset.lst=list()
  for(m in 1:nSim)
  {
    sc.eset.lst[[m]] = genRefscRNAseqData.default2(sc.eset = sc.eset, 
			       cellType.var = cellType.var, 
			       samples.var = samples.var,
      prop.keep.low=prop.keep.low, prop.keep.upp=prop.keep.upp
    )
  } 
  names(sc.eset.lst)=paste("sc.eset.lst", 1:nSim, sep="")

  invisible(sc.eset.lst)
}

###################

# generate data based on algorithm used by 
# Menden et al., Sci. Adv. 2020; 6 : eaba2619 22 July 2020
genRefscRNAseqData.default = function(sc.eset, cellType.var, samples.var)
{
  pDat=colData(sc.eset)
  cellTypes = pDat[, c(cellType.var)]
  cellTypes.u= unique(cellTypes)
  # number of cell types
  nCellTypes = length(cellTypes.u)

  subjs=pDat[, c(samples.var)]
  # number of subjects
  subjs.u = unique(subjs)
  nSubjs = length(subjs.u)

  #rcVec = runif(nCellTypes)
  #fcVec = rcVec/sum(rcVec)

  pos.resample = NULL 
  for(j in 1:nSubjs)
  {
    count = 0
    while(count<10)
    {
      count = count + 1
      # generate 'nCellTypes' random number from Uniform[0, 1]
      rcVec = runif(nCellTypes)
      # get a set of cell type proportions for the 'nCellTypes' cell types
      fcVec = rcVec/sum(rcVec)

      # obtain phenotype data for the j-th subject
      posj = which(subjs == subjs.u[j])
      pDatj = pDat[posj, ,drop=FALSE]

      # the original cell types for j-th sbujects
      cellTypes.j=pDatj[, c(cellType.var)]

      # number of cells for subject j
      Ntotal  = length(posj)
      # number of cells per cell types
      Nc = floor(Ntotal*fcVec)
      if(!any(Nc==0))
      {
        break;
      }
    }
    Nc.s = sum(Nc)
    if(Nc.s < Ntotal)
    {
      Nc[nCellTypes] = Nc[nCellTypes] + (Ntotal - Nc.s)
    }

    # for each cell type in subject j, we do resampling with replacement
    for(k in 1:nCellTypes)
    {
      posjk = which(cellTypes.j == cellTypes.u[k]) 
      if(length(posjk))
      {
        # resample with replacement Nc[k] cells
        posjk.r=sample(posjk, size=Nc[k], replace=TRUE)
        pos.resample = c(pos.resample, posj[posjk.r])
      }
    }
  }
  
  sc.eset.r = sc.eset[, c(pos.resample)]
  colnames(sc.eset.r)= paste("s", 1:ncol(sc.eset.r), sep="")

  invisible(sc.eset.r)
}

# generate data based on algorithm used by 
# Menden et al., Sci. Adv. 2020; 6 : eaba2619 22 July 2020
genRefscRNAseqData = function(sc.eset, cellType.var, samples.var, nSim=5)
{
  sc.eset.lst=list()
  for(m in 1:nSim)
  {
    sc.eset.lst[[m]] = genRefscRNAseqData.default(sc.eset = sc.eset, 
			       cellType.var = cellType.var, 
			       samples.var = samples.var)
  } 
  names(sc.eset.lst)=paste("sc.eset.lst", 1:nSim, sep="")

  invisible(sc.eset.lst)
}