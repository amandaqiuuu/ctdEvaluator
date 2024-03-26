# Nov 4, 2023
# Define a function to convert an ExpressionSet to a SingleCellExperiment

library(Biobase)
library(SingleCellExperiment)


convertES2SCE <- function(es) {
  # Check if the input is an ExpressionSet
  if (!inherits(es, "ExpressionSet")) {
    stop("The input must be an ExpressionSet object.")
  }
  
  # Extract the expression matrix
  exprs_mat <- exprs(es)
  
  # Extract and create colData and rowData
  col_data <- pData(es)
  row_data <- fData(es)
  
  # Create the SingleCellExperiment
  es_sce <- SingleCellExperiment(assays = list(counts = exprs_mat), colData = col_data, rowData = row_data)
  
  # Return the SingleCellExperiment object
  return(es_sce)
}
