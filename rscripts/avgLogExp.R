#' @title Compute average log expression
#'
#' @description Custom function to compute log average expression level
#'
#' @param object A Seurat Object
#' @param features Features to compute
#' @param cells Cells to average on
#' @param base A positive or complex number: the base with respect to which logarithms are computed. Defaults to ee=exp(1).
#' @param pseudocount.use Pseudocount for computing log. Defaults to 1.
#' @return Mean in log space
#' @export
#'
avgLogExp <- function(object, features = NULL, cells = NULL, base = exp(1), pseudocount.use = 1) {
  if(is.null(cells)){
    cells = SeuratObject::Cells(object)  # Default is to compute on all cells
  }

  x = SeuratObject::GetAssayData(object, slot = "data", assay = "RNA")[features,cells]

  if(is.numeric(x)){

    log_mean_exp = log(x = mean(x = expm1(x = x)) + pseudocount.use, base = base)

  } else if (is(x, "sparseMatrix")){

  log_mean_exp = log(x = Matrix::rowMeans(x = expm1(x = x)) + pseudocount.use, base = base)

  }

  return(log_mean_exp)
}
