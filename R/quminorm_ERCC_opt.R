#' ERCC-optimization of quminorm shape
#'
#' Optimizing the shape parameter for the quantile normalization on a per-cell
#' basis, by optimizing the correlation between the observed and the known
#' concentrations of the External RNA Controls Consortium (ERCC) standard
#' transcripts. The underlying assumption is that cells with different skewing
#' of their transcriptomes will need different shape parameters, and as the
#' ERCC should be uniformly PCR-skewed in all cells, they can be used to tune
#' this.
#' @importFrom BiocParallel bplapply
#' @importFrom stats cor
#' @param m a matrix or sparse Matrix of non-negative integers (count data).
#' @param ERCC_conc_data A dataframe or matrix with two columns, where column 1
#' should be a character vector containing the names of the ERCC, and column 2
#' should contain the known concentrations in the input.
#' @param shapeVals The values for the shape that should be evaluated.
#' @param quadpts positive integer, number of quadrature points. Increase for
#' greater precision but slower computation.
#'
#' @return a matrix-like object of same type as m but with the nonzero values
#' normalized to match the target quasi-UMI distribution.
#' @export quminorm_ERCC_opt
quminorm_ERCC_opt <- function(m, ERCC_conc_data,
                              shapeVals = seq(1.4, 2.6, 0.1),
                              quadpts = 1000){
    control_pos <- which(rownames(m) %in% ERCC_conc_data[,1])
    resultList <- bplapply(seq_len(ncol(m)), function(x){
        quminorm_ERCC_obj(m[,x],shapeVals = shapeVals,
                          control_pos = control_pos,
                          control_conc = ERCC_conc_data$conc,
                          quadpts=quadpts)
    })

    fullResult <- do.call("rbind", resultList)
}

quminorm_ERCC_obj <- function(v, shapeVals, quadpts,
                              control_pos, control_conc){
    resList <- lapply(seq_along(shapeVals), function(x){
        resVec <- quminorm_poilog(v,shape = shapeVals[x], quadpts=quadpts)
        corRes <- cor(resVec[control_pos], control_conc, method = "pearson")
        list(resVec, corRes)
    })
    maxPos <- which.max(unlist(lapply(resList, "[[", 2)))
    return(resList[[maxPos]][[1]])

}
