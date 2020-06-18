#' @title Poisson-lognormal maximum likelihood estimation
#' @rdname poilog_mle
#' @description Compute the maximum likelihood estimates for parameters of the
#' Poisson-lognormal distribution for each column of a matrix.
#'
#' @param object A \code{SingleCellExperiment}, \code{SummarizedExperiment},
#'   \code{matrix} or sparse \code{Matrix} of UMI counts (non-negative integers)
#' @param assayName In case object is a SingleCellExperiment or
#'   SummarizedExperiment, the assay containing the UMI counts.
#' @param mc.cores Positive integer indicating the number of cores to use for
#'   parallel processing. See \code{\link[parallel]{mclapply}}.
#'
#' @return a data frame whose rows correspond to the columns of the input data.
#'   The columns contain the MLEs for the parameters (mu and sig),
#'   the log-likelihood, and
#'   the \code{\link[stats]{BIC}} values. NA values indicate numerical
#'   problems with the MLE fit.
#'
#' @details This is essentially a convenience wrapper around
#'   \code{\link[sads]{fitpoilog}} with optional parallelization. WARNING:
#'   only counts from unique molecular identifiers (UMIs) should be used to
#'   fit the Poisson-lognormal, as read counts and other normalized counts
#'   are poorly fit by this distribution.
#'
#' @export
setMethod(f = "poilog_mle",
          signature = signature(object = "Matrix"),
          definition = function(object, mc.cores = 1){
              poilog_mle_sparse(object, mc.cores=mc.cores)
          })

#' @rdname poilog_mle
#' @export
setMethod(f = "poilog_mle",
          signature = signature(object = "matrix"),
          definition = function(object, mc.cores = 1){
              poilog_mle_dense(object, mc.cores=mc.cores)
          })

#' @rdname poilog_mle
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(f = "poilog_mle",
          signature = signature(object = "SummarizedExperiment"),
          definition = function(object, assayName = "counts", mc.cores = 1){
              m <- assay(object, assayName)
              poilog_mle(m, mc.cores=mc.cores)
          })
