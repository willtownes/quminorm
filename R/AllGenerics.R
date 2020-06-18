#' @title Quantile normalization of non-UMI single-cell gene expression
#' @rdname quminorm
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
    name = "quminorm",
    signature = 'object',
    def = function(object, ...) {
        standardGeneric("quminorm")
    }
)

#' @title Poisson-lognormal maximum likelihood estimation
#' @rdname poilog_mle
#' @param ... for the generic, additional arguments to pass to object-specific
#'   methods.
#' @export
setGeneric(
    name = "poilog_mle",
    signature = 'object',
    def = function(object, ...) {
        standardGeneric("poilog_mle")
    }
)
