#' @name arrangements-package
#' @docType package
#' @importFrom methods new
#' @importFrom R6 R6Class
#' @useDynLib arrangements
"_PACKAGE"

Arrangements <- R6::R6Class(
    c("Arrangements", "iter", "abstractiter"),
    private = list(
        state = NULL
    ),
    public = list(
        nextElem = function() {
            out <- self$getnext()
            is.null(out) && stop("StopIteration", call. = FALSE)
            out
        }
    )
)
