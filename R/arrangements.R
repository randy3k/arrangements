#' @name arrangements-package
#' @title Fast Generators and Iterators of Permutations, Combinations and Partitions
#' @docType package
#' @useDynLib arrangements
NULL

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
