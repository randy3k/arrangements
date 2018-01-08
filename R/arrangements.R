#' @name arrangements-package
#' @title Fast Generators and Iterators of Permutations, Combinations and Partitions
#' @description
#' The iterators allow users to generate arrangements in a memory efficient
#' manner and the generated arrangements are in lexicographical (dictionary)
#' order. Permutations and combinations can be drawn with/without replacement and
#' support multisets. It has been demonstrated that 'arrangements' outperforms
#' most of the existing packages of similar kind.
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
