#' Number of compositions
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param bigz an logical to use [gmp::bigz]
#' @seealso [compositions] for generating all compositions and [icompositions] for iterating compositions
#' @examples
#' # number of compositions of 10
#' ncompositions(10)
#' # number of compositions of 10 into 5 parts
#' ncompositions(10, 5)
#'
#' # integer overflow
#' \dontrun{ncompositions(160)}
#' ncompositions(160, bigz = TRUE)
#'
#' # zero sized compositions
#' ncompositions(0)
#' ncompositions(5, 0)
#' ncompositions(5, 6)
#' ncompositions(0, 0)
#' ncompositions(0, 1)
#' @export
ncompositions <- function(n, k = NULL, bigz = FALSE) {
    .Call(C_ncompositions, n, k, bigz)
}


#' Compositions generator
#'
#' This function generates the compositions of an non-negative interger `n` into `k` parts or parts of any sizes.
#' The results are in lexicographical or reversed lexicographical order.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending an logical to use reversed lexicographical order
#' @template param_type
#' @param nitem number of compositions required, usually used with \code{skip}
#' @param skip the number of compositions skipped
#' @param index a vector of indices of the desired compositions
#' @param nsample sampling random compositions
#' @param drop vectorize a matrix or unlist a list
#' @seealso [icompositions] for iterating compositions and [ncompositions] to calculate number of compositions
#' @examples
#' # all compositions of 4
#' compositions(4)
#' # reversed lexicographical order
#' compositions(4, descending = TRUE)
#'
#' # fixed number of parts
#' compositions(6, 3)
#' # reversed lexicographical order
#' compositions(6, 3, descending = TRUE)
#'
#' # column major
#' compositions(4, layout = "column")
#' compositions(6, 3, layout = "column")
#'
#' # list output
#' compositions(4, layout = "list")
#' compositions(6, 3, layout = "list")
#'
#' # zero sized compositions
#' dim(compositions(0))
#' dim(compositions(5, 0))
#' dim(compositions(5, 6))
#' dim(compositions(0, 0))
#' dim(compositions(0, 1))
#'
#' @export
compositions <- function(n, k = NULL, descending = FALSE, layout = NULL,
                       nitem = -1L, skip = NULL, index = NULL, nsample = NULL, drop = NULL) {
    .Call(C_get_compositions,
          n, k, descending, layout, nitem, index, nsample, NULL, skip, drop)
}


#' @details
#' The `Compositions` class can be initialized by using the convenient wrapper `icompositions` or
#' \preformatted{
#' Compositions$new(n, k = NULL, descending = FALSE)
#' }
#' @template iterator_methods
#' @rdname icompositions
#' @export
Compositions <- R6::R6Class(
    "Compositions",
    inherit = Arrangements,
    private = list(
        state = NULL,
        null_pending = FALSE
    ),
    public = list(
        n = NULL,
        k = NULL,
        descending = NULL,
        skip = NULL,
        initialize = function(n, k = NULL, descending = FALSE, skip = NULL) {
            self$n <- as.integer(n)
            if (!is.null(k)) {
                self$k <- as.integer(k)
            }
            self$descending <- descending
            self$skip <- skip
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$state$null_pending <- FALSE
        },
        collect = function(layout = "row") {
            out <- self$getnext(-1L, layout, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, layout = NULL, drop = NULL) {
            if (d > 0 && private$state$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- .Call(C_get_compositions,
                             self$n, self$k, self$descending, layout, d, NULL, NULL,
                             private$state, self$skip, drop)
                is.null(out) && self$reset()
            }
            out
        },
        print = function(...) {
            if (is.null(self$k)) {
                cat("Compositions of", self$n, "\n")
            } else {
                cat("Compositions of", self$n, "into", self$k, "parts\n")
            }
            invisible(self)
        }
    )
)


#' @title Compositions iterator
#' @description
#' This function returns a [Compositions] iterator for iterating
#' compositions of an non-negative integer `n` into `k` parts or parts of any sizes.
#' The iterator allows users to fetch the next partition(s) via the `getnext()` method.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending an logical to use reversed lexicographical order
#' @param skip the number of compositions skipped
#' @seealso [compositions] for generating all compositions and [ncompositions] to calculate number of compositions
#' @examples
#' ipart <- icompositions(4)
#' ipart$getnext()
#' ipart$getnext(2)
#' ipart$getnext(layout = "column", drop = FALSE)
#' # collect remaining compositions
#' ipart$collect()
#'
#' library(foreach)
#' foreach(x = icompositions(6, 2), .combine=c) %do% {
#'   prod(x)
#' }
#' @export
icompositions <- function(n, k = NULL, descending = FALSE, skip = NULL) {
    (n %% 1 == 0  && n >= 0) || stop("expect integer")
    if (!is.null(k)) {
        (k %% 1 == 0 && k >= 0) || stop("expect integer")
    }
    Compositions$new(n, k, descending, skip)
}
