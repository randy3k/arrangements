#' Number of partitions
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param bigz an logical to indicate using [gmp::bigz]
#' @seealso [partitions] for generating all partitions and [ipartitions] for iterating partitions
#' @examples
#' # number of partitions of 10
#' npartitions(10)
#' # number of partitions of 10 into 5 parts
#' npartitions(10, 5)
#'
#' # integer overflow
#' \dontrun{npartitions(160)}
#' npartitions(160, bigz = TRUE)
#'
#' # zero sized partitions
#' npartitions(0)
#' npartitions(5, 0)
#' npartitions(5, 6)
#' npartitions(0, 0)
#' npartitions(0, 1)
#' @export
npartitions <- function(n, k = NULL, bigz = FALSE) {
    .Call("npartitions", PACKAGE = "arrangements", n, k, bigz)
}


#' Partitions generator
#'
#' This function partitions an non-negative interger `n` into `k` parts or any part size.
#' The results are in lexicographical or reversed lexicographical order.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending logical to use reversed lexicographical order
#' @template param_type
#' @seealso [ipartitions] for iterating partitions and [npartitions] to calculate number of partitions
#' @examples
#' # all partitions of 6
#' partitions(6)
#' # reversed lexicographical order
#' partitions(6, descending = TRUE)
#'
#' # fixed number of parts
#' partitions(10, 5)
#' # reversed lexicographical order
#' partitions(10, 5, descending = TRUE)
#'
#' # column major
#' partitions(6, layout = "column")
#' partitions(6, 3, layout = "column")
#'
#' # list output
#' partitions(6, layout = "list")
#' partitions(6, 3, layout = "list")
#'
#' # zero sized partitions
#' dim(partitions(0))
#' dim(partitions(5, 0))
#' dim(partitions(5, 6))
#' dim(partitions(0, 0))
#' dim(partitions(0, 1))
#'
#' @export
partitions <- function(n, k = NULL, descending = FALSE, layout = NULL) {
    .Call("get_partitions", PACKAGE = "arrangements",
          n, k, descending, layout, -1, NULL, NULL, NULL, NULL, FALSE)
}


#' @details
#' The `Partitions` class can be initialized by using the convenient wrapper `ipartitions` or
#' \preformatted{
#' Partitions$new(n, k = NULL, descending = FALSE)
#' }
#' @template iterator_methods
#' @rdname ipartitions
#' @export
Partitions <- R6::R6Class(
    "Partitions",
    inherit = Arrangements,
    private = list(
        state = NULL,
        null_pending = FALSE
    ),
    public = list(
        n = NULL,
        k = NULL,
        descending = NULL,
        initialize = function(n, k = NULL, descending = FALSE) {
            self$n <- as.integer(n)
            if (!is.null(k)) {
                self$k <- as.integer(k)
            }
            self$descending <- descending
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
            if (private$state$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- .Call("get_partitions", PACKAGE = "arrangements",
                             self$n, self$k, self$descending, layout, d, NULL, NULL,
                             private$state, NULL, drop)
                is.null(out) && self$reset()
            }
            out
        },
        print = function(...) {
            if (is.null(self$k)) {
                cat("Partitions of", self$n, "\n")
            } else {
                cat("Partitions of", self$n, "into", self$k, "parts\n")
            }
            invisible(self)
        }
    )
)


#' @title Partitions iterator
#' @description
#' This function returns a [Partitions](Partitions-class.html) iterator for iterating
#' partitions of an non-negative integer `n` into `k` parts or any part size.
#' The iterator allows users to fetch the next partition(s) via the `getnext()` method.
#'
#' @param n an non-negative integer to be partitioned
#' @param k number of parts
#' @param descending logical to use reversed lexicographical order
#' @seealso [partitions] for generating all partitions and [npartitions] to calculate number of partitions
#' @examples
#' ipart <- ipartitions(10)
#' ipart$getnext()
#' ipart$getnext(2)
#' ipart$getnext(layout = "column", drop = FALSE)
#' # collect remaining partitions
#' ipart$collect()
#'
#' library(foreach)
#' foreach(x = ipartitions(6, 2), .combine=c) %do% {
#'   prod(x)
#' }
#' @export
ipartitions <- function(n, k = NULL, descending = FALSE) {
    (n %% 1 == 0  && n >= 0) || stop("expect integer")
    if (!is.null(k)) {
        (k %% 1 == 0 && k >= 0) || stop("expect integer")
    }
    Partitions$new(n, k, descending)
}
