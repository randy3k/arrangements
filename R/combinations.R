#' Number of combinations
#' @template param_pc
#' @param bigz an logical to use [gmp::bigz]
#' @seealso [combinations] for generating all combinations and [icombinations] for iterating combinations
#' @examples
#' ncombinations(5, 2)
#' ncombinations(LETTERS, k = 5)
#'
#' # integer overflow
#' \dontrun{ncombinations(40, 15)}
#' ncombinations(40, 15, bigz = TRUE)
#'
#' # number of combinations of `c("a", "b", "b")`
#' # they are `c("a", "b")` and `c("b", "b")`
#' ncombinations(k = 2, freq = c(1, 2))
#'
#' # zero sized combinations
#' ncombinations(5, 0)
#' ncombinations(5, 6)
#' ncombinations(0, 1)
#' ncombinations(0, 0)
#'
#' @export
ncombinations <- function(x = NULL, k = NULL, n = NULL, v = NULL, freq = NULL, replace = FALSE,
                          bigz = FALSE) {
    .Call(C_ncombinations, x, k, n, v, freq, replace, bigz)
}


#' @title Combinations generator
#'
#' @description
#' This function generates all the combinations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @template param_pc
#' @template param_type
#' @param nitem number of combinations required, usually used with \code{skip}
#' @param skip the number of combinations skipped
#' @param index a vector of indices of the desired combinations
#' @param nsample sampling random combinations
#' @param drop vectorize a matrix or unlist a list
#' @seealso [icombinations] for iterating combinations and [ncombinations] to calculate number of combinations
#' @examples
#' # choose 2 from 4
#' combinations(4, 2)
#' combinations(LETTERS[1:3], k = 2)
#'
#' # multiset with frequencies c(2, 3)
#' combinations(k = 3, freq = c(2, 3))
#'
#' # with replacement
#' combinations(4, 2, replace = TRUE)
#'
#' # column major
#' combinations(4, 2, layout = "column")
#'
#' # list output
#' combinations(4, 2, layout = "list")
#'
#' # specifc range of combinations
#' combinations(4, 2, nitem = 2, skip = 3)
#'
#' # specific combinations
#' combinations(4, 2, index = c(3, 5))
#'
#' # random combinations
#' combinations(4, 2, nsample = 3)
#'
#' # zero sized combinations
#' dim(combinations(5, 0))
#' dim(combinations(5, 6))
#' dim(combinations(0, 0))
#' dim(combinations(0, 1))
#'
#' @export
combinations <- function(x = NULL, k = NULL, n = NULL, v = NULL, freq = NULL, replace = FALSE,
                         layout = NULL, nitem = -1L, skip = NULL, index = NULL, nsample = NULL, drop = NULL) {
    .Call(C_get_combinations,
          x, k, n, v, freq, replace, layout, nitem, index, nsample, NULL, skip, drop)
}


#' @details
#' The `Combinations` class can be initialized by using the convenient wrapper `icombinations` or
#' \preformatted{
#' Combinations$new(n, k, v = NULL, freq = NULL, replace = FALSE)
#' }
#' @template iterator_methods
#' @rdname icombinations
#' @export
Combinations <- R6::R6Class(
    "Combinations",
    inherit = Arrangements,
    private = list(
        state = NULL,
        null_pending = FALSE
    ),
    public = list(
        n = NULL,
        k = NULL,
        v = NULL,
        freq = NULL,
        replace = NULL,
        skip = NULL,
        initialize = function(n, k, v = NULL, freq = NULL, replace = FALSE, skip = NULL) {
            self$n <- as.integer(n)
            if (is.null(k)) {
                k <- if (is.null(freq)) n else sum(freq)
            }
            self$k <- as.integer(k)
            self$v <- v
            self$freq <- freq
            self$replace <- replace
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
                out <- .Call(C_get_combinations,
                             NULL, self$k, self$n, self$v, self$freq, self$replace, layout,
                             d, NULL, NULL, private$state, self$skip, drop)
                is.null(out) && self$reset()
            }
            out
        },
        print = function(...) {
            cat("Combinations from", self$n, "to", self$k, "\n")
            invisible(self)
        }
    )
)


#' @title Combinations iterator
#' @description
#' This function returns a [Combinations] iterator for iterating
#' combinations of `k` items from `n` items. The iterator allows users to fetch the next
#' combination(s) via the `getnext()` method.
#' @template param_pc
#' @param skip the number of combinations skipped
#' @seealso [combinations] for generating all combinations and [ncombinations] to calculate number of combinations
#' @examples
#' icomb <- icombinations(5, 2)
#' icomb$getnext()
#' icomb$getnext(2)
#' icomb$getnext(layout = "column", drop = FALSE)
#' # collect remaining combinations
#' icomb$collect()
#'
#' library(foreach)
#' foreach(x = icombinations(5, 2), .combine=c) %do% {
#'   sum(x)
#' }
#' @export
icombinations <- function(x = NULL, k = NULL, n = NULL, v = NULL, freq = NULL, replace = FALSE, skip = NULL) {
    n <- validate_n_value(x, k, n, v, freq, replace)
    Combinations$new(n, k, v, freq, replace, skip)
}
