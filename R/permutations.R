#' Number of permutations
#' @template param_pc
#' @param bigz an logical to indicate using [gmp::bigz]
#' @seealso [permutations] for generating all permutations and [ipermutations] for iterating permutations
#' @examples
#' npermutations(7)
#' npermutations(LETTERS[1:5])
#' npermutations(5, 2)
#' npermutations(LETTERS, k = 5)
#'
#' # integer overflow
#' \dontrun{npermutations(14, 10)}
#' npermutations(14, 10, bigz = TRUE)
#'
#' # number of permutations of `c("a", "b", "b")`
#' # they are `c("a", "b")`, `c("b", "b")` and `c("b", "b")`
#' npermutations(k = 2, freq = c(1, 2))
#'
#' # zero sized partitions
#' npermutations(0)
#' npermutations(5, 0)
#' npermutations(5, 6)
#' npermutations(0, 1)
#' npermutations(0, 0)
#' @export
npermutations <- function(x = NULL, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE,
                          bigz = FALSE) {
    .Call("npermutations", PACKAGE = "arrangements", x, k, n, v, freq, replace, bigz)
}


#' Permutations generator
#'
#' This function generates all the permutations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @template param_pc
#' @template param_type
#' @param nitem number of permutations required, usually used with \code{skip}
#' @param skip the number of permutations skipped
#' @param index a vector of indices of the desired permutations
#' @param nsample sampling random permutations
#' @param drop vectorize a matrix or unlist a list
#' @seealso [ipermutations] for iterating permutations and [npermutations] to calculate number of permutations
#' @examples
#' permutations(3)
#' permutations(LETTERS[1:3])
#'
#' # choose 2 from 4
#' permutations(4, 2)
#' permutations(LETTERS[1:3], k = 2)
#'
#' # multiset with frequencies c(2, 3)
#' permutations(k = 3, freq = c(2, 3))
#'
#' # with replacement
#' permutations(4, 2, replace = TRUE)
#'
#' # column major
#' permutations(3, layout = "column")
#' permutations(4, 2, layout = "column")
#'
#' # list output
#' permutations(3, layout = "list")
#' permutations(4, 2, layout = "list")
#'
#' # specifc range of permutations
#' permutations(4, 2, nitem = 2, skip = 3)
#'
#' # specific permutations
#' permutations(4, 2, index = c(3, 5))
#'
#' # random permutations
#' permutations(4, 2, nsample = 3)
#'
#' # zero sized permutations
#' dim(permutations(0))
#' dim(permutations(5, 0))
#' dim(permutations(5, 6))
#' dim(permutations(0, 0))
#' dim(permutations(0, 1))
#'
#' @export
permutations <- function(x = NULL, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE,
                         layout = NULL, nitem = -1L, skip = NULL, index = NULL, nsample = NULL, drop = NULL) {
    .Call("get_permutations", PACKAGE = "arrangements",
          x, k, n, v, freq, replace, layout, nitem, index, nsample, NULL, skip, drop)
}


#' @details
#' The `Permutations` class can be initialized by using the convenient wrapper `ipermutations` or
#' \preformatted{
#' Permutations$new(n, k, v = NULL, freq = NULL, replace = FALSE)
#' }
#' @template iterator_methods
#' @rdname ipermutations
#' @export
Permutations <- R6::R6Class(
    "Permutations",
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
        initialize = function(n, k, v = NULL, freq = NULL, replace = FALSE) {
            self$n <- as.integer(n)
            self$k <- as.integer(k)
            self$v <- v
            self$freq <- freq
            self$replace <- replace
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
                out <- .Call("get_permutations", PACKAGE = "arrangements",
                             NULL, self$k, self$n, self$v, self$freq, self$replace, layout,
                             d, NULL, NULL, private$state, NULL, drop)
                is.null(out) && self$reset()
            }
            out
        },
        print = function(...) {
            if (is.null(self$k)) {
                cat("Permutations of", self$n, "items\n")
            } else {
                cat("Permutations of", self$k, "items from", self$n, "items\n")
            }
            invisible(self)
        }
    )
)


#' @title Permutations iterator
#' @description
#' This function returns a [Permutations](Permutations-class.html) iterator for iterating
#' permutations of `k` items from `n` items. The iterator allows users to fetch the next
#' permutation(s) via the `getnext()` method.
#'
#' @template param_pc
#' @seealso [permutations] for generating all permutations and [npermutations] to calculate number of permutations
#' @examples
#' iperm <- ipermutations(5, 2)
#' iperm$getnext()
#' iperm$getnext(2)
#' iperm$getnext(layout = "column", drop = FALSE)
#' # collect remaining permutations
#' iperm$collect()
#'
#' library(foreach)
#' foreach(x = ipermutations(5, 2), .combine=c) %do% {
#'   sum(x)
#' }
#' @export
ipermutations <- function(x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE) {
    if (missing(x)) {
        n <- validate_n_value(n, v, freq, replace)
    } else {
        if (length(x) == 1 && is.numeric(x)) {
            n <- validate_n_value(x, v, freq)
        } else {
            v <- x
            n <- validate_n_value(n, v, freq, replace)
        }
    }
    (k %% 1 == 0 && k >= 0) || stop("expect integer")

    Permutations$new(n, k, v, freq, replace)
}
