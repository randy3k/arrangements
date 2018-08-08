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
            private$null_pending <- FALSE
        },
        collect = function(layout = "row") {
            out <- self$getnext(-1L, layout, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, layout = NULL, drop = d == 1L && is.null(layout)) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_combinations(
                    self$n, self$k, d, private$state, self$v, self$freq, self$replace, layout)
                if (layout == "row" || is.null(layout)){
                    if (nrow(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (nrow(out) < d || ncol(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (layout == "column"){
                    if (ncol(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (ncol(out) < d || nrow(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (layout == "list"){
                    if (length(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (length(out) < d) {
                        private$null_pending <- TRUE
                    }
                    if (length(out) > 1 && drop) {
                        out <- unlist(out)
                    }
                }
            }
            out
        },
        print = function(...) {
            cat("Combinations from", self$n, "to", self$k, "\n")
            invisible(self)
        }
    )
)

next_combinations <- function(n, k, d, state, v, freq, replace, layout) {
    if (k == 0) {
        if (layout == "row" || is.null(layout)) {
            if (is.null(v)) {
                out <- integer(0)
            } else {
                out <- new(typeof(v))
            }
            dim(out) <- c(1, 0)
        } else if (layout == "column") {
            if (is.null(v)) {
                out <- integer(0)
            } else {
                out <- new(typeof(v))
            }
            dim(out) <- c(0, 1)
        } else {
            if (n == 0) {
                if (is.null(freq)) {
                    out <- list(integer(0))
                } else {
                    out <- list(new(typeof(v)))
                }
            } else {
                out <- list()
            }
        }
    } else if (replace) {
        out <- .Call("next_replacement_combinations", PACKAGE = "arrangements",
                        n, k, d, state, v, layout)
    } else if (n < k) {
        if (layout == "row" || is.null(layout)) {
            if (is.null(v)) {
                out <- integer(0)
            } else {
                out <- new(typeof(v))
            }
            dim(out) <- c(0, k)
        } else if (layout == "column") {
            if (is.null(v)) {
                out <- integer(0)
            } else {
                out <- new(typeof(v))
            }
            dim(out) <- c(k, 0)
        } else {
            out <- list()
        }
    } else if (!is.null(freq)) {
        out <- .Call("next_multiset_combinations", PACKAGE = "arrangements",
                        n, k, d, state, v, freq, layout)
    } else {
        out <- .Call("next_combinations", PACKAGE = "arrangements",
                        n, k, d, state, v, layout)
    }
    out
}

#' @title Combinations generator
#'
#' @description
#' This function generates all the combinations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @template param_pc
#' @template param_type
#' @param index a vector of indices of the desired combinations
#' @param nsample sampling random combinations
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
#' # zero sized combinations
#' dim(combinations(5, 0))
#' dim(combinations(5, 6))
#' dim(combinations(0, 0))
#' dim(combinations(0, 1))
#'
#' @export
combinations <- function(
        x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE, layout = "row",
        index = NULL, nsample = NULL) {
    if (missing(x)) {
        n <- validate_n_value(n, v, freq)
    } else {
        if (length(x) == 1 && is.numeric(x)) {
            n <- validate_n_value(x, v, freq)
        } else {
            v <- x
            n <- validate_n_value(n, v, freq)
        }
    }
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    if (is.null(index) && is.null(nsample)) {
        next_combinations(n, k, -1L, NULL, v, freq, replace, layout)
    } else {
        if (gmp::is.bigz(index)) {
            index <- as.character(index)
        } else if (is.numeric(index)) {
            index <- index
        }
        if (replace) {
            .Call("get_replacement_combination", PACKAGE = "arrangements",
                    n, k, v, layout, index, nsample)
        } else if (!is.null(freq)) {
            .Call("get_multiset_combination", PACKAGE = "arrangements",
                    freq, k, v, layout, index, nsample)
        } else {
            .Call("get_combinations", PACKAGE = "arrangements", n, k, v, layout, index, nsample)
        }
    }
}

#' @title Combinations iterator
#' @description
#' This function returns a [Combinations](Combinations-class.html) iterator for iterating
#' combinations of `k` items from `n` items. The iterator allows users to fetch the next
#' combination(s) via the `getnext()` method.
#' @template param_pc
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
icombinations <- function(x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE) {
    if (missing(x)) {
        n <- validate_n_value(n, v, freq)
    } else {
        if (length(x) == 1 && is.numeric(x)) {
            n <- validate_n_value(x, v, freq)
        } else {
            v <- x
            n <- validate_n_value(n, v, freq)
        }
    }
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    Combinations$new(n, k, v, freq, replace)
}

#' Number of combinations
#' @template param_pc
#' @param bigz an logical to indicate using [gmp::bigz]
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
ncombinations <- function(x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE, bigz = FALSE) {
    if (missing(x)) {
        n <- validate_n_value(n, v, freq)
    } else {
        if (length(x) == 1 && is.numeric(x)) {
            n <- validate_n_value(x, v, freq)
        } else {
            v <- x
            n <- validate_n_value(n, v, freq)
        }
    }
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            out <- gmp::chooseZ(n, k)
        } else {
            out <- .Call("num_multiset_combinations_bigz", PACKAGE = "arrangements",
                freq, k)
        }

    } else {
        if (replace) {
            out <- choose(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            out <- choose(n, k)
        } else {
            out <- .Call("num_multiset_combinations", PACKAGE = "arrangements",
                freq, k)
        }
    }
    convertz(out, bigz)
}
