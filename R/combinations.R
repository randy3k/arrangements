#' @details
#' The `Combinations` class can be initialized by using the convenient wrapper `icombinations` or
#' \preformatted{
#' Combinations$new(n, k, x = NULL, freq = NULL, replace = FALSE)
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
        x = NULL,
        freq = NULL,
        replace = NULL,
        initialize = function(n, k, x = NULL, freq = NULL, replace = FALSE) {
            (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
            (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
            self$n <- as.integer(n)
            self$k <- as.integer(k)
            self$x <- x
            self$freq <- as_uint_array(freq)
            self$replace <- replace
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$null_pending <- FALSE
        },
        collect = function(type = "r") {
            out <- self$getnext(-1L, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = NULL, drop = d == 1L && is.null(type)) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_combinations(
                    self$n, self$k, d, private$state, self$x, self$freq, self$replace, type)
                if (type == "r" || is.null(type)){
                    if (nrow(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (nrow(out) < d || ncol(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (type == "c"){
                    if (ncol(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (ncol(out) < d || nrow(out) == 0) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (type == "l"){
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

next_combinations <- function(n, k, d, state, x, freq, replace, type) {
    if (k == 0) {
        if (type == "r" || is.null(type)) {
            if (is.null(x)) {
                out <- integer(0)
            } else {
                out <- new(typeof(x))
            }
            dim(out) <- c(1, 0)
        } else if (type == "c") {
            if (is.null(x)) {
                out <- integer(0)
            } else {
                out <- new(typeof(x))
            }
            dim(out) <- c(0, 1)
        } else {
            if (n == 0) {
                if (is.null(freq)) {
                    out <- list(integer(0))
                } else {
                    out <- list(new(typeof(x)))
                }
            } else {
                out <- list()
            }
        }
    } else if (replace) {
        out <- .Call(
            "next_replace_combinations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            type)
    } else if (n < k) {
        if (type == "r" || is.null(type)) {
            if (is.null(x)) {
                out <- integer(0)
            } else {
                out <- new(typeof(x))
            }
            dim(out) <- c(0, k)
        } else if (type == "c") {
            if (is.null(x)) {
                out <- integer(0)
            } else {
                out <- new(typeof(x))
            }
            dim(out) <- c(k, 0)
        } else {
            out <- list()
        }
    } else if (!is.null(freq)) {
        out <- .Call(
            "next_multiset_combinations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            as_uint_array(freq),
            type)
    } else {
        out <- .Call(
            "next_combinations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            type)
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
#' @seealso [icombinations] for iterating combinations and [ncombinations] to calculate number of combinations
#' @examples
#' # choose 2 from 4
#' combinations(4, 2)
#' combinations(x = LETTERS[1:3], k = 2)
#'
#' # multiset with frequencies c(2, 3)
#' combinations(freq = c(2, 3), k = 3)
#'
#' # with replacement
#' combinations(4, 2, replace = TRUE)
#'
#' # column major
#' combinations(4, 2, type = "c")
#'
#' # list output
#' combinations(4, 2, type = "l")
#'
#' # zero sized combinations
#' dim(combinations(5, 0))
#' dim(combinations(5, 6))
#' dim(combinations(0, 0))
#' dim(combinations(0, 1))
#'
#' @export
combinations <- function(n, k, x = NULL, freq = NULL, replace = FALSE, type = "r") {
    if (!replace && !is.null(freq)) {
        n <- sum(freq)
        is.null(x) || length(freq) == length(x) || stop("length of x and freq should be the same")
    } else if (!is.null(x)) {
        n <- length(x)
    }
    next_combinations(n, k, -1L, NULL, x, freq, replace, type)
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
#' icomb$getnext(type = "c", drop = FALSE)
#' # collect remaining combinations
#' icomb$collect()
#'
#' library(foreach)
#' foreach(x = icombinations(5, 2), .combine=c) %do% {
#'   sum(x)
#' }
#' @export
icombinations <- function(n, k, x = NULL, freq = NULL, replace = FALSE) {
    if (!replace && !is.null(freq)) {
        n <- sum(freq)
        is.null(x) || length(freq) == length(x) || stop("length of x and freq should be the same")
    } else if (!is.null(x)) {
        n <- length(x)
    }
    Combinations$new(n, k, x, freq, replace)
}

#' Number of combinations
#' @template param_pc
#' @param bigz an logical to indicate using [gmp::bigz]
#' @seealso [combinations] for generating all combinations and [icombinations] for iterating combinations
#' @examples
#' ncombinations(5, 2)
#' ncombinations(x = LETTERS, k = 5)
#'
#' # integer overflow
#' \dontrun{ncombinations(40, 15)}
#' ncombinations(40, 15, bigz = TRUE)
#'
#' # number of combinations of `c("a", "b", "b")`
#' # they are `c("a", "b")` and `c("b", "b")`
#' ncombinations(freq = c(1, 2), k = 2)
#'
#' # zero sized combinations
#' ncombinations(5, 0)
#' ncombinations(5, 6)
#' ncombinations(0, 1)
#' ncombinations(0, 0)
#'
#' @export
ncombinations <- function(n, k, x = NULL, freq  =NULL, replace = FALSE, bigz = FALSE) {
    if (!replace && !is.null(freq)) {
        n <- sum(freq)
        is.null(x) || length(freq) == length(x) || stop("length of x and freq should be the same")
    } else if (!is.null(x)) {
        n <- length(x)
    }
    (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
    (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            out <- gmp::chooseZ(n, k)
        } else {
            out <- .Call("ncomb_f_bigz", PACKAGE = "arrangements", as_uint_array(freq), k)
        }

    } else {
        if (replace) {
            out <- choose(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            out <- choose(n, k)
        } else {
            out <- .Call("ncomb_f", PACKAGE = "arrangements", as_uint_array(freq), k)
        }
    }
    convertz(out, bigz)
}
