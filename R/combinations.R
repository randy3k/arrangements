#' Combinations class
#'
#' An R6 class of combinations iterator. [icombinations] is a convenient wrapper for initializing the class.
#'
#' @section Usage:
#' \preformatted{
#' Combinations$new(n, k, x = NULL, f = NULL, replace = FALSE)
#' ...$getnext(d = 1L, type = "r", drop = d == 1L)
#' ...$collect(type = "r")
#' }
#' @name Combinations-class
#' @seealso [icombinations]
NULL

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
        f = NULL,
        replace = NULL,
        initialize = function(n, k, x = NULL, f = NULL, replace = FALSE) {
            (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
            (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
            self$n <- as.integer(n)
            self$k <- as.integer(k)
            self$x <- x
            self$f <- as_uint_array(f)
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
        getnext = function(d = 1L, type = "r", drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_combinations(
                    self$n, self$k, d, private$state, self$x, self$f, self$replace, type)
                if (type == "r"){
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
                    if (length(out) > 0 && drop) {
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

next_combinations <- function(n, k, d, state, x, f, replace, type) {
    if (k == 0) {
        if (type == "r") {
            out <- integer(0)
            dim(out) <- c(1, 0)
        } else if (type == "c") {
            out <- integer(0)
            dim(out) <- c(0, 1)
        } else {
            out <- list()
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
        if (type == "r") {
            out <- integer(0)
            dim(out) <- c(0, k)
        } else if (type == "c") {
            out <- integer(0)
            dim(out) <- c(k, 0)
        } else {
            out <- list()
        }
    } else if (is.null(f)) {
        out <- .Call(
            "next_combinations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            type)
    } else {
        out <- .Call(
            "next_multiset_combinations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            as_uint_array(f),
            type)
    }
    out
}

#' Combinations generator
#'
#' This function generates all the combinations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @inheritParams ncombinations
#' @param type if "r", "c" or "l" is specified, the return results would be a
#'  "row" matrix, "column" matrix or a list respectively
#' @return a matrix if `type` is "r" or "c", a list if `type` is "l".
#' @seealso [icombinations] for iterating combinations and [ncombinations] to calculate number of combinations
#' @export
combinations <- function(n, k, x=NULL, f=NULL, replace=FALSE, type = "r") {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    next_combinations(n, k, -1L, NULL, x, f, replace, type)
}

#' Combinations iterator
#' @inheritParams ncombinations
#' @seealso [combinations] for generating combinations and [ncombinations] to calculate number of combinations
#' @export
icombinations <- function(n, k, x=NULL, f=NULL, replace = FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    Combinations$new(n, k, x, f, replace)
}

#' Number of combinations
#' @param n an integer, would be determined implicitly from `x` or `f` if missing
#' @param k an integer
#' @param x an optional vector indicating item labels
#' @param f an integer vector of item repeat frequencies
#' @param replace an logical to select with replacement
#' @seealso [combinations] for generating combinations and [icombinations] for iterating combinations
#' @export
ncombinations <- function(n, k, x=NULL, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
    (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
    if (!is.null(f)) {
        f <- as_uint_array(f)
    }
    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(f)) {
            out <- gmp::chooseZ(n, k)
        } else {
            out <- .Call("ncomb_f_bigz", PACKAGE = "arrangements", as_uint_array(f), k)
        }

    } else {
        if (replace) {
            out <- choose(n + k - 1, k)
        } else if (n < k) {
            out <- 0
        } else if (is.null(f)) {
            out <- choose(n, k)
        } else {
            out <- .Call("ncomb_f", PACKAGE = "arrangements", as_uint_array(f), k)
        }
    }
    convertz(out, bigz)
}
