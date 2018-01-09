#' Permutations class
#'
#' An R6 class of permutation iterator. [ipermutations] is a convenient wrapper for initializing the class.
#'
#' @section Usage:
#' \preformatted{
#' Permutations$new(n, k, x = NULL, f = NULL, replace = FALSE)
#' ...$getnext(d = 1L, type = "r", drop = d == 1L)
#' ...$collect(type = "r")
#' }
#' @name Permutations-class
#' @seealso [ipermutations]
NULL

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
        x = NULL,
        f = NULL,
        replace = NULL,
        initialize = function(n, k, x=NULL, f=NULL, replace = FALSE) {
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
                out <- next_permutations(
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
            if (is.null(self$k)) {
                cat("Permutations of", self$n, "items\n")
            } else {
                cat("Permutations of", self$k, "items from", self$n, "items\n")
            }
            invisible(self)
        }
    )
)

next_permutations <- function(n, k, d, state, x, f, replace, type) {
    if (k == 0) {
        if (type == "r") {
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
                if (is.null(x)) {
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
            "next_replace_permutations",
            PACKAGE = "arrangements",
            n,
            k,
            d,
            state,
            x,
            type)
    } else if (n < k) {
        if (type == "r") {
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
    } else if (n == k) {
        out <- .Call(
            "next_permutations",
            PACKAGE = "arrangements",
            n,
            d,
            state,
            x,
            as_uint_array(f),
            type)
    } else {
        out <- .Call(
            "next_k_permutations",
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

#' Permutations generator
#'
#' This function generates all the permutations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @inheritParams npermutations
#' @param type if "r", "c" or "l" is specified, the return results would be a
#'  "row" matrix, "column" matrix or a list respectively
#' @return a matrix if `type` is "r" or "c", a list if `type` is "l".
#' @export
permutations <- function(n, k=n, x=NULL, f=NULL, replace=FALSE, type = "r") {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    next_permutations(n, k, -1L, NULL, x, f, replace, type)
}


#' Permutations generator
#' @inheritParams npermutations
#' @export
ipermutations <- function(n, k=n, x=NULL, f=NULL, replace = FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    Permutations$new(n, k, x, f, replace)
}

#' Number of permutations
#' @inheritParams ncombinations
#' @export
npermutations <- function(n, k=n, x=NULL, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    (n %% 1 == 0  && n >= 0) || stop("expect non-negative integer")
    (k %% 1 == 0  && k >= 0) || stop("expect non-negative integer")
    if (bigz) {
        if (replace) {
            out <- gmp::as.bigz(n) ^ k
        } else if (n < k) {
            out <- 0
        } else if (is.null(f)) {
            if (n == k) {
                out <- gmp::factorialZ(n)
            } else {
                out <- out <- .Call("nperm_k_bigz", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("nperm_n_bigz", PACKAGE = "arrangements", as_uint_array(f))
            } else {
                out <- .Call("nperm_f_bigz", PACKAGE = "arrangements", as_uint_array(f), k)
            }
        }

    } else {
        if (replace) {
            out <- n ^ k
        } else if (n < k) {
            out <- 0
        } else if (is.null(f)) {
            if (n == k) {
                out <- factorial(n)
            } else {
                out <- .Call("nperm_k", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("nperm_n", PACKAGE = "arrangements", as_uint_array(f))
            } else {
                out <- .Call("nperm_f", PACKAGE = "arrangements", as_uint_array(f), k)
            }
        }
    }
    convertz(out, bigz)
}
