#' @details
#' The `Permutations` class can be initialized by using the convenient wrapper `ipermutations` or
#' \preformatted{
#' Permutations$new(n, k, x = NULL, freq = NULL, replace = FALSE)
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
                out <- next_permutations(
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

next_permutations <- function(n, k, d, state, x, freq, replace, type) {
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
    } else if (n == k) {
        out <- .Call(
            "next_permutations",
            PACKAGE = "arrangements",
            n,
            d,
            state,
            x,
            as_uint_array(freq),
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
            as_uint_array(freq),
            type)
    }
    out
}

#' Permutations generator
#'
#' This function generates all the permutations of selecting `k` items from `n` items.
#' The results are in lexicographical order.
#'
#' @template param_pc
#' @template param_type
#' @seealso [ipermutations] for iterating permutations and [npermutations] to calculate number of permutations
#' @examples
#' permutations(3)
#' permutations(x = LETTERS[1:3])
#'
#' # choose 2 from 4
#' permutations(4, 2)
#' permutations(x = LETTERS[1:3], k = 2)
#'
#' # multiset with frequencies c(2, 3)
#' permutations(freq = c(2, 3), k = 3)
#'
#' # with replacement
#' permutations(4, 2, replace = TRUE)
#'
#' # column major
#' permutations(3, type = "c")
#' permutations(4, 2, type = "c")
#'
#' # list output
#' permutations(3, type = "l")
#' permutations(4, 2, type = "l")
#'
#' # zero sized permutations
#' dim(permutations(0))
#' dim(permutations(5, 0))
#' dim(permutations(5, 6))
#' dim(permutations(0, 0))
#' dim(permutations(0, 1))
#'
#' @export
permutations <- function(n, k=n, x = NULL, freq = NULL, replace = FALSE, type = "r") {
    if (!replace && !is.null(freq)) {
        n <- sum(freq)
        is.null(x) || length(freq) == length(x) || stop("length of x and freq should be the same")
    } else if (!is.null(x)) {
        n <- length(x)
    }
    next_permutations(n, k, -1L, NULL, x, freq, replace, type)
}


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
#' iperm$getnext(type = "c", drop = FALSE)
#' # collect remaining permutations
#' iperm$collect()
#'
#' library(foreach)
#' foreach(x = ipermutations(5, 2), .combine=c) %do% {
#'   sum(x)
#' }
#' @export
ipermutations <- function(n, k=n, x = NULL, freq = NULL, replace = FALSE) {
    if (!replace && !is.null(freq)) {
        n <- sum(freq)
        is.null(x) || length(freq) == length(x) || stop("length of x and freq should be the same")
    } else if (!is.null(x)) {
        n <- length(x)
    }
    Permutations$new(n, k, x, freq, replace)
}

#' Number of permutations
#' @template param_pc
#' @param bigz an logical to indicate using [gmp::bigz]
#' @seealso [permutations] for generating all permutations and [ipermutations] for iterating permutations
#' @examples
#' npermutations(7)
#' npermutations(x = LETTERS[1:5])
#' npermutations(5, 2)
#' npermutations(x = LETTERS, k = 5)
#'
#' # integer overflow
#' \dontrun{npermutations(14, 10)}
#' npermutations(14, 10, bigz = TRUE)
#'
#' # number of permutations of `c("a", "b", "b")`
#' # they are `c("a", "b")`, `c("b", "b")` and `c("b", "b")`
#' npermutations(freq = c(1, 2), k = 2)
#'
#' # zero sized partitions
#' npermutations(0)
#' npermutations(5, 0)
#' npermutations(5, 6)
#' npermutations(0, 1)
#' npermutations(0, 0)
#' @export
npermutations <- function(n, k=n, x = NULL, freq = NULL, replace = FALSE, bigz = FALSE) {
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
            out <- gmp::as.bigz(n) ^ k
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            if (n == k) {
                out <- gmp::factorialZ(n)
            } else {
                out <- out <- .Call("nperm_k_bigz", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("nperm_n_bigz", PACKAGE = "arrangements", as_uint_array(freq))
            } else {
                out <- .Call("nperm_f_bigz", PACKAGE = "arrangements", as_uint_array(freq), k)
            }
        }

    } else {
        if (replace) {
            out <- n ^ k
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            if (n == k) {
                out <- factorial(n)
            } else {
                out <- .Call("nperm_k", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("nperm_n", PACKAGE = "arrangements", as_uint_array(freq))
            } else {
                out <- .Call("nperm_f", PACKAGE = "arrangements", as_uint_array(freq), k)
            }
        }
    }
    convertz(out, bigz)
}
