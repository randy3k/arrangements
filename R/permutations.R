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
                out <- next_permutations(
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

next_permutations <- function(n, k, d, state, v, freq, replace, layout) {
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
                if (is.null(v)) {
                    out <- list(integer(0))
                } else {
                    out <- list(new(typeof(v)))
                }
            } else {
                out <- list()
            }
        }
    } else if (replace) {
        out <- .Call("next_replacement_permutations", PACKAGE = "arrangements",
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
    } else if (n == k) {
        # next_permutations can also handle multiset with r=n
        out <- .Call("next_ordinary_permutations", PACKAGE = "arrangements",
                        n, d, state, v, freq, layout)
    } else if (!is.null(freq)) {
        out <- .Call("next_multiset_permutations", PACKAGE = "arrangements",
                        n, k, d, state, v, freq, layout)
    } else {
        out <- .Call("next_k_permutations", PACKAGE = "arrangements",
                        n, k, d, state, v, layout)
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
#' @param index a vector of indices of the desired permutations
#' @param nsample sampling random permutations
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
#' # zero sized permutations
#' dim(permutations(0))
#' dim(permutations(5, 0))
#' dim(permutations(5, 6))
#' dim(permutations(0, 0))
#' dim(permutations(0, 1))
#'
#' @export
permutations <- function(
        x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE, layout = "row",
        index = NULL, nsample = NULL) {
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
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    if (is.null(index) && is.null(nsample)) {
        next_permutations(n, k, -1L, NULL, v, freq, replace, layout)
    } else {
        if (gmp::is.bigz(index)) {
            index <- as.character(index)
        }
        if (replace) {
            .Call("get_replacement_permutations", PACKAGE = "arrangements",
                n, k, v, layout, index, nsample)
        } else if (!is.null(freq)) {
            .Call("get_multiset_permutation", PACKAGE = "arrangements",
                freq, k, v, layout, index, nsample)
        } else if (k == n) {
            .Call("get_ordinary_permutations", PACKAGE = "arrangements", n, v, layout, index, nsample)
        } else {
            .Call("get_k_permutations", PACKAGE = "arrangements", n, k, v, layout, index, nsample)
        }
    }
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
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    Permutations$new(n, k, v, freq, replace)
}

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
npermutations <- function(x, k = n, n = NULL, v = NULL, freq = NULL, replace = FALSE, bigz = FALSE) {
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
    (k %% 1 == 0 && k >= 0) || stop("expect non-negative integer")

    if (bigz) {
        if (replace) {
            out <- gmp::as.bigz(n) ^ k
        } else if (n < k) {
            out <- 0
        } else if (is.null(freq)) {
            if (n == k) {
                out <- gmp::factorialZ(n)
            } else {
                out <- out <- .Call("num_k_permutations_bigz", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("num_multiset_n_permutations_bigz", PACKAGE = "arrangements", freq)
            } else {
                out <- .Call("num_multiset_permutations_bigz", PACKAGE = "arrangements", freq, k)
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
                out <- .Call("num_k_permutations", PACKAGE = "arrangements", n, k)
            }
        } else {
            if (n == k) {
                out <- .Call("num_multiset_n_permutations", PACKAGE = "arrangements", freq)
            } else {
                out <- .Call("num_multiset_permutations", PACKAGE = "arrangements", freq, k)
            }
        }
    }
    convertz(out, bigz)
}
