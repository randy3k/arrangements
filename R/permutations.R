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
        r = NULL,
        x = NULL,
        f = NULL,
        replace = NULL,
        initialize = function(n, r=n, x=NULL, f=NULL, replace = FALSE) {
            self$n <- as.integer(n)
            self$r <- r
            self$x <- x
            self$f <- f
            self$replace <- replace
            self$reset()
        },
        reset = function() {
            private$state <- new.env()
            private$null_pending <- FALSE
        },
        collect = function(type = 'r') {
            P <- try(npermutations(self$n, self$r, self$f, self$replace), silent = TRUE)
            if (inherits(P, "try-error")) stop("too many results, use `ipermutations`")
            out <- self$getnext(P, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = 'r', drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
            } else {
                out <- next_permutations(
                    self$n, self$r, d, private$state, self$x, self$f, self$replace, type)
            }
            if (is.null(out) || length(out) == 0) {
                self$reset()
            } else if (d > 1) {
                if (type == 'r' && nrow(out) < d){
                    private$null_pending <- TRUE
                } else if (type == 'c' && ncol(out) < d){
                    private$null_pending <- TRUE
                } else if (type == 'l' && length(out) < d){
                    private$null_pending <- TRUE
                }
            }
            if (!is.null(out) && drop) {
                if (type == 'l') {
                    out <- out[[1]]
                } else {
                    dim(out) <- NULL
                }
            }
            out
        },
        print = function(...) {
            if (is.null(self$r)) {
                cat("Permutations of", self$n, "items\n")
            } else {
                cat("Permutations of", self$r, "items from", self$n, "items\n")
            }
            invisible(self)
        }
    )
)

next_permutations <- function(n, r, d, state, x, f, replace, type) {
    if (d > 1) {
        if ((type != 'l' && is.null(r) && d * r > .Machine$integer.max) ||
                (type == 'l' && d > .Machine$integer.max)) {
            stop("too many results, use `ipermutations`")
        }
    }
    if (!is.null(f)) {
        f <- as.integer(f)
    }

    if (replace) {
        out <- .Call(
            "next_replace_permutations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else if (n == r) {
        out <- .Call(
            "next_permutations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(d),
            state,
            x,
            f,
            type)
    } else {
        out <- .Call(
            "next_k_permutations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            f,
            type)
    }

    if (!is.null(out)) {
        if (type == 'r') {
            dim(out) <- c(length(out) / r, r)
        } else if (type == 'c') {
            dim(out) <- c(r, length(out) / r)
        }
    }
    out
}

#' @export
permutations <- function(n, r=n, x=NULL, f=NULL, replace=FALSE, type = 'r') {
    if (is.null(f) && !is.null(x)) {
        n <- length(x)
    } else if (!is.null(f)) {
        n <- sum(f)
    }
    P <- try(npermutations(n, r, f, replace), silent = TRUE)
    if (inherits(P, "try-error")) stop("too many results, use `ipermutations`")
    next_permutations(n, r, P, NULL, x, f, replace, type)
}


#' @export
ipermutations <- function(n, r=n, x=NULL, f=NULL, replace = FALSE) {
    Permutations$new(n, r, x, f, replace)
}

#' @export
npermutations <- function(n, r=n, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n) && !is.null(f)) {
        n <- sum(f)
    }
    if (bigz) {
        if (replace) {
            out <- gmp::as.bigz(n) ^ r
        } else if (is.null(f)) {
            if (n == r) {
                out <- gmp::factorialZ(n)
            } else {
                out <- prod(gmp::as.bigz(tail(seq_len(n), r)))
            }
        } else {
            if (n == r) {
                out <- .Call("multichoose_bigz", PACKAGE = "arrangements", as.integer(f))
            } else {
                out <- .Call("nperm_f_bigz", PACKAGE = "arrangements", as.integer(f), as.integer(r))
            }
        }

    } else {
        if (replace) {
            out <- n ^ r
        } else if (is.null(f)) {
            if (n == r) {
                out <- factorial(n)
            } else {
                out <- prod(tail(seq_len(n), r))
            }
        } else {
            if (n == r) {
                out <- .Call("multichoose", PACKAGE = "arrangements", as.integer(f))
            } else {
                out <- .Call("nperm_f", PACKAGE = "arrangements", as.integer(f), as.integer(r))
            }
        }
    }
    convertz(out, bigz)
}
