#" @export
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
        initialize = function(n, r, x=NULL, f=NULL, replace = FALSE) {
            self$n <- n
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
        collect = function(type = "r") {
            P <- npermutations(self$n, self$r, self$x, self$f, self$replace)
            out <- self$getnext(P, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = "r", drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_permutations(
                    self$n, self$r, d, private$state, self$x, self$f, self$replace, type)
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
    if (r == 0) {
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
            "next_replace_permutations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else if (n < r) {
        if (type == "r") {
            out <- integer(0)
            dim(out) <- c(0, r)
        } else if (type == "c") {
            out <- integer(0)
            dim(out) <- c(r, 0)
        } else {
            out <- list()
        }
    } else if (n == r) {
        out <- .Call(
            "next_permutations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(d),
            state,
            x,
            as.integer(f),
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
            as.integer(f),
            type)
    }
    out
}

#" @export
permutations <- function(n, r=n, x=NULL, f=NULL, replace=FALSE, type = "r") {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    P <- npermutations(n, r, x, f, replace)
    next_permutations(n, r, P, NULL, x, f, replace, type)
}


#" @export
ipermutations <- function(n, r=n, x=NULL, f=NULL, replace = FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    Permutations$new(n, r, x, f, replace)
}

#" @export
npermutations <- function(n, r=n, x=NULL, f=NULL, replace=FALSE, bigz=FALSE) {
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    }
    if (bigz) {
        if (replace) {
            out <- gmp::as.bigz(n) ^ r
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            if (n == r) {
                out <- gmp::factorialZ(n)
            } else {
                out <- out <- .Call("npr_bigz", PACKAGE = "arrangements", as.integer(n), as.integer(r))
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
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            if (n == r) {
                out <- factorial(n)
            } else {
                out <- .Call("npr", PACKAGE = "arrangements", as.integer(n), as.integer(r))
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
