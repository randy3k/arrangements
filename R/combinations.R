#" @export
Combinations <- R6::R6Class(
    "Combinations",
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
            P <- tryCatch(ncombinations(self$n, self$r, self$x, self$f, self$replace),
                error = function(e) {
                    if (startsWith(e$message, "integer overflow")) {
                        stop("too many results")
                    } else {
                        stop(e)
                    }
            })
            out <- self$getnext(P, type, drop = FALSE)
            self$reset()
            out
        },
        getnext = function(d = 1L, type = "r", drop = d == 1L) {
            if (private$null_pending) {
                out <- NULL
                self$reset()
            } else {
                out <- next_combinations(
                    self$n, self$r, d, private$state, self$x, self$f, self$replace, type)
                if (type == "r"){
                    if (nrow(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (nrow(out) < d) {
                        private$null_pending <- TRUE
                    }
                    if (!is.null(out) && drop) {
                        dim(out) <- NULL
                    }
                } else if (type == "c"){
                    if (ncol(out) == 0) {
                        out <- NULL
                        self$reset()
                    } else if (ncol(out) < d) {
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
            cat("Combinations of", self$r, " items from", self$n, "items\n")
            invisible(self)
        }
    )
)

next_combinations <- function(n, r, d, state, x, f, replace, type) {
    if (d > 1) {
        if ((type != "l" && is.null(r) && d * r > .Machine$integer.max) ||
                (type == "l" && d > .Machine$integer.max)) {
            stop("too many results")
        }
    }

    if (replace) {
        out <- .Call(
            "next_replace_combinations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else if (r == 0 || n == 0 || n < r) {
        if (type == "l") {
            out <- list()
        } else {
            out <- integer(0)
        }
    } else if (is.null(f)) {
        out <- .Call(
            "next_combinations",
            PACKAGE = "arrangements",
            as.integer(n),
            as.integer(r),
            as.integer(d),
            state,
            x,
            type)
    } else {
        if (!is.null(f)) {
            f <- as.integer(f)
        }
        out <- .Call(
            "next_multiset_combinations",
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
        if (type == "r") {
            if (r > 0) {
                dim(out) <- c(length(out) / r, r)
            } else {
                dim(out) <- c(0, 0)
            }
        } else if (type == "c") {
            if (r > 0) {
                dim(out) <- c(r, length(out) / r)
            } else {
                dim(out) <- c(0, 0)
            }
        }
    }
    out
}

#" @export
combinations <- function(n, r, x=NULL, f=NULL, replace=FALSE, type = "r") {
    n <- validate_n(n, x, f)
    r <- validate_r(r, n, replace)
    P <- tryCatch(ncombinations(n, r, x, f, replace),
        error = function(e) {
            if (startsWith(e$message, "integer overflow")) {
                stop("too many results")
            } else {
                stop(e)
            }
    })
    next_combinations(n, r, P, NULL, x, f, replace, type)
}


#" @export
icombinations <- function(n, r, x=NULL, f=NULL, replace = FALSE) {
    n <- validate_n(n, x, f)
    r <- validate_r(r, n, replace)
    Combinations$new(n, r, x, f, replace)
}

#" @export
ncombinations <- function(n, r, x=NULL, f=NULL, replace=FALSE, bigz=FALSE) {
    n <- validate_n(n, x, f)
    r <- validate_r(r, n)
    if (bigz) {
        if (replace) {
            out <- gmp::chooseZ(n + r - 1, r)
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            out <- gmp::chooseZ(n, r)
        } else {
            out <- .Call("ncomb_f_bigz", PACKAGE = "arrangements", as.integer(f), as.integer(r))
        }

    } else {
        if (replace) {
            out <- choose(n + r - 1, r)
        } else if (n < r) {
            out <- 0
        } else if (is.null(f)) {
            out <- choose(n, r)
        } else {
            out <- .Call("ncomb_f", PACKAGE = "arrangements", as.integer(f), as.integer(r))
        }
    }
    convertz(out, bigz)
}
