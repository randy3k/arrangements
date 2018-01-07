convertz <- function(z, bigz){
    if (bigz){
        return(gmp::as.bigz(z))
    } else if (abs(z) < .Machine$integer.max) {
        return(as.integer(round(z)))
    } else {
        stop("integer overflow, consider using big integer")
    }
}

validate_n <- function(n, x, f){
    if (!is.null(f)) {
        if (!all((f >= 0) & (f %% 1 == 0))) {
            stop("f must be non-negative")
        }
    }
    if (missing(n)) {
        if (is.null(f) && !is.null(x)) {
            n <- length(x)
        } else if (!is.null(f)) {
            n <- sum(f)
        }
    } else {
        if (is.null(f) && !is.null(x) && n != length(x)) {
            stop("n does not equal to length(x)")
        } else if (!is.null(f)) {
            if (n != sum(f)) stop("n does not equal to sum(f)")
        }
    }
    (n >= 0 && n %% 1 == 0) || stop("n should be a non-negative integer")
    n
}

validate_r <- function(r, n, replace=TRUE) {
    if (missing(r)) {
        r <- n
    }
    (r >= 0 && r %% 1 == 0) || stop("r should be a non-negative integer")
    r
}
