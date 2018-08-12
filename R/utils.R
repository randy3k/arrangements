validate_n_value <- function(n, v, freq, replace) {
    if (!is.null(freq)) {
        all(freq %% 1 == 0 & freq >= 0) || stop("expect integer")
    }
    if (is.null(v)) {
        if (is.null(freq)) {
            (n %% 1 == 0 && n >= 0) || stop("expect integer")
        } else if (replace) {
            n <- length(freq)
            (length(v) == n) || stop("length(v) != length(freq)")
        } else {
            n <- sum(freq)
        }
    } else {
        if (is.null(freq)) {
            n <- length(v)
        } else if (replace) {
            n <- length(freq)
            (length(v) == n) || stop("length(v) != length(freq)")
        } else {
            n <- sum(freq)
            (length(v) == length(freq)) || stop("length(v) != length(freq)")
        }
    }
    n
}
