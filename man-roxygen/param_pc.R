#' @param x an integer or a vector, will be treated as \code{n} if integer; otherwise, will be treated as \code{v}.
#'          Should not be specified together with \code{n} and \code{v}.
#' @param k an integer, the number of items drawn, defaults to \code{n} if \code{freq} is \code{NULL} else \code{sum(freq)}
#' @param n an integer, the total number of items, its value may be implicitly deduced from \code{length(v)} or \code{length(freq)}
#' @param v a vector to be drawn, defaults to \code{1:n}.
#' @param freq an integer vector of item repeat frequencies
#' @param replace an logical to draw items with replacement
