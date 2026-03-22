context("Combinations")

test_that("Combinations - ncombinations", {
    L = c(LETTERS, letters)
    expect_equal(ncombinations(10, 3), 120)
    expect_equal(ncombinations(x = L[1:10], k = 3), 120)
    expect_error(ncombinations(40, 15), "integer overflow")
    expect_error(ncombinations(x = L[1:40], k = 15), "integer overflow")
    expect_equal(ncombinations(40, 15, bigz = TRUE), gmp::as.bigz("40225345056"))
    expect_equal(ncombinations(10, 0), 1)
    expect_equal(ncombinations(10, 11), 0)
    expect_error(ncombinations(10, -1), "expect integer")
    expect_error(ncombinations(10, 1.5), "expect integer")

    expect_equal(ncombinations(0, 0), 1)
    expect_equal(ncombinations(0, 1), 0)
})

test_that("Combinations - combinations", {
    comb <- combinations(10, 3)
    expect_equal(nrow(comb), 120)
    expect_equal(ncol(comb), 3)
    expect_equal(comb[1, ], 1:3)
    expect_equal(comb[120, ], 8:10)

    comb <- combinations(10, 3, layout = "row")
    expect_equal(nrow(comb), 120)
    expect_equal(ncol(comb), 3)
    expect_equal(comb[1, ], 1:3)
    expect_equal(comb[120, ], 8:10)

    comb <- combinations(10, 3, layout = "column")
    expect_equal(ncol(comb), 120)
    expect_equal(nrow(comb), 3)
    expect_equal(comb[, 1], 1:3)
    expect_equal(comb[, 120], 8:10)

    comb <- combinations(10, 3, layout = "list")
    expect_equal(length(comb), 120)
    expect_equal(comb[[1]], 1:3)
    expect_equal(comb[[120]], 8:10)

    comb <- combinations(x = LETTERS[1:10], k = 3)
    expect_equal(nrow(comb), 120)
    expect_equal(ncol(comb), 3)
    expect_equal(comb[1, ], LETTERS[1:3])
    expect_equal(comb[120, ], LETTERS[8:10])

    expect_error(combinations(40, 15), "too many results")
    expect_error(combinations(10, -1), "expect integer")
    expect_error(combinations(10, 1.5), "expect integer")
    expect_equal(dim(combinations(10, 0)), c(1, 0))
    expect_equal(dim(combinations(10, 11)), c(0, 11))
    expect_equal(dim(combinations(0, 0)), c(1, 0))
    expect_equal(dim(combinations(0, 1)), c(0, 1))
})

test_that("Combinations - icombinations", {
    icomb <- icombinations(10, 3)
    comb <- combinations(10, 3)
    expect_equal(icomb$collect(), comb)
    expect_equal(icomb$getnext(), 1:3)
    expect_equal(icomb$getnext(), c(1, 2, 4))
    icomb$getnext(110)
    expect_equal(nrow(icomb$getnext(10)), 8)
    expect_equal(icomb$getnext(), NULL)

    comb <- combinations(10, 3, layout = "row")
    expect_equal(icomb$collect(layout = "row"), comb)
    expect_equal(icomb$getnext(layout = "row"), t(1:3))
    expect_equal(icomb$getnext(layout = "row"), t(c(1, 2, 4)))
    icomb$getnext(110, layout = "row")
    expect_equal(nrow(icomb$getnext(10, layout = "row")), 8)
    expect_equal(icomb$getnext(layout = "row"), NULL)

    comb <- combinations(10, 3, layout = "column")
    expect_equal(icomb$collect(layout = "column"), comb)
    expect_equal(icomb$getnext(layout = "column"), t(t(1:3)))
    expect_equal(icomb$getnext(layout = "column"), t(t(c(1, 2, 4))))
    icomb$getnext(110, layout = "column")
    expect_equal(ncol(icomb$getnext(10, layout = "column")), 8)
    expect_equal(icomb$getnext(layout = "column"), NULL)

    comb <- combinations(10, 3, layout = "list")
    expect_equal(icomb$collect(layout = "list"), comb)
    expect_equal(icomb$getnext(layout = "list"), list(1:3))
    expect_equal(icomb$getnext(layout = "list"), list(c(1, 2, 4)))
    icomb$getnext(110, layout = "list")
    expect_equal(length(icomb$getnext(10, layout = "list")), 8)
    expect_equal(icomb$getnext(layout = "list"), NULL)

    icomb <- icombinations(10, 0)
    expect_equal(dim(icomb$collect()), c(1, 0))
    expect_equal(icomb$getnext(), integer(0))
    icomb <- icombinations(10, 11)
    expect_equal(dim(icomb$collect()), c(0, 11))
    expect_equal(icomb$getnext(), NULL)
    expect_error(icombinations(10, -1), "expect integer")
    expect_error(icombinations(10, 1.5), "expect integer")
})

test_that("Combinations - index", {
    expect_equal(combinations(5, 3, index = 1:10), combinations(5, 3))
    expect_equal(combinations(5, 3, index = as.numeric(1:10)), combinations(5, 3))
    expect_equal(combinations(5, 3, index = as.character(1:10)), combinations(5, 3))
    expect_equal(combinations(5, 3, index = gmp::as.bigz(1:10)), combinations(5, 3))
    expect_equal(combinations(5, 3, index = 2), c(1, 2, 4))
    expect_equal(combinations(5, 3, index = 10), c(3, 4, 5))
    expect_error(combinations(5, 3, index = -1), "invalid index")
    expect_error(combinations(5, 3, index = 1.5), "invalid index")
    expect_error(combinations(5, 3, index = "-1"), "invalid index")
    expect_error(combinations(5, 3, index = "1.5"), "invalid index")

    expect_equal(combinations(50, 30, index = 2), c(1:29, 31))
    expect_equal(combinations(50, 30, index = ncombinations(50, 30, bigz = TRUE)), 21:50)
    expect_error(combinations(50, 30, index = -1), "invalid index")
    expect_error(combinations(50, 30, index = 1.5), "invalid index")
    expect_error(combinations(50, 30, index = "-1"), "invalid index")
    expect_error(combinations(50, 30, index = "1.5"), "invalid index")

    expect_equal(combinations(5, 3, index = 1:10, layout = "row"), combinations(5, 3, layout = "row"))
    expect_equal(combinations(5, 3, index = as.numeric(1:10), layout = "row"), combinations(5, 3, layout = "row"))
    expect_equal(combinations(5, 3, index = as.character(1:10), layout = "row"), combinations(5, 3, layout = "row"))
    expect_equal(combinations(5, 3, index = gmp::as.bigz(1:10), layout = "row"), combinations(5, 3, layout = "row"))
    expect_equal(combinations(5, 3, index = 2, layout = "row")[1, ], c(1, 2, 4))
    expect_equal(combinations(5, 3, index = 10, layout = "row")[1, ], c(3, 4, 5))
    expect_equal(combinations(50, 30, index = 2, layout = "row")[1, ], c(1:29, 31))

    expect_equal(combinations(5, 3, index = 1:10, layout = "column"), combinations(5, 3, layout = "column"))
    expect_equal(combinations(5, 3, index = as.numeric(1:10), layout = "column"), combinations(5, 3, layout = "column"))
    expect_equal(combinations(5, 3, index = as.character(1:10), layout = "column"), combinations(5, 3, layout = "column"))
    expect_equal(combinations(5, 3, index = gmp::as.bigz(1:10), layout = "column"), combinations(5, 3, layout = "column"))
    expect_equal(combinations(5, 3, index = 2, layout = "column")[, 1], c(1, 2, 4))
    expect_equal(combinations(5, 3, index = 10, layout = "column")[, 1], c(3, 4, 5))
    expect_equal(combinations(50, 30, index = 2, layout = "column")[, 1], c(1:29, 31))

    expect_equal(combinations(5, 3, index = 1:10, layout = "list"), combinations(5, 3, layout = "list"))
    expect_equal(combinations(5, 3, index = as.numeric(1:10), layout = "list"), combinations(5, 3, layout = "list"))
    expect_equal(combinations(5, 3, index = as.character(1:10), layout = "list"), combinations(5, 3, layout = "list"))
    expect_equal(combinations(5, 3, index = gmp::as.bigz(1:10), layout = "list"), combinations(5, 3, layout = "list"))
    expect_equal(combinations(5, 3, index = 2, layout = "list"), list(c(1, 2, 4)))
    expect_equal(combinations(5, 3, index = 10, layout = "list"), list(c(3, 4, 5)))
    expect_equal(combinations(50, 30, index = 2, layout = "list"), list(c(1:29, 31)))

    expect_error(combinations(0, 1, index = 1), "invalid index")
    expect_equal(combinations(0, 0, index = 1), integer(0))
    expect_error(combinations(0, 1, index = gmp::as.bigz(1)), "invalid index")
    expect_equal(combinations(0, 0, index = gmp::as.bigz(1)), integer(0))
})

test_that("Combinations - skip", {
    expect_equal(combinations(5, 3, skip = 10), combinations(5, 3))
    expect_equal(combinations(5, 3, skip = 3), combinations(5, 3)[4:10, ])
    expect_equal(combinations(5, 3, skip = 3, nitem = 4), combinations(5, 3)[4:7, ])
    expect_equal(combinations(5, 3, skip = gmp::as.bigz(3), nitem = 4), combinations(5, 3)[4:7, ])
})


test_that("Combinations - small cases", {
    expect_equal(combinations(0), matrix(0, nr = 1, nc = 0))
    expect_equal(combinations(0, 1), matrix(0, nr = 0, nc = 1))
    expect_equal(combinations(1), matrix(1, nr = 1, nc = 1))
    expect_equal(combinations(1, 0), matrix(1, nr = 1, nc = 0))

    icomb <- icombinations(0, 0)
    expect_equal(icomb$getnext(), integer(0))
    expect_equal(icomb$getnext(), NULL)

    icomb <- icombinations(0, 1)
    expect_equal(icomb$getnext(), NULL)
    expect_equal(icomb$getnext(), NULL)

    icomb <- icombinations(1, 1)
    expect_equal(icomb$getnext(), 1)
    expect_equal(icomb$getnext(), NULL)

    icomb <- icombinations(1, 0)
    expect_equal(icomb$getnext(), integer(0))
    expect_equal(icomb$getnext(), NULL)
})

test_that("Combinations - nsample", {
    expect_equal(nrow(combinations(10, 3, nsample = 5)), 5)
    expect_equal(ncol(combinations(10, 3, nsample = 5)), 3)

    set.seed(1)
    res1 <- combinations(10, 3, nsample = 5)
    set.seed(1)
    res2 <- combinations(10, 3, nsample = 5)
    expect_equal(res1, res2)

    # nsample with multisets
    expect_equal(nrow(combinations(freq = c(2, 2, 2), k = 3, nsample = 5)), 5)
})

test_that("Combinations - drop", {
    expect_equal(combinations(5, 2, nitem = 1, drop = TRUE), 1:2)
    expect_is(combinations(5, 2, nitem = 1, drop = FALSE), "matrix")
    expect_equal(dim(combinations(5, 2, nitem = 1, drop = FALSE)), c(1, 2))

    expect_equal(combinations(5, 2, index = 1, drop = TRUE), 1:2)
    expect_is(combinations(5, 2, index = 1, drop = FALSE), "matrix")

    set.seed(1)
    s1 <- combinations(5, 2, nsample = 1, drop = TRUE)
    set.seed(1)
    s2 <- combinations(5, 2, nsample = 1, drop = FALSE)
    expect_equal(s1, as.vector(s2))
})

test_that("Combinations - Factors", {
    f <- factor(c("a", "b", "c"), levels = c("a", "b", "c"))
    comb <- combinations(x = f, k = 2)
    expect_true(is.matrix(comb))
    expect_true(is.factor(comb))
    expect_equal(levels(comb), levels(f))
})

test_that("Combinations - Empty and edge cases", {
    expect_error(combinations(NULL, 2), "n is missing")
    expect_equal(nrow(combinations(1:3, k = 4)), 0)
    expect_error(combinations(n = 5, v = 1:4, k = 2), "n != length(v)", fixed = TRUE)

    # NA in x
    expect_equal(combinations(x = c(1, NA), k = 1), matrix(c(1, NA), nc = 1))

    # freq with zeros
    expect_equal(combinations(freq = c(1, 0, 1), k = 2), matrix(c(1, 3), nr = 1))

    # nsample with k=0
    expect_equal(combinations(5, 0, nsample = 1), integer(0))
    expect_equal(dim(combinations(5, 0, nsample = 1, drop = FALSE)), c(1, 0))

    # bigz skip in iterator
    n <- ncombinations(50, 25, bigz = TRUE)
    icomb <- icombinations(50, 25, skip = n - 1)
    expect_equal(icomb$getnext(), 26:50)
})

