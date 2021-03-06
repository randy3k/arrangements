context("Permutations")

test_that("Permutations - npermutations", {
    expect_equal(npermutations(5), 120)
    expect_equal(npermutations(x = LETTERS[1:5]), 120)
    expect_error(npermutations(13), "integer overflow")
    expect_error(npermutations(x = LETTERS[1:13]), "integer overflow")
    expect_equal(npermutations(13, bigz = TRUE), gmp::as.bigz("6227020800"))
    expect_equal(npermutations(0), 1)
    expect_error(npermutations(-1), "expect integer")
    expect_error(npermutations(1.5), "expect integer")
})

test_that("Permutations - permutations", {
    perm <- permutations(5)
    expect_equal(nrow(perm), 120)
    expect_equal(ncol(perm), 5)
    expect_equal(perm[1, ], 1:5)
    expect_equal(perm[120, ], 5:1)

    perm <- permutations(5, layout = "row")
    expect_equal(nrow(perm), 120)
    expect_equal(ncol(perm), 5)
    expect_equal(perm[1, ], 1:5)
    expect_equal(perm[120, ], 5:1)

    perm <- permutations(5, layout = "column")
    expect_equal(ncol(perm), 120)
    expect_equal(nrow(perm), 5)
    expect_equal(perm[, 1], 1:5)
    expect_equal(perm[, 120], 5:1)

    perm <- permutations(5, layout = "list")
    expect_equal(length(perm), 120)
    expect_equal(perm[[1]], 1:5)
    expect_equal(perm[[120]], 5:1)

    perm <- permutations(x = LETTERS[1:5])
    expect_equal(nrow(perm), 120)
    expect_equal(ncol(perm), 5)
    expect_equal(perm[1, ], LETTERS[1:5])
    expect_equal(perm[120, ], LETTERS[5:1])

    expect_error(permutations(13), "too many results")
    expect_error(permutations(-1), "expect integer")
    expect_error(permutations(1.5), "expect integer")
    expect_equal(dim(permutations(0)), c(1, 0))
})

test_that("Permutations - ipermutations", {
    iperm <- ipermutations(5)
    perm <- permutations(5)
    expect_equal(iperm$collect(), perm)
    expect_equal(iperm$getnext(), 1:5)
    expect_equal(iperm$getnext(), c(1, 2, 3, 5, 4))
    iperm$getnext(110)
    expect_equal(nrow(iperm$getnext(10)), 8)
    expect_equal(iperm$getnext(), NULL)

    perm <- permutations(5, layout = "row")
    expect_equal(iperm$collect(), perm)
    expect_equal(iperm$getnext(layout = "row"), t(1:5))
    expect_equal(iperm$getnext(layout = "row"), t(c(1, 2, 3, 5, 4)))
    iperm$getnext(110, layout = "row")
    expect_equal(nrow(iperm$getnext(10, layout = "row")), 8)
    expect_equal(iperm$getnext(layout = "row"), NULL)

    perm <- permutations(5, layout = "column")
    expect_equal(iperm$collect(layout = "column"), perm)
    expect_equal(iperm$getnext(layout = "column"), t(t(1:5)))
    expect_equal(iperm$getnext(layout = "column"), t(t(c(1, 2, 3, 5, 4))))
    iperm$getnext(110, layout = "column")
    expect_equal(ncol(iperm$getnext(10, layout = "column")), 8)
    expect_equal(iperm$getnext(layout = "column"), NULL)

    perm <- permutations(5, layout = "list")
    expect_equal(iperm$collect(layout = "list"), perm)
    expect_equal(iperm$getnext(layout = "list"), list(1:5))
    expect_equal(iperm$getnext(layout = "list"), list(c(1, 2, 3, 5, 4)))
    iperm$getnext(110, layout = "list")
    expect_equal(length(iperm$getnext(10, layout = "list")), 8)
    expect_equal(iperm$getnext(layout = "list"), NULL)

    expect_error(ipermutations(-1), "expect integer")
    expect_error(ipermutations(1.5), "expect integer")
})

test_that("Permutations - index", {
    perm <- permutations(5)
    expect_equal(permutations(5, index = 1:120), perm)
    expect_equal(permutations(5, index = as.numeric(1:120)), perm)
    expect_equal(permutations(5, index = as.character(1:120)), perm)
    expect_equal(permutations(5, index = gmp::as.bigz(1:120)), perm)
    expect_equal(permutations(5, index = 2), c(1, 2, 3, 5, 4))
    expect_equal(permutations(5, index = 120), 5:1)
    expect_equal(permutations(13, index = 2), c(1:11, 13, 12))
    expect_equal(permutations(13, index = "6227020800"), 13:1)

    expect_equal(permutations(0, index = 1), integer(0))
})

test_that("Permutations - skip", {
    expect_equal(permutations(5, skip = 120), permutations(5))
    expect_equal(permutations(5, skip = 3), permutations(5)[4:120, ])
    expect_equal(permutations(5, skip = 3, nitem = 4), permutations(5)[4:7, ])
    expect_equal(permutations(5, skip = gmp::as.bigz(3), nitem = 4), permutations(5)[4:7, ])
})


test_that("Permutations - small cases", {
    expect_equal(permutations(0), matrix(0, nr = 1, nc = 0))
    expect_equal(permutations(0, 1), matrix(0, nr = 0, nc = 1))
    expect_equal(permutations(1), matrix(1, nr = 1, nc = 1))
    expect_equal(permutations(1, 0), matrix(1, nr = 1, nc = 0))

    iperm <- ipermutations(0, 0)
    expect_equal(iperm$getnext(), integer(0))
    expect_equal(iperm$getnext(), NULL)

    iperm <- ipermutations(0, 1)
    expect_equal(iperm$getnext(), NULL)
    expect_equal(iperm$getnext(), NULL)

    iperm <- ipermutations(1, 1)
    expect_equal(iperm$getnext(), 1)
    expect_equal(iperm$getnext(), NULL)

    iperm <- ipermutations(1, 0)
    expect_equal(iperm$getnext(), integer(0))
    expect_equal(iperm$getnext(), NULL)
})
