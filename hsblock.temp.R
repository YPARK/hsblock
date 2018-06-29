library(Rcpp)
library(RcppEigen)
library(RcppProgress)
library(Matrix)
library(dplyr)

options(stringsAsFactors = FALSE)

## Read football pairs and construct Adj matrix
pairs <- read.table('data/football.pairs')
nodes <- unlist(pairs) %>% unique() %>% sort()
nodes.tab <- data.frame(v = nodes) %>%
    mutate(i = 1:n())
pairs.tab <- pairs %>%
    left_join(nodes.tab %>% rename(V1 = v, ii = i)) %>%
    left_join(nodes.tab %>% rename(V2 = v, jj = i))
n <- length(nodes)
m <- nrow(pairs.tab)

A <- sparseMatrix(i = pairs.tab$ii,
                  j = pairs.tab$jj,
                  x = rep(1.0, m),
                  dims = c(n, n))

A <- A + t(A)
A[A > 1] <- 1

## A <- A * 3

depth <- 5
K <- 2 ** (depth - 1)
Z <- sparseMatrix(i = sample(K, n, replace = TRUE), j = 1:n, x = rep(1, n),
                  dims = c(K, n))

opt <- list(tree.depth = depth,
            vbiter = 1000,
            record.interval = 1,
            inner.iter = 10,
            rate0 = 0.1,
            decay = -.75,
            delay = 5.,
            economy = FALSE)

Sys.setenv("PKG_CXXFLAGS"="-std=c++14 -Wno-unknown-pragmas")
sourceCpp('src/rcpp_hsblock.cc', verbose = TRUE)

temp = rcpp_hsblock(A, Z, opt)


oo <- order(apply(temp$Z, 2, which.max))

image(A[oo, oo])


image(Matrix(temp$Z[, oo]))




