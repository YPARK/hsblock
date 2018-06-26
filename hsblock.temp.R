library(Rcpp)
library(RcppEigen)
library(RcppProgress)
library(Matrix)
library(dplyr)

options(stringsAsFactors = FALSE)

## Read football pairs and construct Adj matrix
pairs <- read.table('football.pairs')
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


Sys.setenv("PKG_CXXFLAGS"="-std=c++14 -Wno-unknown-pragmas")
sourceCpp('src/rcpp_hsblock.cc', verbose = TRUE)

depth <- 5
K <- 2 ** (depth - 1)
Z <- matrix(runif(n * K), nrow = K, ncol = n)
Z <- sweep(Z, 2, apply(Z, 2, sum), `/`)

temp = rcpp_hsblock(A, Z, depth)


c.check = temp$Z %*% A



oo <- order(apply(temp$Z, 2, which.max))


image(A[oo, oo])

