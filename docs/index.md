# Estimation of Hierarchical Stochastic Block model

Example:
```{r}
## Read pairs
pairs.file <- system.file('extdata', 'football.pairs', package = 'hsblock')
pairs <- read.table(pairs.file, stringsAsFactors = FALSE)

## Construct sparse matrix
net.data <- pairs.to.sparse.matrix(pairs)
A <- net.data$A

## Estimate the block model
out <- fit.hsblock(A, distrib = 'bernoulli', tree.depth = 5, vbiter = 1000, rate = 0.01)

library(Matrix)
oo <- order(apply(out$Z, 2, which.max))
image(A[oo, oo])
```

## Citation

Park and Bader, `https://arxiv.org/abs/1711.05150`

## Release notes

- __0.1.0__ `Rcpp` porting and improved stability
