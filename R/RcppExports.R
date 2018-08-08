#' Fit a hierarchical stochastic block on the network data
#'
#' @param A [n x n] adjacency matrix
#' @param Z   [K x n] latent assignment matrix#'
#' @param distrib type of distribution for edge weights (bernoulli, poisson)
#' @param degree.corrected fit degree-corrected model or not (default: FALSE)
#'
#' @param Opt a list of options to feed to set the following options simultaneously
#' @param tree.depth depth of binary tree model should be greater than 1 (default: 2)
#' @param rseed random seed for Monte Carlo steps (default: 19937)
#' @param vbiter number of variational Bayes iterations (default: 1000)
#' @param inner.iter number of Monte Carlo sampling (default: 100)
#' @param final.inner.iter number of Monte Carlo sampling at the final step (default: 1000)
#' @param burnin.iter number of Monte Carlo sampling to discard (default: 100)
#' @param record.interval recording interval (default: 10)
#' @param verbose verbosity (default: TRUE)
#' @param rate base rate for stochastic gradient update (default: 0.01)
#' @param decay rate decay parameter (default: -0.55)
#' @param delay rate delay parameter (default: 1)
#'
#'
#' @return a list of variatinal inference results
#'
#' @author Yongjin Park, \email{ypp@@csail.mit.edu}, \email{yongjin.peter.park@@gmail.com}
#' @keywords network clustering EM
#' @seealso \code{pairs.to.sparse.matrix}
#' @seealso \code{build.random.sparse.Z}
#'
#' @examples
#'
#' # Build a network from a list of pairs
#' pairs.file <- system.file('extdata', 'football.pairs', package = 'hsblock')
#' pairs <- read.table(pairs.file, stringsAsFactors = FALSE)
#' net.data <- pairs.to.sparse.matrix(pairs)
#' A <- net.data$A
#'
#' out <- fit.hsblock(A, distrib = 'bernoulli', tree.depth = 5, vbiter = 1000, rate = 0.01)
#'
#' library(Matrix)
#' oo <- order(apply(out$Z, 2, which.max))
#' image(A[oo, oo])
#'
#' y.lim <- range(c(out$llik.variational, out$llik.empirical))
#' plot(out$llik.empirical, ylim = y.lim, type = 'b')
#' points(out$llik.variational, col = 2, type = 'b', pch = 19, cex = .5)
#'
#' @export
fit.hsblock <- function(A,
                        distrib = c('bernoulli', 'poisson'),
                        Opt = list(),
                        Z = NULL,
                        degree.corrected = FALSE,
                        tree.depth = NULL,
                        rseed = NULL,
                        vbiter = NULL,
                        inner.iter = NULL,
                        final.inner.iter = NULL,
                        burnin.iter = NULL,
                        record.interval = NULL,
                        verbose = TRUE,
                        rate = 0.01,
                        decay = -0.55,
                        delay = 1) {

    distrib <- match.arg(distrib)

    n <- nrow(A)
    stopifnot(n == ncol(A))

    default.depth <- 2

    ## update Opt settings if provided
    opt.vars <- c('tree.depth', 'rseed', 'vbiter',
                  'inner.iter', 'final.inner.iter', 'burnin.iter', 'record.interval',
                  'verbose', 'rate', 'decay', 'delay')

    .eval <- function(txt) eval(parse(text = txt))

    for(v in opt.vars) {
        val <- .eval(v)
        if(!(v %in% names(Opt)) && !is.null(val)) {
            Opt[[v]] <- val
        }
    }

    if(!('tree.depth' %in% names(Opt))) {
        Opt$tree.depth <- default.depth
    }

    stopifnot(Opt$tree.depth > 1)

    if(is.null(Z)) {
        K <- 2 ** (Opt$tree.depth - 1)
        Z <- build.random.sparse.Z(K, n)
    }

    if(distrib == 'bernoulli') {
        ret <- .Call('vem_hsb_bern', A, Z, Opt, PACKAGE = 'hsblock')
    } else if (distrib == 'poisson') {
        if(degree.corrected) {
            ret <- .Call('vem_dhsb_pois', A, Z, Opt, PACKAGE = 'hsblock')
        } else {
            ret <- .Call('vem_hsb_pois', A, Z, Opt, PACKAGE = 'hsblock')
        }
    }
    return(ret)
}

#' Random sparse matrix of latent assignment matrix
#' @param n number of vertices
#' @param K number of groups
#' @return [K x n] matrix
#'
#' @examples
#' Z <- build.random.sparse.Z(20, 10)
#' library(Matrix)
#' image(Matrix(Z))
#'
#' @export
build.random.sparse.Z <- function(K, n) {
    ret <- sparseMatrix(i = sample(K, n, replace = TRUE),
                        j = 1:n,
                        x = rep(1, n),
                        dims = c(K, n))
    return(ret)
}

#' Build a sparse matrix from a list pairs or weighted pairs
#'
#' @param pairs a data.frame of network edges (m x 2) or (m x 3) where m is number of edges
#' @param vertices a set of vertices
#' @param default.weight default edge weight (default: 1)
#' @param symmetrize build a symmetric adjacency matrix (default: TRUE)
#'
#' @return a list that contains (A = Matrix::sparseMatrix, V = vertices)
#'
#' @examples
#'
#' pairs.file <- system.file('extdata', 'football.pairs', package = 'hsblock')
#' pairs <- read.table(pairs.file, stringsAsFactors = FALSE)
#' net.data <- pairs.to.sparse.matrix(pairs)
#'
#' library(Matrix)
#' image(net.data$A[1:10, 1:10])
#'
#' @export
pairs.to.sparse.matrix <- function(pairs, vertices = NULL, default.weight = 1, symmetrize = TRUE) {

    stopifnot(ncol(pairs) >= 2)

    pairs <- unique(pairs)

    .df <- function(...) data.frame(..., stringsAsFactors = FALSE)

    if(symmetrize) {
        nc <- ncol(pairs)
        if(nc <= 2) {
            pairs <- rbind(.df(u = pairs[, 1], v = pairs[, 2]),
                           .df(u = pairs[, 2], v = pairs[, 1]))
        } else {
            pairs <- rbind(.df(u = pairs[, 1], v = pairs[, 2], pairs[, -(1:2)]),
                           .df(u = pairs[, 2], v = pairs[, 1], pairs[, -(1:2)]))
        }
    }

    if(is.null(vertices)) {
        vertices <- sort(unique(unlist(pairs[, 1:2])))
    } else {
        valid <- (pairs[, 1] %in% vertices) & (pairs[, 2] %in% vertices)
        pairs <- pairs[valid, , drop = FALSE]
    }

    m <- nrow(pairs)
    n <- length(vertices)

    ii <- match(pairs[, 1], vertices)
    jj <- match(pairs[, 2], vertices)

    if(ncol(pairs) < 3) {
        weights <- rep(default.weight, m)
    } else {
        weights <- pairs[, 3]
    }

    A <- Matrix::sparseMatrix(i = ii,
                              j = jj,
                              x = weights,
                              dims = c(n, n))
    ret <- list(A = A, V = as.character(vertices))
    return(ret)
}
