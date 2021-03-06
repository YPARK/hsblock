% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{build.random.sparse.Z}
\alias{build.random.sparse.Z}
\title{Random sparse matrix of latent assignment matrix}
\usage{
build.random.sparse.Z(K, n)
}
\arguments{
\item{K}{number of groups}

\item{n}{number of vertices}
}
\value{
[K x n] matrix
}
\description{
Random sparse matrix of latent assignment matrix
}
\examples{
Z <- build.random.sparse.Z(20, 10)
library(Matrix)
image(Matrix(Z))

}
