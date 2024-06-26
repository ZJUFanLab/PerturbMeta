% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perturbmeta.R
\name{perturbmeta}
\alias{perturbmeta}
\title{Calculate PLS-DA and NMF for each cluster}
\usage{
perturbmeta(
  count,
  sgrna,
  factor = 20,
  nfeatures = 3000,
  npcs = 20,
  filter = 0.5,
  ctrl = "NT",
  VIP = 1,
  cluster = 5,
  ncores = 0.5
)
}
\arguments{
\item{count}{A gene count matrix.}

\item{sgrna}{A cell-perturbation data frame. The colname of sgrna data frame must be 'cell' and 'perturbation'.}

\item{factor}{The number of factor in perturbmeta Default: 20.}

\item{nfeatures}{A parameter for extracting high-variable genes. Default: 3000.}

\item{npcs}{A number of PCs in PCA. Default: 20 (Same as number of factor).}

\item{filter}{The ratio of filtering cells that are perturbed non-significant. Default: 0.5.}

\item{ctrl}{The name of control perturbation. Default: NT.}

\item{VIP}{The threshold of extracting features from PLS-DA. Default: 1.}

\item{cluster}{The number of clusters that used for perturbmeta.}

\item{ncores}{The number of cores used for calculating NMF. Default: 0.5.}
}
\value{
A list contain the result of perturbmeta for each cluster
}
\description{
Filter the cells that are non-significant perturbed, then carry out NMF, extract the cell-embedding and feature-loading for
re-clustering analysis. And perform PLS-DA for each cluster, then extract high-contribution genes for further NMF.
}
\examples{
res <- perturbmeta(count, sgrna,factor=20, nfeatures=3000, ctrl='NT', cluster=5)
}
