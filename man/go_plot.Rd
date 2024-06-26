% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/go_plot.R
\name{go_plot}
\alias{go_plot}
\title{Plot the result of GO enrichment analysis for interested cluster}
\usage{
go_plot(go_res, score, p, cluster = 5, top = 5, factor = 20, plot = 1)
}
\arguments{
\item{go_res}{The result from function enrich_go.}

\item{score}{The first result of function perturbmeta_score.}

\item{p}{The second result of function perturbmeta_score .}

\item{cluster}{The number of clusters that used for perturbmeta.}

\item{top}{The number of top GO pathways. Default: 1.}

\item{factor}{The number of factor in perturbmeta Default: 20.}

\item{plot}{A parameter for which cluster to plot dotplot.}
}
\value{
A dotplot for GO enrichment analysis for interested cluster
}
\description{
Use dotplot to plot the result of GO enrichment analysis. It can plot the GO enrichment results (p < 0.05) for interested cluster
by setting the plot parameter  to the number of interested cluster.
}
\examples{
p<-go_plot(go_res,res_score[[1]],res_score[[2]],top=5,plot=1)
}
