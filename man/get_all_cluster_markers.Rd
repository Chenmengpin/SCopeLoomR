% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loom.R
\name{get_all_cluster_markers}
\alias{get_all_cluster_markers}
\title{get_all_cluster_markers}
\usage{
get_all_cluster_markers(
  loom,
  log.fc.threshold = 1.5,
  adj.p.value = 0.05,
  n.signif.digits = NULL
)
}
\arguments{
\item{loom}{The loom file handler.}

\item{log.fc.threshold}{The log fold change threshold}

\item{n.signif.digits}{Number of significant digits in the output.}

\item{adj.p.value.threshold}{The adjusted p-value threshold}
}
\value{
A list of data.frames containing all the clusterings with their cluster markers (adj. p-value and log fold changes). The list names refer to the different clustering names.
}
\description{
Get cluster markers from all clusterings fromthe given loom file with filter applied by the given log.fc.threshold and adj.p.value.threshold.
}
