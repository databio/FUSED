% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FUSED.R
\name{constructORIEN}
\alias{constructORIEN}
\title{Modify a valid ORIEN fusion data.table for downstream analysis. ORIEN uses
a merged output derived from both STAR and ARRIBA fusion predictions.}
\usage{
constructORIEN(x, window = 99)
}
\arguments{
\item{x}{A valid data.table object containing ORIEN predicted fusions.}

\item{window}{A numeric value by which to expand the region around the
breakpoint}
}
\value{
A fusion data.table object
}
\description{
Modify a valid ORIEN fusion data.table for downstream analysis. ORIEN uses
a merged output derived from both STAR and ARRIBA fusion predictions.
}
