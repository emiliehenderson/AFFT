% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{GetBandIndices}
\alias{GetBandIndices}
\title{Copies source airphotos to intermediate directory, calculates extra indices to add to image. Resolution unchanged.}
\usage{
GetBandIndices(
  filelist,
  indfuns = indexFuns,
  outpath = "1_intermediate",
  ncpu = 4
)
}
\arguments{
\item{filelist}{character vector with tif image names. If full path not included, GetMetrics will assume that files are in 0_raw folder.}

\item{indfuns}{file path pointer where R will store intermediate (full resolution) image tiles, with calculated indices as well as original bands (although not blue currently)}
}
\value{
filepaths pointing to single-band rasters nested in subfolders inside 1_intermediate.
}
\description{
Copies source airphotos to intermediate directory, calculates extra indices to add to image. Resolution unchanged.
}
