% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_file.R
\name{get_file}
\alias{get_file}
\title{Get path to playdata file in inst/extdata}
\usage{
get_file(file = NULL)
}
\arguments{
\item{file}{string. Name of file. If `NULL`, the example files will
be listed.}
}
\description{
`playdata` comes bundled with a number of non-RDA files in its `inst/extdata`
directory. This function make them easy to access. This function was
taken from tidyverse/readr.
}
\examples{
get_file()
get_file("kinase_substrates_kea.gmt")
}
