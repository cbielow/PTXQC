% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mzQC.R
\name{writeMZQC}
\alias{writeMZQC}
\title{Writes a full mzQC object to disk}
\usage{
writeMZQC(filepath, mzqc_obj)
}
\arguments{
\item{filepath}{A filename (with path) to write to. Should have '.mzQC' as suffix.}

\item{mzqc_obj}{An mzQC object, which is serialized to JSON and then written to disk}
}
\description{
Writes a full mzQC object to disk
}
