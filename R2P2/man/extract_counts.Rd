% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_counts.R
\name{extract_counts}
\alias{extract_counts}
\title{Count data extraction}
\usage{
extract_counts(seqdata, sampleinfo, rownameColumn)
}
\arguments{
\item{seqdata}{The seqdata dataframe}

\item{sampleinfo}{The metadata dataframe}

\item{rownameColumn}{The name of the column containing the gene identifiers,
NB! CANNOT be gene symbols to avoid duplication}
}
\value{
processed dataframe of gene counts
}
\description{
Extracts count data from the seqdata dataframe
}
\examples{
countdata <- extract_counts(seqdata=seqdata,sampleinfo=sampleinfo,rownameColumn="ENSEMBL")

}
\keyword{}
