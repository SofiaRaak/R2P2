#' Import metadata
#'
#'
#' Imports metadata in .csv format, needs a minimum of 3 columns with the following headings:
#'
#' 1: ID: Sample names corresponding to the sample names in the count matrix.
#'
#' 2: Genotype
#' May be left blank. Genotype of sample, e.g. WT or SOD1
#'
#' 3: Treatment
#' May be left blank. Treatment of sample, e.g. drug or placebo
#'
#'
#' @param sourcedir Path to directory containing metadata in .csv format
#' @param metadata Name of .csv file containing metadata
#'
#' @keywords metadata
#'
#' @returns a dataframe containing metadata of the experiment
#'
#' @examples
#' metadata <- import_metadata(sourcedir="C:/Documents/",metadata="metadata.csv")
#'
#' @export


import_metadata <- function(sourcedir,metadata){
  temp <- read.csv(paste0(sourcedir,metadata))
  return(temp)
}
