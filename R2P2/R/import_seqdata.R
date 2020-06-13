#' Import seqdata
#'
#' Imports RNA-seq count matrix in .csv form
#' NB! Make sure the column with the gene identifiers (e.g. ENSEMBL ID, EntrezID,
#' NCBI accession number) is named!
#'
#' @param countMatrix Path to directory containing your count matrix in .csv format
#' @param seqdata Name of your .csv file
#'
#' @keywords Count matrix
#'
#' @return a dataframe of count data from a .csv file
#'
#' @examples
#' seqdata <- import_seqdata(sourcedir="C:/Documents/",countMatrix="countMatrix.csv")
#'
#'  @export
#'


import_seqdata <- function(sourcedir,countMatrix){
  temp <- read.csv(paste0(sourcedir,countMatrix))
  return(temp)
}
