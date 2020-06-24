#' Count data extraction
#'
#' Extracts count data from the seqdata dataframe
#'
#' @param seqdata The seqdata dataframe
#' @param sampleinfo The metadata dataframe
#' @param rownameColumn The name of the column containing the gene identifiers,
#' NB! CANNOT be gene symbols to avoid duplication
#'
#' @keywords
#'
#' @returns processed dataframe of gene counts
#'
#' @examples
#' countdata <- extract_counts(seqdata=seqdata,sampleinfo=sampleinfo,rownameColumn="ENSEMBL")
#'
#' @export


extract_counts <- function(seqdata,sampleinfo,rownameColumn){
  data <- cbind(NA)
  IDs <- c(paste(sampleinfo$ID))
  for (i in 1:length(IDs)){
    data <- cbind(data,seqdata[grep(paste0(IDs[i]),colnames(seqdata))])
  }
  data <- cbind(data[2:length(data)])
  rownames(data) <- seqdata[,grep(paste0(rownameColumn),colnames(seqdata))]
  data[,1:length(data)] <- sapply(data[,1:length(data)],as.numeric)
  return(data)
}