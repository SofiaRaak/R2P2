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
  vector <- c()
  for (i in 1:length(sampleinfo$ID)){
    vector <- c(vector,paste(sampleinfo$ID[i]))
  }
  d <- (length(colnames(seqdata))-length(vector))+1
  data <- seqdata[,d:length(seqdata)]
  rownames(data) <- seqdata[,grep(paste0(rownameColumn),colnames(seqdata))]
  return(data)
}
