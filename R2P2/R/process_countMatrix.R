#' Complete pipeline for processing RNA-seq data in count matrix
#'
#' Differential expression analysis with limma-voom of RNA-seq data in count matrix
#' Based on pipeline by Phipson et al., 2016 (RNA-seq analysis in R)
#'
#' @param sourcedir Path to directory containing count matrix and metadata in .csv format
#' @param countMatrix Name of .csv file containing count matrix
#' NB! In the count matrix file, column of gene identifier (e.g. ENSEMBL ID, EntrezID,
#' NCBI accession number) needs to be named.
#' @param metadata Name of .csv contianing metadata.
#' NB! Metadata file requires a minimum of the following three columns:
#' 1: ID: Sample names corresponding to the sample names in the count matrix.
#' 2: Genotype: May be left blank. Genotype of sample, e.g. WT or SOD1
#' 3: Treatment: May be left blank. Treatment of sample, e.g. drug or placebo
#' @param rownameColumn Name of column in count matrix file contianing gene identifiers
#' @param threshold Positive real number; threshold value of minimum counts per million
#' included in analysis
#' @param contrasts Logical. Determines whether logfold changes between two groups are calculated
#' with empirical Bayes shrinkage of variances and estimation of moderated t-statistic and p-values
#' @param figures Logical. Determines whether .png files of plots are saved
#' @param libsize Logical. Determines whether a .csv file of the library sizes for each sample in the
#' data set is saved
#' @param destdir Path to desired output directory
#'
#' @return A .csv file of normalised Log2 counts per million, plus optional contrast .csv files and figures
#'
#' @references
#'
#' Phipson, B., Trigos, A., Ritchie, M, Doyle, M, Dashnow, H, Law, C. (2016).
#' RNA-seq analysis in R. \url{https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html}
#'
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015).
#' limma powers differential expression analyses for RNA-sequencing and microarray studies.
#' Nucleic Acids Research 43(7), e47.
#'
#' Su, S., Law, C. W., Ah-Cann, C., Asselin-Labat, M. L., Blewitt, M. E., & Ritchie, M. E. (2017).
#' Glimma: interactive graphics for gene expression analysis. Bioinformatics, 33(13), 2050-2052.
#'
#' Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber, Andy Liaw,
#' Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2020).
#' gplots: Various R Programming Tools for Plotting Data. R package version 3.0.3.
#' https://CRAN.R-project.org/package=gplots
#'
#' Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes. R package version 1.1-2.
#' https://CRAN.R-project.org/package=RColorBrewer
#'
#' Ron Ammar (2019). randomcoloR: Generate Attractive Random Colors. R package version 1.1.0.1.
#' https://CRAN.R-project.org/package=randomcoloR
#'
#' Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2020).
#' dplyr: A Grammar of Data Manipulation. R package version 1.0.0. https://CRAN.R-project.org/package=dplyr
#'
#' Robinson MD, McCarthy DJ and Smyth GK (2010).
#' edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.
#' Bioinformatics 26, 139-140
#'
#' McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation.
#' Nucleic Acids Research 40, 4288-4297
#'
#' @keywords
#'
#' @examples
#' process_countMatrix(sourcedir="C:/Documents/",countMatrix="countMatrix.csv",
#' metadata="metadata.csv",rownameColumn="ENSEMBL",threshold=0.5,contrasts=TRUE,
#' figures=TRUE,libsize=FALSE,destdir="C:/Documents/Analysed/")
#'
#' @import limma
#' @import Glimma
#' @import gplots
#' @import RColorBrewer
#' @import randomcoloR
#' @import dplyr
#' @import edgeR
#'
#' @export



process_countMatrix <- function(sourcedir,countMatrix,metadata,rownameColumn,
                                threshold,contrasts=TRUE,figures=TRUE,libsize=FALSE,
                                destdir){

  seqdata <- tryCatch(import_seqdata(sourcedir=sourcedir,countMatrix=countMatrix),
                      error=function(e) print(paste("Importation of count matrix failed, please try again.")))

  sampleinfo <- tryCatch(import_metadata(sourcedir=sourcedir,metadata=metadata),
                         error=function(e) print(paste("Importation of metadata failed, please try again")))

  countdata <- extract_counts(seqdata=seqdata,sampleinfo=sampleinfo,
                              rownameColumn=rownameColumn)

  myCPM <- cpm(countdata)
  thresh <- myCPM > threshold
  keep <- rowSums(thresh) >=2
  counts.keep <- countdata[keep,]

  y <- DGEList(counts.keep)
  logcounts <- cpm(y,log=TRUE)
  y <- calcNormFactors(y)


  sampleinfo$group <- paste0(sampleinfo$Genotype,"_",sampleinfo$Treatment)
  if (contrasts==TRUE){
    group <- sampleinfo$group
    design <- tryCatch(model.matrix(~0+group),
                       error=function(e) print(paste("Contrasts requires at least two distinct entries in
                                                     sampleinfo$Genotype or sampleinfo$Treatment. Please try again.")))
    v <- voom(y,design,plot=FALSE)
    fit <- lmFit(v)

    comb <- c(unique(group))
    temp <- c()
    for (i in 1:length(comb)){
      for (k in i+1:length(comb)){
        if (k <= length(comb)){
          temp <- c(temp,paste0(comb[i],"-",comb[k]))
        }
      }
    }

    cont.matrix <- makeContrasts(contrasts=temp,levels=comb)
    fit.cont <- contrasts.fit(fit,cont.matrix)
    fit.cont <- eBayes(fit.cont)
    summa.fit <- decideTests(fit.cont)

    for (i in 1:length(temp)){
      limma.res <- topTable(fit.cont,coef=temp[i],sort.by="p",n="Inf")
      write.csv(limma.res,file=paste0(destdir,temp[i],".csv"),row.names=TRUE)
    }

    if (figures==TRUE){
      for (i in 1:length(temp)){
        png(filename=paste0(destdir,temp[i],"_MD.png"))
        plotMD(fit.cont,coef=1,status=summa.fit[,paste0(temp[i])],values=c(-1,1),main=paste(temp[i]))
        dev.off()
      }
    }
  }

  if (libsize==TRUE){
    write.csv(y$samples$lib.size,file=paste0(destdir,"libsize.csv"),row.names=TRUE)
  }

  write.csv(logcounts,file=paste0(destdir,"Log2cpm.csv"),row.names=TRUE)

  if (figures==TRUE){

    png(filename=paste0(destdir,"libsize.png"))
    barplot(y$samples$lib.size,las=2,main="Barplot of library sizes")
    dev.off()

    png(filename=paste0(destdir,"boxplot_datanorm.png"))
    par(mfrow=c(1,2))
    boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
    abline(h=median(logcounts),col="blue")
    boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
    abline(h=median(v$E),col="blue")
    dev.off

    png(filename=paste0(destdir,"volcano.png"))
    volcanoplot(fit.cont,coef=1,highlight=100,names=rownames(fit.cont))
    dev.off

    var_genes <- apply(logcounts,1,var)
    select_var <- names(sort(var_genes,decreasing=TRUE))[1:500]
    highly_variable_lcpm <- logcounts[select_var,]
    mypalette <- brewer.pal(11,"RdYlBu")
    morecols <- colorRampPalette(mypalette)
    cols <- c()
    col.cell <- c()
    for (i in 1:length(group)){
      cols <- c(paste(distinctColorPalette(k=3,altCol=FALSE,runTsne=FALSE)))
      for (j in 1:length(cols)){
        col.cell <- c(col.cell,paste(replicate(3,cols[j])))
      }
      return(col.cell)
    }


    png(filename=paste0(destdir,"heatmap.png"))
    heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",
              main="Top 500 most variable genes across samples",
              ColSideColors=col.cell,scale="row",
              Colv=FALSE, dendrogram="row")
    dev.off
  }

}
