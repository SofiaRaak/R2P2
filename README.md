# RNA-seq Read count Processing Pipeline: R2P2

RNA-seq analysis includes many often labour intensive steps. This pipeline offers a one-stop solution to differential gene expression analysis of annotated count data through a single function. The function takes a limma-vloom approach to modelling followed by empirical Bayes moderation. The package was developed and tested in R 3.6.1.

## Getting started

### Input files
This pipeline assumes raw RNA-seq data has already been mapped an annotated with a gene identifier, e.g. ENSEMBLE, EntrezID, or NCBI accession (gene symbol probably will not work in case of duplications). The count data should be in a .csv file with column names for gene identifiers and sample names:

```
ENSEMBL   SRR5079936  SRR5079937  SRR5079938  SRR5079939  SRR5079940  SRR5079941
ENSMUSG00000000001  90.485747 104.407157  96.732121 38.5074519  76.0860922  84.1125044
ENSMUSG00000000003  0.000000   0.000000   0.000000  0.0000000   0.0000000   0.0000000
ENSMUSG00000000028  1.074343   4.309391   3.461556  2.1408018   1.0712002   2.1430561
ENSMUSG00000000037  7.521970   2.666617   6.954252  0.7552044   0.8874543   0.7569267
ENSMUSG00000000049  0.000000   0.000000   0.000000  0.0000000   0.9869715   0.0000000
ENSMUSG00000000056  159.414126 116.012337 314.379830 61.9002918 140.0571995 109.7252262
```

The pipeline also requires some metadata, also in .csv format. The metadata should contain a minimum of 3 columns with the following headings: ID, containing the sample numbers of the RNA-seq experiment; Genotype, containing the genotype of each sample; and Treatment, containing any treatment of the cells used prior to RNA extraction:

```
ID  Genotype  Treatment
SRR5079936  WT  FTY720
SRR5079937  WT  FTY720
SRR5079938  WT  FTY720
SRR5079939  WT  FTY720_LPS
SRR5079940  WT  FTY720_LPS
SRR5079941  WT  FTY720_LPS
```

In the metadata file, the Genotype and Treatment columns may be left empty, but if comparisons between groups are necessary at least one of the columns must contain two or more distinct values.
NB! Values in the metadata file CANNOT start with a number!


### Prerequisites

This pipeline requires several R-packages to successfully run. They can be installed with the following code:

```
install.packages(c("devtools","limma","Glimma","glpots","RColorBrewer","randomcoloR","dplyr","edgeR"))
```

### Installing

This pipeline can be installed through directly from github as follows:

```
devtools::install_github("SofiaRaak/R2P2")
```

Alternatively, the pipeline can be downloaded and installed directly:

```
devtools::install("{PATH_TO\\R2P2-master\\R2P2}")
```

The pipeline can then be loaded as a normal R package.

## Using the pipeline


### Options

The pipeline has a number of options. `sourcedir ` is the path to the directory containig the input files. `countMatrix` and `metadata` are both the filenames and extension of the count data and metadata files respectively. `rownameColumn` is the name of the column containing gene identifiers.  `threshold` is the minimum count per million included in downstream analysis, *e.g.* `threshold=0.5` only includes transcripts with counts exceeding 500 000 counts per million. `contrasts`, `figures` and `libsize` are all logicals determining whether comparisons between groups are performed, figures of plots and library sizes are saved. `destdir` is the path to the desired output directory.

Example pipeline options:

```
process_countMatrix(sourcedir="C:\\Documents\\SampleData\\",
                    countMatrix="counts.csv",
                    metadata="metadata.csv",
                    rownameColumn="ENSEMBL",
                    threshold=0.5,
                    contrasts=TRUE,
                    figures=TRUE,
                    libsize=FALSE,
                    destdir="C:\\Documents\\SampleData\\DGEanalysis\\"
```


### Outputs

The pipeline always outputs a .csv file containing the Log2 counts per million (log2cpm.csv). If `contrasts=TRUE` is used, the pipeline also outputs .csv files of comparisons between two groups (in all possible combinations) with log2FC, estimated t-values and associate p-values. If `figures=TRUE` is used, the pipeline also outputs the following .pngs:

* Barplot of library sizes of each sample
* Boxplots of unnormalised and normalised count data
* A volcano plot highlighting expression differences over 1 times log2fold change
* A heat map of the top 500 most variable genes.

NB! To group the samples in the heatmap on Genotype or Treatment, be sure to organise the metadata file with the samples grouped together.

Additionally, if `contrasts=TRUE` is used `figures=TRUE` also outputs .pngs of mean difference (MD) plots.

`libsize=TRUE` outputs a .csv file of the library size of each sample.



## Built With

* [roxygen2](https://github.com/r-lib/roxygen2)
* [devtools](https://github.com/r-lib/devtools)

Developed and tested in R v. 3.6.1 for Windows (64-bit)

## Authors

* **Sofia Raak**

This pipeline was developed for use in research producing a dissertation in partial fulfillment of the requirement for the Masters of Research in Stem Cell Neurobiology at School of Biosciences, Cardiff University (2020). Provisional title: *Investigating the role of microglia in excitotoxicity and neurodegeneration*

### Citation

If used in research publications, please cite

Raak, S.B. 2020. *RNA-seq Read count Processing Pipeline: R2P2*. [https://github.com/SofiaRaak/R2P2](https://github.com/SofiaRaak/R2P2)

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details

## Acknowledgments

This pipeline is based on the excellent guide to differential gene expression analysis by [Phipson *et al.* (2016)](https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html).
