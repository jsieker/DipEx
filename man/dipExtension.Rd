\name{dipExtension}
\alias{dipExtension}

\title{
Assigns 'Region' values to break up the dip distribution
}
\description{
Calculates the dip value for each gene's distribution and then separates those genes into user-defined groups. If using normalized RNA expression values, the corresponding raw RNA counts can also be used to apply a 'minimum counts' filter, exempting those genes that do not pass the user-defined threshold.
}
\usage{
dipExtension(breaks, RNAdata, rawRNAdata, minimumCounts)
}
\arguments{
  \item{breaks}{
'breaks' is an optional input vector that defines the number and location of the lines separating regions of dip values. There will always be one more region than break points. Similarly, the labels vector must be one value longer than the breaks vector. If left blank, the outputs will all be calculated (Zero-excluded median, Zero-included median, Dip value, Dip P value, and GeneID), but the genes will not be split up into separate regions.
}

  \item{RNAdata}{
This argument specifies the RNAseq dataset to be analyzed. It may be in the form of either raw or normalized data. If the analysis is to involve a 'minimumCounts' screen to filter out low-expression genes, then 'RNAdata' should specify the normalized expression data. The rows must be genes (with gene names as row names) and the columns must be different samples (ideally, with the sample names as the column names, but this specific exemption will not disable the program). Columns with non-numerical data (or containing data not relating to a sample) should be specifically exempted before any analysis is attempted.
}
  \item{rawRNAdata}{
This is an optional argument used only when a 'minimumCounts' filter is to be applied. Each gene's highest expression level is extracted from 'rawRNAdata'. If that maximum expression does not exceed the value supplied by 'minimumCounts', then that gene will be exempted from the analysis of the normalized counts. It is crucial to note that this is not the dataset to be analyzed. This set serves as part of an optional filter. The rows must be genes (with gene names as row names) and the columns must be the samples, both of which should correspond directly with the rows and columns of the normalized data supplied as 'RNAdata'. Columns with non-numerical data (or containing data not relating to a sample) should be specifically exempted before any analysis is attempted.
}
  \item{minimumCounts}{
If 'rawRNAdata' is supplied, 'minimumCounts' is the threshold that each gene's maximum raw expression value must exceed to remain in the normalized RNA data for the analysis.
}
}

\value{
The output of this function is a dataframe containing each gene's GeneID, dip information, region assignment, zero-excluded median, and zero-included median.

The zero-excluded median (ZXM) is a value that takes the average of every expression value (raw or normalized) that exceeds zero. It is a measure developed as a metric for how intensely expression is elevated when highly bimodal genes are activated. The zero-included median includes zeroes in that calculation and therefore is not different from the median.
}
\references{
Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest
}
\author{
Software authors: Jeremy Sieker, Sohyon Lee, Kristin Baldwin
}
\examples{
x <- paste("https://www.ebi.ac.uk/gxa/experiments-content",
  "/E-GEOD-70484/resources/BaselineProfilesWriterService.RnaSeq/tsv", sep = "")
ss <- read.table(url(x),sep = '\t', header = TRUE)
ss <- ss[(as.numeric(row.names(ss)) \%\% 2 == 1),]
#if the previous line fails, remove the slashes from the modulo division and try again
row.names(ss) <- ss$Gene.ID #cutting it in half to speed up the examples
ss <- ss[,-c(1:2)] #removing everything that is not expression data

mod2 <- dipExtension(RNAdata = ss, breaks = c(-Inf, 0.05, 0.11, Inf))


#in many datasets, genes with extremely low dip values (or high dip p values)
#will appear in their expression plots as normal distributions with a mean around zero.
#These are typically just genes that don't have registered counts in any of the samples.
#To remove these genes, there are a few options.
#One can apply both normalized and raw counts (i.e.- supplying both the RNAdata
# and rawRNA arguments)
#and employ the 'minimumCounts' filter.
#Alternatively, pre-filter your data for genes that do not pass your
# desired expression threshold, then simply
#use that data for your RNAdata argument and leave the rawRNAdata
# and minimumCounts arguments blank.

#Due to the difficulty of finding publicly available datasets that have paired
# raw and normalized counts, the filter will not be
#demonstrated in this example. However, if you were to have a raw counts set
# called raws and a normalized counts set called ss, the
#code would be along the lines of

#mod2 <- dipExtension(RNAdata = ss, rawRNAdata <- raws, minimumCounts = 50,
#                      breaks = c(-Inf, 0.05, 0.11, Inf))

}
