\name{plotSamplesByDipRegion}
\alias{plotSamplesByDipRegion}

\title{
Plots expression distributions of genes selected randomly from each dip output region
}
\description{
Selects a user-defined number of sample genes from across the spectrum of dip values and plots their RNAseq expression distributions.
}
\usage{
plotSamplesByDipRegion(minimumCounts, breaks,
rawRNAdata, RNAdata, xlab, samples)
}

\arguments{
  \item{minimumCounts}{
If 'rawRNAdata' is supplied, 'minimumCounts' is the threshold that each gene's maximum raw expression value must exceed to remain in the normalized RNA data for the analysis.
}
  \item{breaks}{
'breaks' is an optional input vector that defines the number and location of the lines separating regions of dip values. There will always be one more region than break points. Similarly, the labels vector must be one value longer than the breaks vector. If left blank, the samples will be selected from all remaining genes with no respect to dip value. A maximum of four break points and five labels may be specified.
}

  \item{rawRNAdata}{
This is an optional argument used only when a 'minimumCounts' filter is to be applied. Each gene's highest expression level is extracted from 'rawRNAdata'. If that maximum expression does not exceed the value supplied by 'minimumCounts', then that gene will be exempted from the analysis of the normalized counts. It is crucial to note that this is not the dataset to be analyzed. This set serves as part of an optional filter. The rows must be genes (with gene names as row names) and the columns must be the samples, both of which should correspond directly with the rows and columns of the normalized data supplied as 'RNAdata'.
}
  \item{RNAdata}{
This argument specifies the RNAseq dataset to be analyzed. It may be in the form of either raw or normalized data. If the analysis is to involve a 'minimumCounts' screen to filter out low-expression genes, then 'RNAdata' should specify the normalized expression data. The rows must be genes (with gene names as row names) and the columns must be different samples (ideally, with the sample names as the column names, but this specific exemption will not disable the program, just make it hard to analyze meaningfully). Columns with non-numerical data (or containing data not relating to a sample) should be specifically exempted before any analysis is attempted.
}

  \item{xlab}{
Provides a standardized x-axis label for the output plots.
}
  \item{samples}{
This specifies the number of gene samples to be selected from each region for plotting and may range from 1-4, with 4 as the default.
}
}
\details{
The sampling starts to break down at the upper end of the dip distribution. Typically, the plots become unreliable when the number of samples exceeds 4 and when the dip values of the samples in the region exceed 0.12.

It is important to note than in the dataframe this function outputs, any gene whose name started with a number (0-9) will have an 'X' placed in front of its name in the output GeneID. Many hours were put into avoiding this outcome, but it was the best solution we could find in the context of ggplot2's treatment of string variables.
}
\value{
The output of this function is a list containing the elements DF (Dip and region output dataframe) and plots for the samples from regions 1-5 (Plot1, Plot2, Plot3, Plot4, Plot5). If fewer than five regions are sampled, the upper regions' plot variables will simply read '0'.
}
\references{
Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
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

mod <- plotSamplesByDipRegion(RNAdata = ss, breaks = c(-Inf, 0.05, 0.11, Inf),
                              samples = 4, xlab = "Sample RNA Expression")

#in many datasets, genes with extremely low dip values (or high dip p values)
#will appear in their expression plots as normal distributions with a mean around zero.
#These are typically just genes that don't have registered counts in any of the samples.
#To remove these genes, there are a few options.
#One can apply both normalized and raw counts (i.e.- supplying both the RNAdata
# and rawRNA arguments)
#and employ the 'minimumCounts' filter.
#Alternatively, pre-filter your data for genes that do not pass your desired
# expression threshold, then simply
#use that data for your RNAdata argument and leave the rawRNAdata and
# minimumCounts arguments blank.

#Due to the difficulty of finding publicly available datasets that have paired
# raw and normalized counts, the filter will not be
#demonstrated in this example. However, if you were to have a raw counts set called
# raws and a normalized counts set called ss, the
#code would be along the lines of...

#mod <- plotSamplesByDipRegion(RNAdata = ss, rawRNAdata = raws, minimumCounts = 50,
#                              breaks = c(-Inf, 0.05, 0.11, Inf), samples = 4,
#                              xlab = "Sample RNA Expression")



}
