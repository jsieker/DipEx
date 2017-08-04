\name{previewDipDistribution}
\alias{previewDipDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the distribution of the gene set's dip values and dip p-values.
}
\description{
Allows the viewing of where the dataset's dip values and dip p-values are distributed to aid in the assignment of region break lines and the number of regions.
}
\usage{
previewDipDistribution(RNAdata, rawRNAdata, minimumCounts, barLine)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RNAdata}{
This argument specifies the RNAseq dataset to be analyzed. It may be in the form of either raw or normalized data. If the analysis is to involve a 'minimumCounts' screen to filter out low-expression genes, then 'RNAdata' should specify the normalized expression data. The rows must be genes (with gene names as row names) and the columns must be different samples (ideally, with the sample names as the column names, but this specific exemption will not disable the program).
}
  \item{rawRNAdata}{
This is an optional argument used only when a 'minimumCounts' filter is to be applied. Each gene's highest expression level is extracted from 'rawRNAdata'. If that maximum expression does not exceed the value supplied by 'minimumCounts', then that gene will be exempted from the analysis of the normalized counts. It is crucial to note that this is not the dataset to be analyzed. This set serves as part of an optional filter. The rows must be genes (with gene names as row names) and the columns must be the samples, both of which should correspond directly with the rows and columns of the normalized data supplied as 'RNAdata'.
}
  \item{minimumCounts}{
If 'rawRNAdata' is supplied, 'minimumCounts' is the threshold that each gene's maximum raw expression value must exceed to remain in the normalized RNA data for the analysis.
}
  \item{barLine}{
This selects the x-intercept of the bar that can be overlaid on the graph. It will default to 0 if no value is supplied.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a dataframe (DipOutputDF) with the each gene's dip values and dip p-values, a density plot of the genes' dip values (x1), and a density plot of the genes' dip p-values (x2).
}
\references{
Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
}
\author{
Software authors: Sohyon Lee, Jeremy Sieker, Kristin Baldwin
}
