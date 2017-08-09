\name{dipExtension}
\alias{dipExtension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assigns 'region' values to break up the dip distribution
}
\description{
Calculates the dip value for each gene's distribution and then separates those genes into user-defined groups. If using normalized RNA expression values, the corresponding raw RNA counts can also be used to apply a 'minimum counts' filter, exempting those genes that do not pass the user-defined threshold.
}
\usage{
dipExtension(breaks, labels, RNAdata, rawRNAdata, minimumCounts)
}
\arguments{
  \item{breaks}{
'breaks' is an optional input vector that defines the number and location of the lines separating regions of dip values. There will always be one more region than break points. Similarly, the labels vector must be one value longer than the breaks vector. If left blank, the outputs will all be calculated (Zero-excluded median, Zero-included median, Dip value, Dip P value, and GeneID), but the genes will not be split up into separate regions.
}
  \item{labels}{
'labels' is an input vector that gives names to each of the regions. If 'breaks' is unspecified, it need not be defined. However, if 'breaks' is supplied, the resulting regions will require names, which are supplied by 'labels'. 'labels' must be of length one greater than 'breaks'. It may take simple forms like numbers (ex: labelInput <- c(1:8)), which are recommended due to their inherent ordinality, or more complex forms like strings (ex: labelInput <- c("Low Dip", "Moderate Dip", "High Dip")).
}
  \item{RNAdata}{
This argument specifies the RNAseq dataset to be analyzed. It may be in the form of either raw or normalized data. If the analysis is to involve a 'minimumCounts' screen to filter out low-expression genes, then 'RNAdata' should specify the normalized expression data. The rows must be genes (with gene names as row names) and the columns must be different samples (ideally, with the sample names as the column names, but this specific exemption will not disable the program).
}
  \item{rawRNAdata}{
This is an optional argument used only when a 'minimumCounts' filter is to be applied. Each gene's highest expression level is extracted from 'rawRNAdata'. If that maximum expression does not exceed the value supplied by 'minimumCounts', then that gene will be exempted from the analysis of the normalized counts. It is crucial to note that this is not the dataset to be analyzed. This set serves as part of an optional filter. The rows must be genes (with gene names as row names) and the columns must be the samples, both of which should correspond directly with the rows and columns of the normalized data supplied as 'RNAdata'.
}
  \item{minimumCounts}{
If 'rawRNAdata' is supplied, 'minimumCounts' is the threshold that each gene's maximum raw expression value must exceed to remain in the normalized RNA data for the analysis.
}
}

\value{
The output of this function is a dataframe containing each gene's GeneID, dip information, region assignment, zero-excluded median, and zero-included median.

The zero-excluded median (ZXM) is a value that takes the average of every expression value (raw or normalized) that exceeds zero. It is a measure developed as a metric for how intensely expression is elevated when highly bimodal genes are activated. The zero-included median includes zeroes in that calculation.
}
\references{
Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest
}
\author{
Software authors: Jeremy Sieker, Sohyon Lee, Kristin Baldwin

}
