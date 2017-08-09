\name{plotSamplesByDipRegion}
\alias{plotSamplesByDipRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots expression distributions of genes selected randomly from each dip output region
}
\description{
Selects a user-defined number of sample genes from across the spectrum of dip values and plots their RNAseq expression distributions.
}
\usage{
plotSamplesByDipRegion(minimumCounts, breaks, labels,
rawRNAdata, RNAdata, xlab, samples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{minimumCounts}{
If 'rawRNAdata' is supplied, 'minimumCounts' is the threshold that each gene's maximum raw expression value must exceed to remain in the normalized RNA data for the analysis.
}
  \item{breaks}{
'breaks' is an optional input vector that defines the number and location of the lines separating regions of dip values. There will always be one more region than break points. Similarly, the labels vector must be one value longer than the breaks vector. If left blank, the samples will be selected from all remaining genes with no respect to dip value. A maximum of four break points and five labels may be specified.
}
  \item{labels}{
'labels' is an input vector that gives names to each of the regions. If 'breaks' is unspecified, it need not be defined. However, if 'breaks' is supplied, the resulting regions will require names, which are supplied by 'labels'. 'labels' must be of length one greater than 'breaks'. It may take simple forms like numbers (ex: labelInput <- c(1:8)), which are recommended due to their inherent ordinality, or more complex forms like strings (ex: labelInput <- c("Low Dip", "Moderate Dip", "High Dip")). A maximum of four break points and five labels may be specified.
}
  \item{rawRNAdata}{
This is an optional argument used only when a 'minimumCounts' filter is to be applied. Each gene's highest expression level is extracted from 'rawRNAdata'. If that maximum expression does not exceed the value supplied by 'minimumCounts', then that gene will be exempted from the analysis of the normalized counts. It is crucial to note that this is not the dataset to be analyzed. This set serves as part of an optional filter. The rows must be genes (with gene names as row names) and the columns must be the samples, both of which should correspond directly with the rows and columns of the normalized data supplied as 'RNAdata'.
}
  \item{RNAdata}{
This argument specifies the RNAseq dataset to be analyzed. It may be in the form of either raw or normalized data. If the analysis is to involve a 'minimumCounts' screen to filter out low-expression genes, then 'RNAdata' should specify the normalized expression data. The rows must be genes (with gene names as row names) and the columns must be different samples (ideally, with the sample names as the column names, but this specific exemption will not disable the program).
}

  \item{xlab}{
Provides a standardized x-axis label for the output plots.
}
  \item{samples}{
This specifies the number of samples to be selected from each region and may range from 1-6. Using a smaller number of samples helps with reliability and '4' is the recommended, and default, value for 'samples'.
}
}
\details{
The sampling starts to break down at the upper end of the dip distribution. Typically, the plots become unreliable when the number of samples exceeds 4 and when the dip values of the samples in the region exceed 0.12.
}
\value{
Plots are returned from this function grouped by their region of gene origin. The number of regions determines how many sets of plots are returned. If only one region is used, the output plot is called 'x1'. If two, 'x1' and x2'. The same applied up through five regions, at which point the five plots will be titled 'x1', 'x2', 'x3', 'x4', and 'x5'. The number on the plot name indicates which region number the samples are pulled from: dip value ascending with title number. Use 'plot(x1)' to view for each respective plot.
}
\references{
Martin Maechler (2016). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.75-7. https://CRAN.R-project.org/package=diptest

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
}
\author{
Software authors: Jeremy Sieker, Sohyon Lee, Kristin Baldwin

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
