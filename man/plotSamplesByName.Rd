\name{plotSamplesByName}
\alias{plotSamplesByName}

\title{
Plots expression distributions of genes selected by name
}
\description{
Allows the user to plot the expression distributions of up to 8 genes of their choice.
}
\usage{
plotSamplesByName(nameList, RNAdata)
}

\arguments{
  \item{nameList}{
The names of the desired genes should be concatenated into this variable. If more than eight gene names are specified, function will typically fail. (ex. li <- c("Ascl3", "Bambi", "Pols2")) ). Pay very close attention to proper capitalization in specifying gene names. Any error in capitalization or erroneously-placed spaces will cause the function to fail.
}
  \item{RNAdata}{
This argument specifies the RNAseq dataset from which the specified genes' expressional data will be taken. It may take the form of raw or normalized data. Columns with non-numerical data (or containing data not relating to a sample) should be specifically exempted before any analysis is attempted.
}
}

\value{
A single plot is the output of this function.
}
\references{
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

geneList <- row.names(ss[sample(nrow(ss), 4), ]) #get some random gene names
mod4 <- plotSamplesByName(nameList = geneList, RNAdata = ss)
}

