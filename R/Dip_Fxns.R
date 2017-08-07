
dipExtension <- function(breaks, labels, RNAdata, rawRNAdata, minimumCounts){

  x.1 <- 0
  x.2 <- 0
  x.3 <- 0
  x.4 <- 0
  if(missing(breaks)) {
    x.1 <- 0; print("cutting = FALSE")
  } else {
    x.1 <-1
  }
  if(missing(labels)) {
    x.2 <- 0
  } else {
    x.2 <-1
  }
  if(missing(rawRNAdata)) {
    x.3 <- 0; print("filtering = FALSE")
  } else {
    x.3 <-1
  }
  if(missing(minimumCounts)) {
    minimumCounts <- 0
  }
  print("reading in data")
  if(x.1==1) {if(length(breaks) != (length(labels) + 1)) {print("invalid dimensions of breaks and labels")}}
  if(mode(RNAdata)=="character") {
  	RNAdata <- as.matrix(read.table(RNAdata))
  	
  }
  
 
  RNAdataDF <- data.frame(RNAdata)


  #filtering of RNAdataDF and DipOutputDF. If the DipOutputDF one fails, move this filter in front of the DipOutputDF creation.
  if(x.3 == 1) {	print("filtering")
  	rawRNAdata <- read.table(rawRNAdata)
    #row.names(DipOutputDF) <- row.names(RNAdataDF)
    rawRNAdata$max <- 0
    rawRNAdata$max <- apply(rawRNAdata, 1, max)
    RNArawbelow <- rawRNAdata[which(rawRNAdata$max<minimumCounts),]
    rr <- row.names(RNArawbelow)
    RNAdataDF_R <- RNAdataDF[!(row.names(RNAdataDF) %in% rr),]
    RNAdataMat <- as.matrix(RNAdataDF_R)
    DipOutput <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
    DipOutputPvals <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
    #DipOutputDF <- DipOutputDF[!(row.names(DipOutputDF)) %in% rr,]
    #row.names(DipOutputDF) <- c(1:nrow(DipOutputDF))
  }
  if(x.3 == 0) {
  	#consider removing/changing redundant variable names
  	RNAdataDF_R <- RNAdataDF
  	RNAdataMat <- as.matrix(RNAdataDF_R)
  	DipOutput <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
    DipOutputPvals <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  }


  print("conducting dip analysis")
  for (k in 1:nrow(RNAdataMat)){	#dip analysis

    e <- diptest::dip(RNAdataMat[k,], full.result = FALSE, min.is.0 = FALSE, debug = FALSE)
    f <- diptest::dip.test(RNAdataMat[k,], simulate.p.value = FALSE, B = 2000)
    DipOutput[k,1] <- e #save output into pre-prepared matrix
    DipOutputPvals[k,1] <- f$p.value
  }

  DipOutputDF <- data.frame(DipOutput)
  DipOutputPDF <- data.frame(DipOutputPvals)
  DipOutputDF$GeneID <- rownames(RNAdataDF_R)
  DipOutputDF$Dip_p_Vals <- DipOutputPDF[,1]


  #cutting
  if(x.1 == 1 && x.2 == 1) {
    print("cutting")
    DipOutputDF$Region <- cut(DipOutputDF$DipOutput,
                              breaks = breaks,
                              labels = labels,
                              right = FALSE)
    DipOutputDF$Region <- as.numeric(as.character(DipOutputDF$Region))

  }	#breaks

  ZeroXMedian <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  ZeroIMedian <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  row.names(ZeroXMedian) <- row.names(RNAdataDF_R)
  colnames(ZeroXMedian) <- c("Zero-excluded Median")
  row.names(ZeroIMedian) <- row.names(RNAdataDF_R)
  colnames(ZeroIMedian) <- c("Zero-included Median")

  print("calculating ZXM and ZIM")
  for(k in 1:nrow(RNAdataDF_R)) {
    shelf <- RNAdataDF_R[k,]	#which(RNAdataDF[k,i] >= 0)
    row.names(shelf) <- c(1)
    shelf2 <- shelf[,(shelf[1,]>0)]
    ZeroXMedian[k,] <- median(as.numeric(shelf2))
    ZeroIMedian[k,] <- median(as.numeric(shelf))
  }


  DipOutputDF$ZeroXMedian <- as.numeric(ZeroXMedian)
  DipOutputDF$ZeroIMedian <- as.numeric(ZeroIMedian)
  DipOutputDF <- dplyr::mutate(DipOutputDF, Divergence = ZeroXMedian-ZeroIMedian)
  DipOutputDF <<- DipOutputDF[,c(2, 1, 3, 4, 5, 6)]
  return(DipOutputDF)

}

plotSamplesByDipRegion <- function(minimumCounts, breaks, labels, rawRNAdata,
                                  RNAdata, xlab, samples) {
  # parameters: filter, break locations (1-4 OR "NONE"), labels,
  # rawRNAdata, normRNAdata, xlab(xx) (make it work from 1 to 5. No default), samples(number to be sampled)

  xx <- xlab
  x.1 <- 0
  if(missing(filter)) {print("no filter applied")}
  if(missing(RNAdata) && missing(rawRNAdata)) {print("RNA data missing"); x.1 <- 1}
  if(missing(breaks)) {print("no breaks selected");  regions <- 1}
  else(regions <- (length(breaks) + 1))
  if(missing(RNAdata) && !missing(rawRNAdata)) {RNAdata <- rawRNAdata}
  if(missing(samples)) {samples <- 4}
  

  if(mode(RNAdata)=="character") {
    RNAdata <- as.matrix(read.table(RNAdata))
  }
  RNAdataDF <- data.frame(RNAdata)
  
  if(mode(rawRNAdata)=="character") {
    rawRNAdata <- as.matrix(read.table(rawRNAdata))
  }
  
  if(!missing(rawRNAdata && !missing(RNAdata))){
  RNArawcounts <- rawRNAdata
  RNArawcounts$max <- 0
  RNArawcounts$max <- apply(RNArawcounts, 1, max)
  RNArawbelow50 <- RNArawcounts[which(RNArawcounts$max<minimumCounts),]
  rr <- row.names(RNArawbelow50)
  RNAdataDF_R <- RNAdataDF[(!row.names(RNAdataDF) %in% rr), ]
  RNAdataMat <- as.matrix(RNAdataDF_R)
  DipOutput <- matrix(nrow=nrow(RNAdataDF_R), ncol=1) #preparing a matrix to receive output in the function step
  DipOutputPvals <- matrix(nrow=nrow(RNAdataDF_R), ncol=1) #preparing a matrix to receive output in the function step
  } else {
  RNAdataDF_R <- RNAdataDF
  RNAdataMat <- as.matrix(RNAdataDF_R)
  DipOutput <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  DipOutputPvals <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
}

  

  #setup and calculation
  for (k in 1:nrow(RNAdataMat)){

    e <- diptest::dip(RNAdataMat[k,], full.result = FALSE, min.is.0 = FALSE, debug = FALSE) #actual dip statistical test
    f <- diptest::dip.test(RNAdataMat[k,], simulate.p.value = FALSE, B = 2000)
    DipOutput[k,1] <- e #save output into pre-prepared matrix
    DipOutputPvals[k,1] <- f$p_value
    ##consider extracting the rownames off RNAdata as a vector, which you can then reimpute here
  }

  DipOutputDF <- data.frame(DipOutput) # matrix -> dataframe
  DipOutputPDF <- data.frame(DipOutputPvals)
  DipOutputDF$GeneID <- rownames(RNAdataDF_R)

  if(regions==1){labels = "Undivided Sample"}
  if(regions==2){labels = c("Region A", "Region B")}
  if(regions==3){labels = c("Region A", "Region B", "Region C")}
  if(regions==4){labels = c("Region A", "Region B", "Region C", "Region D")}
  if(regions==5){labels = c("Region A", "Region B", "Region C", "Region D", "Region E")}


  if(breaks!="NONE" | (is.numeric(breaks) && breaks > 0)){
    DipOutputDF$Region <- cut(DipOutputDF$DipOutput,
                              breaks = breaks,
                              labels = labels,
                              right = FALSE)
    DipOutputDF$Region <- as.numeric(as.character(DipOutputDF$Region))
  }

  ZeroXMedian <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  ZeroIMedian <- matrix(nrow=nrow(RNAdataDF_R), ncol=1)
  row.names(ZeroXMedian) <- row.names(RNAdataDF_R)
  colnames(ZeroXMedian) <- c("Zero-excluded Median")
  row.names(ZeroIMedian) <- row.names(RNAdataDF_R)
  colnames(ZeroIMedian) <- c("Zero-included Median")

  for(k in 1:nrow(RNAdataDF_R)) {
    shelf <- RNAdataDF_R[k,]
    row.names(shelf) <- c(1)
    shelf2 <- shelf[,(shelf[1,]>0)]
    ZeroXMedian[k,] <- median(as.numeric(shelf2))
    ZeroIMedian[k,] <- median(as.numeric(shelf))
  }


  DipOutputDF$ZeroXMedian <- as.numeric(ZeroXMedian)
  DipOutputDF$ZeroIMedian <- as.numeric(ZeroIMedian)
  DipOutputDF <<- dplyr::mutate(DipOutputDF, Divergence = ZeroXMedian-ZeroIMedian)

  return(DipOutputDF)

  if(samples==1){sublabels = "Sample A"}
  if(samples==2){sublabels = c("Sample A", "Sample B")}
  if(samples==3){sublabels = c("Sample A", "Sample B", "Sample C")}
  if(samples==4){sublabels = c("Sample A", "Sample B", "Sample C", "Sample D")}
  if(samples==5){sublabels = c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E")}
  if(samples==6){sublabels = c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E", "Sample F")}

  #plotting (executables)

  if(regions>=1){	#pick four from whole distribution
    ifelse(regions==1, DF.1 <- DipOutputDF, DF.1 <- DipOutputDF[DipOutputDF$Region==labels[1],])
    ifelse(nrow(DF.1)<samples, print("error: too few observations
                                     remaining in region 1 for declared sampling number to be selected"), print("Region 1 samples selected"))
    DF.1Samples <- DF.1[sample(nrow(DF.1), size = samples, replace = FALSE),]
    DF.1Names <- DF.1Samples[,"GeneID"]
    DF.1SampleGeneData <- data.frame(t(RNAdataDF_R[row.names(RNAdataDF_R) %in% DF.1Names, ]))
    SampleNames <- DF.1Names
    R1A <- paste(SampleNames[1], collapse = NULL)
    pp1 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1A)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1A) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    if(samples>=2){
      R1B <- paste(SampleNames[2], collapse = NULL)
      pp2 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1B)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::xlab(xx) + ggplot2::ggtitle(R1B) + ggplot2::ylab("Density")
    }
    if(samples>=3){
      R1C <- paste(SampleNames[3], collapse = NULL)
      pp3 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1C)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1C) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=4){
      R1D <- paste(SampleNames[4], collapse = NULL)
      pp4 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1D)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1D) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=5){
      R1E <- paste(SampleNames[5], collapse = NULL)
      pp5 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1E)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1E) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=6){
      R1F <- paste(SampleNames[6], collapse = NULL)
      pp6 <- ggplot2::ggplot(DF.1SampleGeneData, aes_string(x=R1F)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1F) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples==1){x1 <<- pp1}
    if(samples==2){x1 <<- gridExtra::grid.arrange(pp1, pp2, ncol=2)}
    if(samples==3){x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol=3)}
    if(samples==4){x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4,  ncol=2)}
    if(samples==5){x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol=2)}
    if(samples==6){x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol=2)}
  }

  if(regions>=2){
    DF.2 <- DipOutputDF[DipOutputDF$Region==labels[2],]
    ifelse(nrow(DF.2)<samples, print("error: too few observations
                                     remaining in region 2 for declared number of samples to be selected"), print("Region 2 samples selected"))
    DF.2Samples <- DF.2[sample(nrow(DF.2), size = samples, replace = FALSE),]
    DF.2Names <- DF.2Samples[,"GeneID"]
    DF.2SampleGeneData <- data.frame(t(RNAdataDF_R[row.names(RNAdataDF_R) %in% DF.2Names, ]))
    SampleNames2 <- DF.2Names
    R2A <- paste(SampleNames2[1], collapse = NULL)
    pp1 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2A)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1A) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    if(samples>=2){
      R2B <- paste(SampleNames2[2], collapse = NULL)
      pp2 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2B)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::xlab(xx) + ggplot2::ggtitle(R1B) + ggplot2::ylab("Density")
    }
    if(samples>=3){
      R2C <- paste(SampleNames2[3], collapse = NULL)
      pp3 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2C)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1C) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=4){
      R2D <- paste(SampleNames2[4], collapse = NULL)
      pp4 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2D)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1D) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=5){
      R2E <- paste(SampleNames2[5], collapse = NULL)
      pp5 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2E)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1E) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=6){
      R2F <- paste(SampleNames2[6], collapse = NULL)
      pp6 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R2F)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R1F) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples==1){x2 <<- pp1}
    if(samples==2){x2 <<- gridExtra::grid.arrange(pp1, pp2, ncol=2)}
    if(samples==3){x2 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol=3)}
    if(samples==4){x2 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4,  ncol=2)}
    if(samples==5){x2 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol=2)}
    if(samples==6){x2 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol=2)}
  }

  if(regions>=3){
    DF.3 <- DipOutputDF[DipOutputDF$Region==labels[3],]
    ifelse(nrow(DF.3)<samples, print("error: too few observations
                                     remaining in region 3 for declared number of samples to be selected"), print("Region 3 samples selected"))
    DF.3Samples <- DF.3[sample(nrow(DF.3), size = samples, replace = FALSE),]
    DF.3Names <- DF.3Samples[,"GeneID"]
    DF.3SampleGeneData <- data.frame(t(RNAdataDF_R[row.names(RNAdataDF_R) %in% DF.3Names, ]))
    SampleNames3 <- DF.3Names
    R3A <- paste(SampleNames3[1], collapse = NULL)
    pp1 <- ggplot2::ggplot(DF.3SampleGeneData, aes_string(x=R3A)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R3A) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    if(samples>=2){
      R3B <- paste(SampleNames3[2], collapse = NULL)
      pp2 <- ggplot2::ggplot(DF.3SampleGeneData, aes_string(x=R3B)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::xlab(xx) + ggplot2::ggtitle(R3B) + ggplot2::ylab("Density")
    }
    if(samples>=3){
      R3C <- paste(SampleNames3[3], collapse = NULL)
      pp3 <- ggplot2::ggplot(DF.3SampleGeneData, aes_string(x=R3C)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R3C) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=4){
      R3D <- paste(SampleNames3[4], collapse = NULL)
      pp4 <- ggplot2::ggplot(DF.3SampleGeneData, aes_string(x=R3D)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R3D) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=5){
      R3E <- paste(SampleNames3[5], collapse = NULL)
      pp5 <- ggplot2::ggplot(DF.3SampleGeneData, aes_string(x=R3E)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R3E) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=6){
      R3F <- paste(SampleNames3[6], collapse = NULL)
      pp6 <- ggplot2::ggplot(DF.2SampleGeneData, aes_string(x=R3F)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R3F) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples==1){x3 <<- pp1}
    if(samples==2){x3 <<- gridExtra::grid.arrange(pp1, pp2, ncol=2)}
    if(samples==3){x3 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol=3)}
    if(samples==4){x3 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4,  ncol=2)}
    if(samples==5){x3 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol=2)}
    if(samples==6){x3 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol=2)}
  }

  if(regions>=4){
    DF.4 <- DipOutputDF[DipOutputDF$Region==labels[4],]
    ifelse(nrow(DF.4)<samples, print("error: too few observations
                                     remaining in region 4 for declared number of samples to be selected"), print("Region 4 samples selected"))
    DF.4Samples <- DF.4[sample(nrow(DF.4), size = samples, replace = FALSE),]
    DF.4Names <- DF.4Samples[,"GeneID"]
    DF.4SampleGeneData <- data.frame(t(RNAdataDF_R[row.names(RNAdataDF_R) %in% DF.4Names, ]))
    SampleNames4 <- DF.4Names
    R4A <- paste(SampleNames4[1], collapse = NULL)
    pp1 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4A)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R4A) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    if(samples>=2){
      R4B <- paste(SampleNames4[2], collapse = NULL)
      pp2 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4B)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::xlab(xx) + ggplot2::ggtitle(R4B) + ggplot2::ylab("Density")
    }
    if(samples>=3){
      R4C <- paste(SampleNames4[3], collapse = NULL)
      pp3 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4C)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R4C) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=4){
      R4D <- paste(SampleNames4[4], collapse = NULL)
      pp4 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4D)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R4D) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=5){
      R4E <- paste(SampleNames4[5], collapse = NULL)
      pp5 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4E)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R4E) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=6){
      R4F <- paste(SampleNames4[6], collapse = NULL)
      pp6 <- ggplot2::ggplot(DF.4SampleGeneData, aes_string(x=R4F)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R4F) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples==1){x4 <<- pp1}
    if(samples==2){x4 <<- gridExtra::grid.arrange(pp1, pp2, ncol=2)}
    if(samples==3){x4 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol=3)}
    if(samples==4){x4 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4,  ncol=2)}
    if(samples==5){x4 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol=2)}
    if(samples==6){x4 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol=2)}
  }

  if(regions>=5){
    DF.5 <- DipOutputDF[DipOutputDF$Region==labels[5],]
    ifelse(nrow(DF.5)<samples, print("error: too few observations
                                     remaining in region 5 for declared number of samples to be selected"), print("Region 5 samples selected"))
    DF.5Samples <- DF.5[sample(nrow(DF.5), size = samples, replace = FALSE),]
    DF.5Names <- DF.5Samples[,"GeneID"]
    DF.5SampleGeneData <- data.frame(t(RNAdataDF_R[row.names(RNAdataDF_R) %in% DF.5Names, ]))
    SampleNames5 <- DF.5Names
    R5A <- paste(SampleNames5[1], collapse = NULL)
    pp1 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5A)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R5A) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    if(samples>=2){
      R2B <- paste(SampleNames5[2], collapse = NULL)
      pp2 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5B)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::xlab(xx) + ggplot2::ggtitle(R5B) + ggplot2::ylab("Density")
    }
    if(samples>=3){
      R5C <- paste(SampleNames5[3], collapse = NULL)
      pp3 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5C)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R5C) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=4){
      R2D <- paste(SampleNames5[4], collapse = NULL)
      pp4 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5D)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R5D) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=5){
      R5E <- paste(SampleNames5[5], collapse = NULL)
      pp5 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5E)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R5E) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples>=6){
      R5F <- paste(SampleNames5[6], collapse = NULL)
      pp6 <- ggplot2::ggplot(DF.5SampleGeneData, aes_string(x=R5F)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() + ggplot2::ggtitle(R5F) + ggplot2::ylab("Density") + ggplot2::xlab(xx)
    }
    if(samples==1){x5 <<- pp1}
    if(samples==2){x5 <<- gridExtra::grid.arrange(pp1, pp2, ncol=2)}
    if(samples==3){x5 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol=3)}
    if(samples==4){x5 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4,  ncol=2)}
    if(samples==5){x5 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol=2)}
    if(samples==6){x5 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol=2)}
  }

  return(x1)
  try(return(x2))
  try(return(x3))
  try(return(x4))
  try(return(x5))

}

plotSamplesByName <- function(nameList, RNAdata) {

  if(missing(ncol) && length(nameList)==any(c(2, 4, 6))) {
    ncol <- 2
  }

  if(missing(ncol) && length(nameList)==any(c(3, 5, 7))) {
    ncol <- 3
  }

  if(missing(ncol) && length(nameList)==1) {
    ncol <- 1
  }

  if(missing(ncol) && length(nameList)==8) {
    ncol <- 4
  }


  RNAdataDF <- data.frame(read.table(RNAdata))
  GeneData <- data.frame(t(RNAdataDF[row.nameList(RNAdataDF) %in% nameList, ]))
  R1 <- paste(nameList[1], collapse = NULL)
  try(R2 <- paste(nameList[2], collapse = NULL))
  try(R3 <- paste(nameList[3], collapse = NULL))
  try(R4 <- paste(nameList[4], collapse = NULL))
  try(R5 <- paste(nameList[5], collapse = NULL))
  try(R6 <- paste(nameList[6], collapse = NULL))
  try(R7 <- paste(nameList[7], collapse = NULL))
  try(R8 <- paste(nameList[8], collapse = NULL))
  pp1 <- ggplot2::ggplot(GeneData, aes_string(x=R1)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
    ggplot2::ggtitle(R1) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")

  if(!is.na(nameList[2])) {
    pp2 <- ggplot2::ggplot(GeneData, aes_string(x=R2)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R2) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[3])) {
    pp3 <- ggplot2::ggplot(GeneData, aes_string(x=R3)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R3) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[4])) {
    pp4 <- ggplot2::ggplot(GeneData, aes_string(x=R4)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R4) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[5])) {
    pp5 <- ggplot2::ggplot(GeneData, aes_string(x=R5)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R5) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[6])) {
    pp6 <- ggplot2::ggplot(GeneData, aes_string(x=R6)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R6) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[7])) {
    pp7 <- ggplot2::ggplot(GeneData, aes_string(x=R7)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R7) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(!is.na(nameList[8])) {
    pp8 <- ggplot2::ggplot(GeneData, aes_string(x=R8)) + ggplot2::geom_density(adjust = 1/5, color = "navy blue") + ggplot2::theme_gray() +
      ggplot2::ggtitle(R8) + ggplot2::ylab("Density") + ggplot2::xlab("Relative Expression")
  }

  if(length(nameList)==1){
    x1 <<- pp1
  }

  if(length(nameList)==2){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, ncol = ncol)
  }

  if(length(nameList)==3){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, ncol = ncol)
  }

  if(length(nameList)==4){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, ncol = ncol)
  }

  if(length(nameList)==5){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, ncol = ncol)
  }

  if(length(nameList)==6){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, ncol = ncol)
  }

  if(length(nameList)==7){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, pp7, ncol = ncol)
  }

  if(length(nameList)==8){
    x1 <<- gridExtra::grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, ncol = ncol)
  }
}

previewDipDistribution <- function(RNAdata, rawRNAdata, minimumCounts, barLine){

  RNAdataMat <- as.matrix(RNAdata)
  if(missing(barLine)){barLine <- 0}

  if(missing(rawRNAdata)){
    RNAdataMat <- as.matrix(RNAdata)
  } else{
    RNAdataMat <- as.matrix(RNAdata)
    RNAdataDF <- data.frame(RNAdataMat)
    rawRNAdata$max <- 0
    rawRNAdata$max <- apply(rawRNAdata, 1, max)
    RNArawbelow <- rawRNAdata[which(rawRNAdata$max<minimumCounts),]
    rr <- row.names(RNArawbelow)
    RNAdataDF_R <- RNAdataDF[!(row.names(RNAdataDF) %in% rr),]
    RNAdataMat <- as.matrix(RNAdataDF_R)
  }

  DipOutput <- matrix(nrow=nrow(RNAdataMat), ncol=2) #preparing a matrix to receive output in the function step

  for (k in 1:nrow(RNAdataMat)){

    e <- diptest::dip(RNAdataMat[k,], full.result = FALSE, min.is.0 = FALSE, debug = FALSE) #actual dip statistical test
    f <- diptest::dip.test(RNAdataMat[k,], simulate.p.value = FALSE, B = 2000)
    DipOutput[k,1] <- e #save output into pre-prepared matrix
    DipOutput[k,2] <- f$p_value
    ##consider extracting the rownames off RNAdata as a vector, which you can then reimpute here
  }

  DipOutputDF <- data.frame(DipOutput)
  colnames(DipOutputDF) <- c("Dip", "p-value")
  row.names(DipOutputDF) <- row.names(RNAdataMat)
  DipOutputDF <<- DipOutputDF

  x1 <<- ggplot2::ggplot(DipOutputDF, aes(x=Dip)) +
    ggplot2::geom_density() + ggplot2::geom_vline(xintercept=barLine, colour="#FF9999") +
    ggplot2::scale_x_continuous(breaks = round(seq(min(DipOutputDF$Dip), max(DipOutputDF$Dip), by = 0.01),2)) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::labs(title= "Bimodality among RNA expression patterns",x="Hartigan's Dip Score")

  x2 <<- ggplot2::ggplot(DipOutputDF, aes(x = p-value)) +
    ggplot2::geom_density() + ggplot2::geom_vline(xintercept = barLine, colour = "#FF9999") +
    ggplot2::scale_x_continuous(breaks = round(seq(min(DipOutputDF$p-value),
                                          max(DipOutputDF$p-value), by = 0.01), 2)) + ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::labs(title = "Bimodality among RNA expression patterns",
         x = "Hartigan's Dip Score")
  return(x1)
  return(x2)
}

#make the adjust into an argument
