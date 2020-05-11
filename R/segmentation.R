# file:   segmentation.R
# author: Gavin Ha, Ph.D.
#         Justin Rhoades
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   Oct 26, 2016
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param validInd PARAM_DESCRIPTION, Default: NULL
#' @param dataType PARAM_DESCRIPTION, Default: 'copy'
#' @param param PARAM_DESCRIPTION, Default: NULL
#' @param chrTrain PARAM_DESCRIPTION, Default: c(1:22)
#' @param maxiter PARAM_DESCRIPTION, Default: 50
#' @param estimateNormal PARAM_DESCRIPTION, Default: TRUE
#' @param estimatePloidy PARAM_DESCRIPTION, Default: TRUE
#' @param estimatePrecision PARAM_DESCRIPTION, Default: TRUE
#' @param estimateSubclone PARAM_DESCRIPTION, Default: TRUE
#' @param estimateTransition PARAM_DESCRIPTION, Default: TRUE
#' @param estimateInitDist PARAM_DESCRIPTION, Default: TRUE
#' @param logTransform PARAM_DESCRIPTION, Default: FALSE
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname HMMsegment

HMMsegment <- function(x, validInd = NULL, dataType = "copy", param = NULL,
    chrTrain = c(1:22), maxiter = 50, estimateNormal = TRUE, estimatePloidy = TRUE,
    estimatePrecision = TRUE, estimateSubclone = TRUE, estimateTransition = TRUE,
    estimateInitDist = TRUE, logTransform = FALSE, verbose = TRUE) {
  chr <- as.factor(seqnames(x[[1]]))
	# setup columns for multiple samples #
	dataMat <- as.matrix(as.data.frame(lapply(x, function(y) { mcols(x[[1]])[, dataType] })))

	# normalize by median and log data #
	if (logTransform){
    dataMat <- apply(dataMat, 2, function(x){ log(x / median(x, na.rm = TRUE)) })
	}else{
	  dataMat <- log(2^dataMat)
	}
	## update variable x with loge instead of log2
  for (i in 1:length(x)){
    mcols(x[[i]])[, dataType] <- dataMat[, i]
  }
  if (!is.null(chrTrain)) {
		chrInd <- chr %in% chrTrain
  }else{
  	chrInd <- !logical(length(chr))
  }
  if (!is.null(validInd)){
    chrInd <- chrInd & validInd
  }

	if (is.null(param)){
		param <- getDefaultParameters(dataMat[chrInd])
	}
	#if (param$n_0 == 0){
	#	param$n_0 <- .Machine$double.eps
	#}
	####### RUN EM ##########
  convergedParams <- runEM(dataMat, chr, chrInd, param, maxiter,
      verbose, estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
      estimateSubclone = estimateSubclone, estimatePrecision = estimatePrecision,
      estimateTransition = estimateTransition, estimateInitDist = estimateInitDist)
  # Calculate likelihood using converged params
 # S <- param$numberSamples
 # K <- length(param$ct)
 # KS <- K ^ S
 # py <- matrix(0, KS, nrow(dataMat))
 # iter <- convergedParams$iter
  # lambdasKS <- as.matrix(expand.grid(as.data.frame(convergedParams$lambda[, , iter])))
  # for (ks in 1:KS) {
  #   probs <- tdistPDF(dataMat, convergedParams$mus[ks, , iter], lambdasKS[ks, ], param$nu)
  #   py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  # }
  #
  viterbiResults <- runViterbi(convergedParams, chr)

  # setup columns for multiple samples #
  segs <- segmentData(x, validInd, viterbiResults$states, convergedParams)
  #output$segs <- processSegments(output$segs, chr, start(x), end(x), x$DataToUse)
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  #if (c(0) %in% param$ct){ #if state 0 HOMD is IN params#
  	#names <- c("HOMD", names)
  	# shift states to start at 2 (HETD)
    #tmp <- lapply(segs, function(x){ x$state <- x$state + 1; x})
    #viterbiResults$states <- as.numeric(viterbiResults$states) + 1
	#}
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  cnaList <- list()
  S <- length(x)
  for (s in 1:S){
    id <- names(x)[s]
    copyNumber <- param$jointCNstates[viterbiResults$state, s]
    subclone.status <- param$jointSCstatus[viterbiResults$state, s]
  	cnaList[[id]] <- data.frame(cbind(sample = as.character(id),
                  chr = as.character(seqnames(x[[s]])),
                  start = start(x[[s]]), end = end(x[[s]]),
                  copy.number = copyNumber,
                  event = names[copyNumber + 1],
                  logR = round(log2(exp(dataMat[,s])), digits = 4),
                  subclone.status = as.numeric(subclone.status)
  	))

    cnaList[[id]] <- transform(cnaList[[id]],
                              start = as.integer(as.character(start)),
                              end = as.integer(as.character(end)),
                              copy.number = as.numeric(copy.number),
                              logR = as.numeric(as.character(logR)),
                              subclone.status = as.numeric(subclone.status))

  	## order by chromosome ##
  	chrOrder <- unique(chr) #c(1:22,"X","Y")
  	cnaList[[id]] <- cnaList[[id]][order(match(cnaList[[id]][, "chr"],chrOrder)),]
  	## remove MT chr ##
    cnaList[[id]] <- cnaList[[id]][cnaList[[id]][,"chr"] %in% chrOrder, ]

    ## segment mean loge -> log2
    #segs[[s]]$median.logR <- log2(exp(segs[[s]]$median.logR))
    segs[[s]]$median <- log2(exp(segs[[s]]$median))
    ## add subclone status
    segs[[s]]$subclone.status <-  param$jointSCstatus[segs[[s]]$state, s]
  }
  convergedParams$segs <- segs
  return(list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param K PARAM_DESCRIPTION
#' @param e PARAM_DESCRIPTION
#' @param strength PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname getTransitionMatrix

getTransitionMatrix <- function(K, e, strength){
  A <- matrix(0, K, K)
  for (j in 1:K) {
    A[j, ] <- (1 - e[1]) / (K - 1)
    A[j, j] <- e[1]
  }
  A <- normalize(A)
  A_prior <- A
  dirPrior <- A * strength[1]
  return(list(A=A, dirPrior=dirPrior))
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param maxCN PARAM_DESCRIPTION, Default: 5
#' @param ct.sc PARAM_DESCRIPTION, Default: NULL
#' @param ploidy PARAM_DESCRIPTION, Default: 2
#' @param e PARAM_DESCRIPTION, Default: 0.9999999
#' @param e.sameState PARAM_DESCRIPTION, Default: 10
#' @param strength PARAM_DESCRIPTION, Default: 1e+07
#' @param includeHOMD PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname getDefaultParameters

getDefaultParameters <- function(x, maxCN = 5, ct.sc = NULL, ploidy = 2, e = 0.9999999, e.sameState = 10, strength = 10000000, includeHOMD = FALSE){
  if (includeHOMD){
    ct <- 0:maxCN
  }else{
    ct <- 1:maxCN
  }
	param <- list(
		strength = strength, e = e,
		ct = c(ct, ct.sc),
		ct.sc.status = c(rep(FALSE, length(ct)), rep(TRUE, length(ct.sc))),
		phi_0 = 2, alphaPhi = 4, betaPhi = 1.5,
		n_0 = 0.5, alphaN = 2, betaN = 2,
		sp_0 = 0.5, alphaSp = 2, betaSp = 2,
		lambda = as.matrix(rep(100, length(ct)+length(ct.sc)), ncol=1),
		nu = 2.1,
		kappa = rep(75, length(ct)),
		alphaLambda = 5
	)
	K <- length(param$ct)
  ## initialize hyperparameters for precision using observed data ##
	if (!is.null(dim(x))){ # multiple samples (columns)
    param$numberSamples <- ncol(x)
    #betaLambdaVal <- ((apply(x, 2, function(x){ sd(diff(x), na.rm=TRUE) }) / sqrt(length(param$ct))) ^ 2)
    betaLambdaVal <- ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	}else{ # only 1 sample
	  param$numberSamples <- 1
	  betaLambdaVal <- ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	}
	param$betaLambda <- matrix(betaLambdaVal, ncol = param$numberSamples, nrow = length(param$ct), byrow = TRUE)
  param$alphaLambda <- rep(param$alphaLambda, K)

	# increase prior precision for -1, 0, 1 copies at ploidy
	#param$lambda[param$ct %in% c(1,2,3)] <- 1000 # HETD, NEUT, GAIN
	#param$lambda[param$ct == 4] <- 100
	#param$lambda[which.max(param$ct)] <- 50 #highest CN
	#param$lambda[param$ct == 0] <- 1 #HOMD
	S <- param$numberSamples
	logR.var <- 1 / ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	if (!is.null(dim(x))){ # multiple samples (columns)
		param$lambda <- matrix(logR.var, nrow=K, ncol=S, byrow=T, dimnames=list(c(),colnames(x)))
	}else{ # only 1 sample
		#logR.var <- 1 / ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
    param$lambda <- matrix(logR.var, length(param$ct))
    param$lambda[param$ct %in% c(2)] <- logR.var
    param$lambda[param$ct %in% c(1,3)] <- logR.var
    param$lambda[param$ct >= 4] <- logR.var / 5
    param$lambda[param$ct == max(param$ct)] <- logR.var / 15
    param$lambda[param$ct.sc.status] <- logR.var / 10
  }
  # define joint copy number states #
  param$jointCNstates <- expand.grid(rep(list(param$ct), S))
  param$jointSCstatus <- expand.grid(rep(list(param$ct.sc.status), S))
  colnames(param$jointCNstates) <- paste0("Sample.", 1:param$numberSamples)
  colnames(param$jointSCstatus) <- paste0("Sample.", 1:param$numberSamples)

	# Initialize transition matrix to the prior
	txn <- getTransitionMatrix(K ^ S, e, strength)
  ## set higher transition probs for same CN states across samples ##
  # joint states where at least "tol" fraction of samples with the same CN state
	#apply(param$jointCNstates, 1, function(x){ sum(duplicated(as.numeric(x))) > 0 })
  cnStateDiff <- apply(param$jointCNstates, 1, function(x){ (abs(max(x) - min(x)))})
  if (e.sameState > 0 & S > 1){
		txn$A[, cnStateDiff == 0] <- txn$A[, cnStateDiff == 0] * e.sameState * K
		txn$A[, cnStateDiff >= 3] <- txn$A[, cnStateDiff >=3]  / e.sameState / K
	}
  for (i in 1:nrow(txn$A)){
    for (j in 1:ncol(txn$A)){
      if (i == j){
        txn$A[i, j] <- e
      }
    }
  }
  txn$A <- normalize(txn$A)
	param$A <- txn$A
	param$dirPrior <- txn$A * strength[1]
  param$A[, param$ct.sc.status] <- param$A[, param$ct.sc.status] / 10
  param$A <- normalize(param$A)
  param$dirPrior[, param$ct.sc.status] <- param$dirPrior[, param$ct.sc.status] / 10

  if (includeHOMD){
    K <- length(param$ct)
    param$A[1, 2:K] <- param$A[1, 2:K] * 1e-5; param$A[2:K, 1] <- param$A[2:K, 1] * 1e-5;
    param$A[1, 1] <- param$A[1, 1] * 1e-5
    param$A <- normalize(param$A); param$dirPrior <- param$A * param$strength
  }

  param$kappa <- rep(75, K ^ S)
  param$kappa[cnStateDiff == 0] <- param$kappa[cnStateDiff == 0] + 125
	param$kappa[cnStateDiff >=3] <- param$kappa[cnStateDiff >=3] - 50
	param$kappa[which(rowSums(param$jointCNstates==2) == S)] <- 800

  return(param)
}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dataGR PARAM_DESCRIPTION
#' @param validInd PARAM_DESCRIPTION
#' @param states PARAM_DESCRIPTION
#' @param convergedParams PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname segmentData

segmentData <- function(dataGR, validInd, states, convergedParams){
  if (sum(convergedParams$param$ct == 0) ==0){
  	includeHOMD <- FALSE
  }else{
  	includeHOMD <- TRUE
  }
  if (!includeHOMD){
    names <- c("HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }else{
    names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }
  states <- states[validInd]
  S <- length(dataGR)
  jointStates <- convergedParams$param$jointCNstates
  jointSCstatus <- convergedParams$param$jointSCstatus
  colNames <- c("seqnames", "start", "end", "copy")
  segList <- list()
  for (i in 1:S){
  	id <- names(dataGR)[i]
  	dataIn <- dataGR[[i]][validInd, ]
    rleResults <- t(sapply(runValue(seqnames(dataIn)), function(x){
      ind <- as.character(seqnames(dataIn)) == x
      r <- rle(states[ind])
    }))
    rleLengths <- unlist(rleResults[, "lengths"])
    rleValues <- unlist(rleResults[, "values"])
    sampleDF <- as.data.frame(dataIn)
    numSegs <- length(rleLengths)
    segs <- as.data.frame(matrix(NA, ncol = 7, nrow = numSegs,
                   dimnames = list(c(), c("chr", "start", "end", "state", "event", "median", "copy.number"))))
    prevInd <- 0
    for (j in 1:numSegs){
      start <- prevInd + 1
      end <- prevInd + rleLengths[j]
      segDF <- sampleDF[start:end, colNames]
      prevInd <- end
      numR <- nrow(segDF)
      segs[j, "chr"] <- as.character(segDF[1, "seqnames"])
      segs[j, "start"] <- segDF[1, "start"]
      segs[j, "state"] <- rleValues[j]
      segs[j, "copy.number"] <- jointStates[rleValues[j], i]
      if (segDF[1, "seqnames"] == segDF[numR, "seqnames"]){
        segs[j, "end"] <- segDF[numR, "end"]
        segs[j, "median"] <- round(median(segDF$copy, na.rm = TRUE), digits = 6)
        if (includeHOMD){
        	segs[j, "event"] <- names[segs[j, "copy.number"] + 1]
        }else{
        	segs[j, "event"] <- names[segs[j, "copy.number"]]
        }
      }else{ # segDF contains 2 different chromosomes
        print(j)
      }
    }
    segList[[id]] <- segs
  }
  return(segList)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param convergedParams PARAM_DESCRIPTION
#' @param chr PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname runViterbi

runViterbi <- function(convergedParams, chr){
  message("runViterbi: Segmenting and classifying")
  chrs <- levels(chr)
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  segs <- vector('list', length(chrs))
  py <- convergedParams$py
  N <- ncol(py)
  Z <- rep(0, N)
  convergeIter <- convergedParams$iter
  piG <- convergedParams$pi[, convergeIter]
  A <- convergedParams$A


  for(c in 1:length(chrsI)) {
    I <- chrsI[[c]]
    output <- .Call("viterbi", log(piG), log(A), log(py[, I]), PACKAGE = "HMMcopy")
    Z[I] <- output$path
    segs[[c]] <- output$seg
  }
  return(list(segs=segs, states=Z))
}

# Normalize a given array to sum to 1
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param A PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname normalize

normalize <- function(A) {
	vectorNormalize <- function(x){ x / (sum(x) + (sum(x) == 0)) }
	if (length(dim(A)) < 2){
  	M <- vectorNormalize(A)
  }else{
  	M <- t(apply(A, 1, vectorNormalize))
  }
  return(M);
}


# processSegments <- function(seg, chr, start, end, copy) {
#   segment <- data.frame()
#   chromosomes <- levels(chr)
#   for (i in 1:length(chromosomes)) {
#     seg_length = dim(seg[[i]])[1]
#     chr_name <- rep(chromosomes[i], seg_length)
#     chr_index <- which(chr == chromosomes[i])
#     chr_start <- start[chr_index][seg[[i]][, 1]]
#     chr_stop <- end[chr_index][seg[[i]][, 2]]
#     chr_state <- seg[[i]][, 3]
#     chr_median <- rep(0, seg_length)
#     for(j in 1:seg_length) {
#       chr_median[j] <-
#         median(na.rm = TRUE, log2(exp(copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])))
#     }
#     segment <- rbind(segment, cbind(chr = chr_name,
#       start = as.numeric(chr_start), end = chr_stop, state = chr_state,
#       median = chr_median))
#   }
#   segment <- transform(segment, start = as.numeric(as.character(start)),
#     end = as.numeric(as.character(end)), as.numeric(as.character(state)),
#     median = as.numeric(as.character(median)))
#   return(segment)
# }

#' Run IchorCNA
#'
#' @param WIG
#' @param NORMWIG
#' @param gcWig
#' @param mapWig
#' @param normalPanel
#' @param exons.bed
#' @param patientID
#' @param centromere
#' @param minMapScore
#' @param rmCentromereFlankLength
#' @param normal
#' @param scStates
#' @param coverage
#' @param lambda
#' @param lambdaScaleHyperParam
#' @param ploidy
#' @param maxCN
#' @param estimateNormal
#' @param estimateScPrevalence
#' @param estimatePloidy
#' @param maxFracCNASubclone
#' @param maxFracGenomeSubclone
#' @param minSegmentBins
#' @param altFracThreshold
#' @param chrNormalize
#' @param chrTrain
#' @param chrs
#' @param genomeBuild
#' @param genomeStyle
#' @param normalizeMaleX
#' @param minTumFracToCorrect
#' @param fracReadsInChrYForMale
#' @param includeHOMD
#' @param txnE
#' @param txnStrength
#' @param plotFileType
#' @param plotYLim
#' @param outDir
#' @param libdir
#'
#' @return
#' @export
#'
#' @examples
runIchorCNA <- function(WIG = NULL, NORMWIG = NULL, gcWig = NULL, mapWig = NULL, normalPanel = NULL,
                        exons.bed = NULL, patientID = "test", centromere = NULL, minMapScore = 0.9,
                        rmCentromereFlankLength = 1e+05, normal = 0.5, scStates = NULL, coverage = NULL,
                        lambda = NULL, lambdaScaleHyperParam = 3, ploidy = 2, maxCN = 7, estimateNormal = TRUE,
                        estimateScPrevalence = TRUE, estimatePloidy = TRUE, maxFracCNASubclone = 0.7,
                        maxFracGenomeSubclone = 0.5, minSegmentBins = 50, altFracThreshold = 0.05,
                        chrNormalize = c(1:22), chrTrain = c(1:22), chrs = c(1:22, "X"), genomeBuild = "hg19",
                        genomeStyle = "NCBI", normalizeMaleX = TRUE, minTumFracToCorrect = 0.1,
                        fracReadsInChrYForMale = 0.001, includeHOMD = FALSE, txnE = 0.9999999, txnStrength = 1e+07,
                        plotFileType = "pdf", plotYLim = c(-2, 2), outDir = "test", libdir = NULL) {

  if (!file.exists(outdir)){
    dir.create(outDir)
  }

  # browser()

  # set params------------------------------

  # ------------------------------
  tumour_file <- WIG
  normal_file <- NORMWIG
  gcWig <- gcWig
  mapWig <- mapWig
  normal_panel <- normalPanel
  exons.bed <- exons.bed # "0" if none specified
  centromere <- centromere
  minMapScore <- minMapScore
  flankLength <- rmCentromereFlankLength
  normal <- normal
  scStates <- scStates
  lambda <- lambda
  # normal <- eval(parse(text = normal))
  # scStates <- eval(parse(text = scStates))
  # lambda <- eval(parse(text = lambda))
  lambdaScaleHyperParam <- lambdaScaleHyperParam
  estimateNormal <- estimateNormal
  estimatePloidy <- estimatePloidy
  estimateScPrevalence <- estimateScPrevalence
  maxFracCNASubclone <- maxFracCNASubclone
  maxFracGenomeSubclone <- maxFracGenomeSubclone
  minSegmentBins <- minSegmentBins
  altFracThreshold <- altFracThreshold
  ploidy <- ploidy
  # ploidy <- eval(parse(text = ploidy))
  coverage <- coverage
  maxCN <- maxCN
  txnE <- txnE
  txnStrength <- txnStrength
  normalizeMaleX <- as.logical(normalizeMaleX)
  includeHOMD <- as.logical(includeHOMD)
  minTumFracToCorrect <- minTumFracToCorrect
  fracReadsInChrYForMale <- fracReadsInChrYForMale
  chrXMedianForMale <- -0.1
  outDir <- outDir
  libdir <- libdir
  plotFileType <- plotFileType
  plotYLim <- plotYLim
  # plotYLim <- eval(parse(text=plotYLim))
  gender <- NULL
  outImage <- paste0(outDir, "/", patientID, ".RData")
  genomeBuild <- genomeBuild
  genomeStyle <- genomeStyle
  chrs <- chrs
  # chrs <- as.character(eval(parse(text = chrs)))
  chrTrain <- as.character(chrTrain)
  # chrTrain <- as.character(eval(parse(text=chrTrain)));
  chrNormalize <- as.character(chrNormalize)
  # chrNormalize <- as.character(eval(parse(text=chrNormalize)));
  seqlevelsStyle(chrs) <- genomeStyle
  seqlevelsStyle(chrNormalize) <- genomeStyle
  seqlevelsStyle(chrTrain) <- genomeStyle


  ## load seqinfo
  # browser()
  seqinfo <- getSeqInfo(genomeBuild, genomeStyle)

  if (substr(tumour_file, nchar(tumour_file) - 2, nchar(tumour_file)) == "wig") {
    wigFiles <- data.frame(cbind(patientID, tumour_file))
  } else {
    wigFiles <- read.delim(tumour_file, header = F, as.is = T)
  }

  ## FILTER BY EXONS IF PROVIDED ##
  ## add gc and map to GRanges object ##
  if (is.null(exons.bed) || exons.bed == "None" || exons.bed == "NULL") {
    targetedSequences <- NULL
  } else {
    targetedSequences <- read.delim(exons.bed, header = T, sep = "\t")
  }

  ## load PoN
  if (is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL") {
    normal_panel <- NULL
  }

  if (is.null(centromere) || centromere == "None" || centromere == "NULL") { # no centromere file provided
    centromere <- system.file("extdata", "GRCh37.p13_centromere_UCSC-gapTable.txt",
      package = "ichorCNA"
    )
  }

  # ------------------------------
  centromere <- read.delim(centromere, header = T, stringsAsFactors = F, sep = "\t")
  save.image(outImage)
  ## LOAD IN WIG FILES ##
  numSamples <- nrow(wigFiles)

  tumour_copy <- list()
  for (i in 1:numSamples) {
    id <- wigFiles[i, 1]
    ## create output directories for each sample ##
    dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
    ### LOAD TUMOUR AND NORMAL FILES ###
    message("Loading tumour file:", wigFiles[i, 1])
    tumour_reads <- wigToGRanges(wigFiles[i, 2])

    ## LOAD GC/MAP WIG FILES ###
    # find the bin size and load corresponding wig files #
    binSize <- as.data.frame(tumour_reads[1, ])$width
    message("Reading GC and mappability files")
    if (is.null(gcWig) || gcWig == "None" || gcWig == "NULL") {
      stop("GC wig file is required")
    }
    gc <- wigToGRanges(gcWig)
    if (is.null(mapWig) || mapWig == "None" || mapWig == "NULL") {
      message("No mappability wig file input, excluding from correction")
      map <- NULL
    } else {
      map <- wigToGRanges(mapWig)
    }
    message("Correcting Tumour")

    counts <- loadReadCountsFromWig(tumour_reads,
      chrs = chrs, gc = gc, map = map,
      centromere = centromere, flankLength = flankLength,
      targetedSequences = targetedSequences, chrXMedianForMale = chrXMedianForMale,
      genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
      chrNormalize = chrNormalize, mapScoreThres = minMapScore
    )
    tumour_copy[[id]] <- counts$counts # as(counts$counts, "GRanges")
    gender <- counts$gender
    ## load in normal file if provided
    if (!is.null(normal_file) && normal_file != "None" && normal_file != "NULL") {
      message("Loading normal file:", normal_file)
      normal_reads <- wigToGRanges(normal_file)
      message("Correcting Normal")
      counts <- loadReadCountsFromWig(normal_reads,
        chrs = chrs, gc = gc, map = map,
        centromere = centromere, flankLength = flankLength, targetedSequences = targetedSequences,
        genomeStyle = genomeStyle, chrNormalize = chrNormalize, mapScoreThres = minMapScore
      )
      normal_copy <- counts$counts # as(counts$counts, "GRanges")
      gender.normal <- counts$gender
    } else {
      normal_copy <- NULL
    }

    ### DETERMINE GENDER ###
    ## if normal file not given, use chrY, else use chrX
    message("Determining gender...", appendLF = FALSE)
    gender.mismatch <- FALSE
    if (!is.null(normal_copy)) {
      if (gender$gender != gender.normal$gender) { # use tumour # use normal if given
        # check if normal is same gender as tumour
        gender.mismatch <- TRUE
      }
    }
    message("Gender ", gender$gender)

    ## NORMALIZE GENOME-WIDE BY MATCHED NORMAL OR NORMAL PANEL (MEDIAN) ##
    tumour_copy[[id]] <- normalizeByPanelOrMatchedNormal(tumour_copy[[id]],
      chrs = chrs,
      normal_panel = normal_panel, normal_copy = normal_copy,
      gender = gender$gender, normalizeMaleX = normalizeMaleX
    )

    ### OUTPUT FILE ###
    ### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
    outMat <- as.data.frame(tumour_copy[[id]])
    # outMat <- outMat[,c(1,2,3,12)]
    outMat <- outMat[, c("seqnames", "start", "end", "copy")]
    colnames(outMat) <- c("chr", "start", "end", "log2_TNratio_corrected")
    outFile <- paste0(outDir, "/", id, ".correctedDepth.txt")
    message(paste("Outputting to:", outFile))
    write.table(outMat, file = outFile, row.names = F, col.names = T, quote = F, sep = "\t")
  } ## end of for each sample

  chrInd <- as.character(seqnames(tumour_copy[[1]])) %in% chrTrain
  ## get positions that are valid
  valid <- tumour_copy[[1]]$valid
  if (length(tumour_copy) >= 2) {
    for (i in 2:length(tumour_copy)) {
      valid <- valid & tumour_copy[[i]]$valid
    }
  }
  save.image(outImage)

  ### RUN HMM ###
  ## store the results for different normal and ploidy solutions ##
  ptmTotalSolutions <- proc.time() # start total timer
  results <- list()
  loglik <- as.data.frame(matrix(NA,
    nrow = length(normal) * length(ploidy), ncol = 7,
    dimnames = list(c(), c(
      "init", "n_est", "phi_est", "BIC",
      "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"
    ))
  ))
  counter <- 1
  compNames <- rep(NA, nrow(loglik))
  mainName <- rep(NA, length(normal) * length(ploidy))
  #### restart for purity and ploidy values ####
  for (n in normal) {
    for (p in ploidy) {
      if (n == 0.95 & p != 2) {
        next
      }
      logR <- as.data.frame(lapply(tumour_copy, function(x) {
        x$copy
      })) # NEED TO EXCLUDE CHR X #
      param <- getDefaultParameters(logR[valid & chrInd, , drop = F],
        maxCN = maxCN, includeHOMD = includeHOMD,
        ct.sc = scStates, ploidy = floor(p), e = txnE, e.same = 50, strength = txnStrength
      )
      param$phi_0 <- rep(p, numSamples)
      param$n_0 <- rep(n, numSamples)

      ############################################
      ######## CUSTOM PARAMETER SETTINGS #########
      ############################################
      # 0.1x cfDNA #
      if (is.null(lambda)) {
        logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct)))^2)
        param$lambda <- rep(logR.var, length(param$ct))
        param$lambda[param$ct %in% c(2)] <- logR.var
        param$lambda[param$ct %in% c(1, 3)] <- logR.var
        param$lambda[param$ct >= 4] <- logR.var / 5
        param$lambda[param$ct == max(param$ct)] <- logR.var / 15
        param$lambda[param$ct.sc.status] <- logR.var / 10
      } else {
        param$lambda[param$ct %in% c(2)] <- lambda[2]
        param$lambda[param$ct %in% c(1)] <- lambda[1]
        param$lambda[param$ct %in% c(3)] <- lambda[3]
        param$lambda[param$ct >= 4] <- lambda[4]
        param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
        param$lambda[param$ct.sc.status] <- lambda[2] / 10
      }
      param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))
      # 1x bulk tumors #
      # param$lambda[param$ct %in% c(2)] <- 2000
      # param$lambda[param$ct %in% c(1)] <- 1750
      # param$lambda[param$ct %in% c(3)] <- 1750
      # param$lambda[param$ct >= 4] <- 1500
      # param$lambda[param$ct == max(param$ct)] <- 1000 / 25
      # param$lambda[param$ct.sc.status] <- 1000 / 75
      # param$alphaLambda[param$ct.sc.status] <- 4
      # param$alphaLambda[param$ct %in% c(1,3)] <- 5
      # param$alphaLambda[param$ct %in% c(2)] <- 5
      # param$alphaLambda[param$ct == max(param$ct)] <- 4

      #############################################
      ################ RUN HMM ####################
      #############################################
      hmmResults.cor <- HMMsegment(tumour_copy, valid,
        dataType = "copy",
        param = param, chrTrain = chrTrain, maxiter = 50,
        estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
        estimateSubclone = estimateScPrevalence, verbose = TRUE
      )

      for (s in 1:numSamples) {
        iter <- hmmResults.cor$results$iter
        id <- names(hmmResults.cor$cna)[s]

        ## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
        ## check if there is an altered segment that has at least a minimum # of bins
        segsS <- hmmResults.cor$results$segs[[s]]
        segsS <- segsS[segsS$chr %in% chrTrain, ]
        segAltInd <- which(segsS$event != "NEUT")
        maxBinLength <- -Inf
        if (sum(segAltInd) > 0) {
          maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
          maxSegRD <- GRanges(
            seqnames = segsS$chr[segAltInd[maxInd]],
            ranges = IRanges(start = segsS$start[segAltInd[maxInd]], end = segsS$end[segAltInd[maxInd]])
          )
          hits <- findOverlaps(query = maxSegRD, subject = tumour_copy[[s]][valid, ])
          maxBinLength <- length(subjectHits(hits))
        }
        ## check if there are proportion of total bins altered
        # if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
        cnaS <- hmmResults.cor$cna[[s]]
        altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
        altFrac <- sum(!altInd, na.rm = TRUE) / length(altInd)
        if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)) {
          hmmResults.cor$results$n[s, iter] <- 1.0
        }

        # correct integer copy number based on estimated purity and ploidy
        correctedResults <- correctIntegerCN(
          cn = hmmResults.cor$cna[[s]],
          segs = hmmResults.cor$results$segs[[s]],
          purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
          cellPrev = 1 - hmmResults.cor$results$sp[s, iter],
          maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = minTumFracToCorrect,
          gender = gender$gender, chrs = chrs, correctHOMD = includeHOMD
        )
        hmmResults.cor$results$segs[[s]] <- correctedResults$segs
        hmmResults.cor$cna[[s]] <- correctedResults$cn

        ## plot solution ##
        outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p)
        mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
        plotGWSolution(hmmResults.cor,
          s = s, outPlotFile = outPlotFile, plotFileType = plotFileType,
          logR.column = "logR", call.column = "Corrected_Call",
          plotYLim = plotYLim, estimateScPrevalence = estimateScPrevalence, seqinfo = seqinfo, main = mainName[counter]
        )
      }
      iter <- hmmResults.cor$results$iter
      results[[counter]] <- hmmResults.cor
      loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
      subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x) {
        sum(x$subclone.status)
      }))
      fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x) {
        nrow(x)
      }))
      fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x) {
        sum(x$copy.number != 2)
      }))
      fracAltSub <- lapply(fracAltSub, function(x) {
        if (is.na(x)) {
          0
        } else {
          x
        }
      })
      loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits = 2), collapse = ",")
      loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits = 2), collapse = ",")
      loglik[counter, "init"] <- paste0("n", n, "-p", p)
      loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
      loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

      counter <- counter + 1
    }
  }
  ## get total time for all solutions ##
  elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
  message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

  ### SAVE R IMAGE ###
  save.image(outImage)
  # save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

  ### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
  loglik <- loglik[!is.na(loglik$init), ]
  if (estimateScPrevalence) { ## sort but excluding solutions with too large % subclonal
    fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone &
      loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
    if (length(fracInd) > 0) { ## if there is a solution satisfying % subclonal
      ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing = TRUE)]
    } else { # otherwise just take largest likelihood
      ind <- order(as.numeric(loglik[, "loglik"]), decreasing = TRUE)
    }
  } else { # sort by likelihood only
    ind <- order(as.numeric(loglik[, "loglik"]), decreasing = TRUE)
  }

  # new loop by order of solutions (ind)
  outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_all_sols")
  for (i in 1:length(ind)) {
    hmmResults.cor <- results[[ind[i]]]
    turnDevOff <- FALSE
    turnDevOn <- FALSE
    if (i == 1) {
      turnDevOn <- TRUE
    }
    if (i == length(ind)) {
      turnDevOff <- TRUE
    }
    plotGWSolution(hmmResults.cor,
      s = s, outPlotFile = outPlotFile, plotFileType = "pdf",
      logR.column = "logR", call.column = "Corrected_Call",
      plotYLim = plotYLim, estimateScPrevalence = estimateScPrevalence,
      seqinfo = seqinfo,
      turnDevOn = turnDevOn, turnDevOff = turnDevOff, main = mainName[ind[i]]
    )
  }

  hmmResults.cor <- results[[ind[1]]]
  hmmResults.cor$results$loglik <- as.data.frame(loglik)
  hmmResults.cor$results$gender <- gender$gender
  hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
  hmmResults.cor$results$chrXMedian <- gender$chrXMedian
  hmmResults.cor$results$coverage <- coverage

  outputHMM(
    cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs,
    results = hmmResults.cor$results, patientID = patientID, outDir = outDir
  )
  outFile <- paste0(outDir, "/", patientID, ".params.txt")
  print(hmmResults.cor)
  outputParametersToFile(hmmResults.cor, file = outFile)

  ## plot solutions for all samples
  plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir,
    numSamples = numSamples,
    logR.column = "logR", call.column = "Corrected_Call",
    plotFileType = plotFileType, plotYLim = plotYLim, seqinfo = seqinfo,
    estimateScPrevalence = estimateScPrevalence, maxCN = maxCN
  )

  return(hmmResults.cor)

}
