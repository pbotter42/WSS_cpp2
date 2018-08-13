## Requried cpp functions
Rcpp::sourceCpp("calculateTraceLines.cpp")
Rcpp::sourceCpp("lw1.cpp")
Rcpp::sourceCpp("calculateposterior.cpp")
Rcpp::sourceCpp("calculateInformation.cpp")
##
## Required R packages
library(dplyr)
##
## Required R Functions
# aggregateScores4 #
# normal_distribution #
# USS_mtr
# uncollapsed_wEAP_cal #
# rpPatternExtract #
# fullTest #
##

### normal_distribution ###
normal_distribution <- function(nQuad,     
                                thetaMin,     
                                thetaMax) {     
  qPoints <- seq(from = thetaMin,     
                 to = thetaMax,     
                 by = (thetaMax-thetaMin)/(nQuad-1))     
  
  popDist <- exp(-(qPoints^2)/2)     
  popDist <- popDist/sum(popDist)     
  return(popDist)     
}

### uncollapsed_wEAP_cal ###
uncollapsed_wEAP_cal <- function(x,
                                 distSim,
                                 qPoints) {
  
  uncollapsedWeights <- as.matrix(x[,1])
  uncollapsedWithOutWeights <- as.matrix(x[,-1])
  
  uncollapsedPosterior <- calculatePosterior(summedTraceLines = uncollapsedWithOutWeights,
                                             dist = distSim)
  
  uncollapsedProbabilities <- rowSums(uncollapsedPosterior)
  
  # calculate wEAP... will try to calculate in arma later
  uncollapsed_wEAP <- numeric()
  for(i in 1:length(uncollapsedProbabilities)) {
    uncollapsed_wEAP[i] <- sum((uncollapsedPosterior[i,]*qPoints)/uncollapsedProbabilities[i])
  }
  
  return(uncollapsed_wEAP)
}

### rpPatternExtract ###
rpPatternExtract <- function(x) {
  out <- character(length = length(x))
  # looping through all responses in LW iter 1 through n
  for(i in 1:length(x)) {
    for(j in 1:length(x[[1]])) {
      if (x[[i]][j] == 0) {
        out[i] <- paste(out[i],"0",sep = "")
      } 
      else 
        out[i] <- paste(out[i],"1",sep = "")
    }
  }
  return(out)
}

### aggregateScores4 ###
aggregateScores4 <- function(x) {
  #df <- data.frame(nrow = 8,
  #                 ncol = 22)
  df <- data.frame()
  for(i in 1:length(x)) {
    df[i,1] <- unlist(x[[i]][2])
    for(qp in 1:length(unlist(x[[i]][1]))) {
      df[i, qp+1] <- unlist(x[[i]][1])[qp]
    }
  }
  
  df_aggregated_collapsed <- df %>% 
    dplyr::group_by(V1) %>% 
    dplyr::summarize_all(sum) %>% 
    as.data.frame()
  
  df_aggregated_uncollapsed <- df %>% 
    dplyr::group_by(V1) %>% 
    dplyr::mutate_all(sum) %>% 
    as.data.frame()  
  
  aggregateScoresOut <- list()
  aggregateScoresOut[["df_aggregated_collapsed"]] <- df_aggregated_collapsed
  aggregateScoresOut[["df_aggregated_uncollapsed"]] <- df_aggregated_uncollapsed
  aggregateScoresOut[["allScores"]] <- df[,1]
  aggregateScoresOut[["unsummedScores"]] <- df[,-1]  
  aggregateScoresOut[["allLik"]] <- df
  return(aggregateScoresOut)
}

### USS_mtr ###
USS_mtr <- function(theta,
                    thr,
                    sl) {
  
  wt <- rep(1, length(sl))
  nItems = length(thr)
  nQuad = length(theta)
  
  allTlines = calculateTraceLines(theta = theta,
                                  thr = thr,
                                  sl = sl,
                                  wt = wt,
                                  nItems = nItems,
                                  nQuad = nQuad)
  
  lwFinalIter = lw1(allTlines = allTlines,
                    nItems = nItems)[[1]]
  
  aggregateScoresOut <- aggregateScores4(lwFinalIter)
  summedScores <- aggregateScoresOut[["df_aggregated_collapsed"]]
  unsummedScores <- aggregateScoresOut[["unsummedScores"]]  
  allScores <- aggregateScoresOut[["allScores"]]
  
  # Numver of Quadrature points     
  nQuad <- 21     
  # theta range      
  thetaMin <- -5     
  thetaMax <- 5    
  
  qPoints <- seq(from = thetaMin,     
                 to = thetaMax,     
                 by = (thetaMax-thetaMin)/(nQuad-1))     
  
  distSim <- normal_distribution(nQuad,     
                                 thetaMin,     
                                 thetaMax)  
  
  weights <- as.matrix(summedScores[,1])
  withOutWeights <- as.matrix(summedScores[,-1])
  
  collapsedPosterior <- calculatePosterior(summedTraceLines = withOutWeights,
                                           dist = distSim)
  
  collapsedProbabilities <- rowSums(collapsedPosterior)
  
  # calculate wEAP... will try to calculate in arma later
  wEAP <- numeric()
  for(i in 1:length(collapsedProbabilities)) {
    wEAP[i] <- sum((collapsedPosterior[i,]*qPoints)/collapsedProbabilities[i])
  }
  # calculate SE and item marginal info
  SE <- numeric()
  itemMarginalInfo <- numeric()
  for(i in 1:length(collapsedProbabilities)) {
    SE[i] <- sqrt(sum(collapsedPosterior[i,]*(qPoints - wEAP[i])^2)/collapsedProbabilities[i])
    itemMarginalInfo[i] <- 1 - SE[i]^2
  }
  
  marginalTestReliability <- 1 - sum(collapsedProbabilities*SE)^2
  
  return(marginalTestReliability)
}

### fullTest ###
fullTest <- function(theta,
                     thr,
                     sl,
                     wt) {
  
  nItems = length(thr);
  nQuad = length(theta);
  
  allTlines = calculateTraceLines(theta = theta,
                                  thr = thr,
                                  sl = sl,
                                  wt = wt,
                                  nItems = nItems,
                                  nQuad = nQuad)
  
  lwFinalIter = lw1(allTlines = allTlines,
                    nItems = nItems)[[1]]
  
  # numeric matrix of only liklihoods
  likMat <- matrix(nrow = length(lwFinalIter),
                   ncol = length(theta))
  
  for(j in 1:length(lwFinalIter)) {
    likMat[j,] <- lwFinalIter[[j]][[1]]
  }
  
  #aggregateScoresOut <- aggregateScores2(lwFinalIter)
  #summedScores <- aggregateScoresOut[["summedScores"]]
  #allScores <- aggregateScoresOut[["allScores"]]
  
  # Numver of Quadrature points     
  nQuad <- 21     
  # theta range      
  thetaMin <- -5     
  thetaMax <- 5    
  
  qPoints <- seq(from = thetaMin,     
                 to = thetaMax,     
                 by = (thetaMax-thetaMin)/(nQuad-1))     
  
  distSim <- normal_distribution(nQuad,     
                                 thetaMin,     
                                 thetaMax)  
  
  #weights <- as.matrix(summedScores[,1])
  #withOutWeights <- as.matrix(summedScores[,-1])
  
  posterior <- calculatePosterior(summedTraceLines = likMat,
                                  dist = distSim)
  
  probabilities <- rowSums(posterior)
  
  # calculate wEAP... will try to calculate in arma later
  wEAP <- numeric()
  for(i in 1:length(probabilities)) {
    wEAP[i] <- sum((posterior[i,]*qPoints)/probabilities[i])
  }
  # calculate SE and item marginal info
  SE <- numeric()
  itemMarginalInfo <- numeric()
  for(i in 1:length(probabilities)) {
    SE[i] <- sqrt(sum(posterior[i,]*(qPoints - wEAP[i])^2)/probabilities[i])
    itemMarginalInfo[i] <- 1 - SE[i]^2
  }
  
  marginalTestReliability <- 1 - sum(probabilities*SE)^2
  
  # calculating Item information
  itemInfo = calculateInformation(traceLines = allTlines,
                                  slope = sl,
                                  qPoints = qPoints)
  
  imr <- rowSums(itemInfo) # calculating item marginal info
  
  
  return(marginalTestReliability)
}
