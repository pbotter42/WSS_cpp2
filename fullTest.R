Rcpp::sourceCpp("calculateTraceLines.cpp")
Rcpp::sourceCpp("lwAlgorithmV1.cpp")
Rcpp::sourceCpp("calculateposterior.cpp")
Rcpp::sourceCpp("calculateInformation.cpp")
library(dplyr)
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
  
  lwFinalIter = lwAlgorithmV1(allTlines = allTlines,
                              nItems = nItems)
  
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
