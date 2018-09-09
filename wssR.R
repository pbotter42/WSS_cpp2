## Requried cpp functions
Rcpp::sourceCpp("calculateTraceLines.cpp")
Rcpp::sourceCpp("lw1.cpp")
Rcpp::sourceCpp("calculateposterior.cpp")
Rcpp::sourceCpp("calculateInformation.cpp")
Rcpp::sourceCpp("rpLikelihoodCalculate.cpp")
##
## Required R packages
library(dplyr)
##
## Required R Functions
# aggregateScores4
# normal_distribution
# USS_mtr
# uncollapsed_wEAP_cal
# rpPatternExtract
# imrCal
##
wssR <- function(theta,
                 thr,
                 sl,
                 wt,
                 type = c("LW"),
                 nSim
) {
  
  nItems <- length(thr)
  nQuad <- length(theta)
  
  allTlines <- calculateTraceLines(theta = theta,
                                   thr = thr,
                                   sl = sl,
                                   wt = wt,
                                   nItems = nItems,
                                   nQuad = nQuad)
  
  if(type == "LW") {  
    
    lwOut <- lw1(allTlines = allTlines,
                 nItems = nItems)
    lwFinalIter <- lwOut[[1]]
    rpRaw <- lwOut[[2]][[length(lwOut[[2]])]] # extract the last iteration of the raw rp pattern. so includes non 1s
    rpPattern <- rpPatternExtract(rpRaw, type = "LW") # convert raw pattern to 1s and 0s	  
    
  } else if(type == "rpSim") {
    simDat <- mirt::simdata(a = sl,
                            d = thr,
                            N = nSim,
                            itemtype = "2PL")
    lwFinalIter <- calculateRpLik(allTlines = allTlines,
                         testeeRp = simDat,
                         nItems = nItems)
    rpPattern <- rpPatternExtract(simDat, type = "rpSim") # convert raw pattern to 1s and 0s	  
  }
  
  ### Above I am tring to keep everything that differs based on "type" ###
  
  aggregateScoresOut <- aggregateScores4(lwFinalIter)
  df_aggregated_collapsed <- aggregateScoresOut[["df_aggregated_collapsed"]]
  df_aggregated_uncollapsed <- aggregateScoresOut[["df_aggregated_uncollapsed"]]
  unsummedScores <- aggregateScoresOut[["unsummedScores"]]  
  allScores <- aggregateScoresOut[["allScores"]]
  allLik <- aggregateScoresOut[["allLik"]]
  
  # Numver of Quadrature points     
  #nQuad <- 21     
  # theta range      
  thetaMin <- -5    
  thetaMax <- 5    
  
  qPoints <- theta     
  
  distSim <- normal_distribution(nQuad,     
                                 thetaMin,     
                                 thetaMax)  
  
  weights <- as.matrix(df_aggregated_collapsed[,1])
  withOutWeights <- as.matrix(df_aggregated_collapsed[,-1])
  
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
  
  # calculating Item information
  itemInfo = calculateInformation(traceLines = allTlines,
                                  slope = sl,
                                  qPoints = qPoints)
  
  imr <- imrCal(itemInfo, distSim) # calculating item marginal info
  #imr_norm <- rowSums(itemInfo) # calculating item marginal info  
  fullMarginalTestRelibility <- fullTest(theta = theta,
                                         thr = thr,
                                         sl = sl,
                                         wt = weights)
  
  uncollapsedPosterior <- calculatePosterior(summedTraceLines = as.matrix(unsummedScores),
                                             dist = distSim)
  
  uncollapsedProbabilities <- rowSums(uncollapsedPosterior)
  
  # calculate wEAP... will try to calculate in arma later
  rpEAP <- numeric()
  for(i in 1:length(uncollapsedProbabilities)) {
    rpEAP[i] <- sum((uncollapsedPosterior[i,]*qPoints)/uncollapsedProbabilities[i])
  }
  
  uncollapsed_wEAP <- uncollapsed_wEAP_cal(df_aggregated_uncollapsed,
                                           distSim,
                                           qPoints)
  
  outList <- list()
  outList[["allTlines"]] <- allTlines
  outList[["df_aggregated_collapsed"]] <- df_aggregated_collapsed
  outList[["df_aggregated_uncollapsed"]] <- df_aggregated_uncollapsed
  
  outList[["unsummedScores"]] <- unsummedScores 
  outList[["allLik"]] <- allLik
  outList[["allScores"]] <- allScores
  outList[["qPoints"]] <- qPoints
  outList[["distSim"]] <- distSim
  outList[["collapsedPosterior"]] <- collapsedPosterior
  outList[["wEAP"]] <- wEAP
  outList[["rpEAP"]] <- rpEAP
  outList[["uncollapsed_wEAP"]] <- uncollapsed_wEAP
  outList[["allScores"]] <- allScores
  outList[["marginalTestReliability"]] <- marginalTestReliability
  outList[["lwFinalIter"]] <- lwFinalIter
  outList[["imr"]] <- imr
  outList[["allScores"]] <- allScores
  outList[["weights"]] <- weights
  outList[["itemInfo"]] <- itemInfo
  outList[["fullMarginalTestRelibility"]] <- fullMarginalTestRelibility
  outList[["marginalTestRelibilityDif"]] <- fullMarginalTestRelibility - marginalTestReliability
  #outList[["rpRaw"]] <- rpRaw
  outList[["rpPattern"]] <- rpPattern  
  #outList[["partialTest"]] <- partialTest
  out <- list()
  ## weighting scheme info
  `Weighting Scheme info` <- data.frame()
  `Weighting Scheme info`[1,"RP mtr"] <- fullMarginalTestRelibility
  `Weighting Scheme info`[1,"USS mtr"] <- USS_mtr(theta = theta,
                                                  thr = thr,
                                                  sl = sl)  
  `Weighting Scheme info`[1,"WSS mtr"] <- marginalTestReliability
  `Weighting Scheme info`[1,"mtr diff"] <- fullMarginalTestRelibility - marginalTestReliability
  `Weighting Scheme info`[1,"# of unique WSS"] <- length(unique(allScores))
  #  `Weighting Scheme info`[1,"sp_weights_imr"] <- cor(wt, imr, method = "spearman")
  `Weighting Scheme info`[1,"sp_wss_wEAP"] <- cor(weights, wEAP, method = "spearman") 
  `Weighting Scheme info`[1,"cor_wss_wEAP"] <- cor(weights, wEAP, method = "pearson") 
  `Weighting Scheme info`[1,"sp_wss_rpEAP"] <- cor(allScores, rpEAP, method = "spearman") 
  `Weighting Scheme info`[1,"cor_wss_rpEAP"] <- cor(allScores, rpEAP, method = "pearson")  
  `Weighting Scheme info`[1,"sp_wEAP_rpEAP"] <- cor(uncollapsed_wEAP, rpEAP, method = "spearman") 
  `Weighting Scheme info`[1,"cor_wEAP_rpEAP"] <- cor(uncollapsed_wEAP, rpEAP, method = "pearson")  
  `Weighting Scheme info`[1,"sp_weight_imr"] <- cor(wt, imr, method = "spearman") 
  `Weighting Scheme info`[1,"cor_weight_imr"] <- cor(wt, imr, method = "pearson") 
  `Weighting Scheme info`[1,"sp_weight_a"] <- cor(wt, sl, method = "spearman") 
  `Weighting Scheme info`[1,"cor_weight_a"] <- cor(wt, sl, method = "pearson") 
  # Item params
  `Item Parameters` <- data.frame(WSS = weights,
                                  collapsedProbabilities = collapsedProbabilities,
                                  wEAP = wEAP,
                                  SE = SE,
                                  `item marginal information` = itemMarginalInfo)
  
  out[["Weighting Scheme info"]] <- `Weighting Scheme info`
  out[["Item Parameters"]] <- `Item Parameters`
  out[["outList"]] <- outList
  return(lwFinalIter)
}


# To Do
# 1. in the LW version I have the number of numque RP.. i may have to do this differently in the rpSim version
