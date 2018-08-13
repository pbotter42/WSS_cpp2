simFun1 <- function(x) {
  out <- list()
  for(i in 1:length(x)) {
    out[[i]] <- list()
    out[[i]][["threshold"]] <- x[[i]][["threshold"]]
    out[[i]][["slope"]] <- x[[i]][["slope"]]    
    out[[i]][["weights"]] <- x[[i]][["weights"]] 
    out[[i]][["simOutput"]] <- wssR(theta = theta,
                                    thr = x[[i]][["threshold"]],
                                    sl = x[[i]][["slope"]],
                                    wt = x[[i]][["weights"]])
    
    # extractiong the Weighting Scheme info names
    `Weighting Scheme info names` <- names(out[[i]][["simOutput"]]$`Weighting Scheme info`)
    
    # temp vectors containing 
    thrNames <- character()
    slopeNames <- character()
    weightNames <- character()
    imrNames <- character()
    for(j in 1:length(x[[i]][["slope"]])) {
      thrNames[j] <- paste("b ",j, sep = "")
      slopeNames[j] <- paste("a ",j, sep = "")
      weightNames[j] <- paste("weight ",j, sep = "")
      imrNames[j] <- paste("imr ",j, sep = "")      
    }
    
    # names for output for brian
    testLevelOutNames <- c(`Weighting Scheme info names`,
                           thrNames,
                           slopeNames,
                           weightNames,
                           imrNames)
    
    # output for brian
    testLevelOut <- c(as.numeric(out[[i]][["simOutput"]]$`Weighting Scheme info`),
                      x[[i]][["threshold"]],
                      x[[i]][["slope"]],
                      x[[i]][["weights"]],
                      as.numeric(out[[i]][["simOutput"]][["outList"]][["imr"]]))
    names(testLevelOut) <- testLevelOutNames
    out[[i]][["testLevelOut"]] <- testLevelOut
    
    # rp data bind 
    rpOut <- c(out[[i]][["simOutput"]][["outList"]][["allScores"]],
               out[[i]][["simOutput"]][["outList"]][["uncollapsed_wEAP"]],
               out[[i]][["simOutput"]][["outList"]][["rpEAP"]])
    out[[i]][["rpOut"]] <- rpOut
    
    # rp names
    wss_Names <- character()
    wEAP_Names <- character()
    rpEAP_Names <- character()
    for(j in 1:length(out[[i]][["simOutput"]][["outList"]][["allScores"]])) {
      wss_Names[j] <- paste("wss ",j, sep = "")
      wEAP_Names[j] <- paste("wEAP ",j, sep = "")
      rpEAP_Names[j] <- paste("rpEAP ",j, sep = "")      
    }   
    
    rpNames <- c(wss_Names,
                 wEAP_Names,
                 rpEAP_Names)
    out[[i]][["rpNames"]] <- rpNames
    ### END INITIAL FOR LOOP
  }
  
  outMatrix <- matrix(nrow = length(x),
                      ncol = length(out[[1]][["testLevelOut"]]))
  # getting all testLevelOuts from the main list and putting it into a mat
  for(i in 1:length(x)) {
    outMatrix[i,] <- as.numeric(out[[i]][["testLevelOut"]])
  }
  outDf <- as.data.frame(outMatrix)
  names(outDf) <- names(out[[1]][["testLevelOut"]])
  
  infoLoc <- character()
  for(i in 1:length(x)) {
    infoLoc[i] <- x[[i]][["infoLoc"]]
  }
  
  # createing matrix that contains rp info. so wss, wEAP, and rpEAP
  rpMatrix <- matrix(nrow = length(x),
                     ncol = length(out[[1]][["rpOut"]]))
  # puting rp info into matrix
  for(i in 1:length(x)) {
    rpMatrix[i,] <- out[[i]][["rpOut"]]
  }
  rpDf <- as.data.frame(rpMatrix)
  names(rpDf) <- out[[1]][["rpNames"]]
  
  outDf2 <- cbind(infoLoc,
                  outDf,
                  rpDf)
  out2 <- list()
  out2[[1]] <- out
  out2[[2]] <- outDf2
  return(out2)
}