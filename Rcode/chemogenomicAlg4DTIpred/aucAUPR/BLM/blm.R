

blm <- function(Yfold = Yfold,
                Kd = Kd,
                Kt = Kt,
                cc = 1) {
  ### INPUT:
  ##
  
  
  ### OUTPUT:
  ##
  
  
  
  ## prediction by target. Voi moi Target thi du doan score cho toan bo drug.
  ## Do do, voi moi target, se co nDrug sample (duoc bieu dien bang kernel matrix Kd)
  ## number of target
  nDrug <- nrow(Yfold)
  nDrug
  nTarget <- ncol(Yfold)
  nTarget
  YpredByTarget <- matrix(NA, nrow = nDrug, ncol = nTarget)

  ## pre-computing
  trainK <- as.kernelMatrix(Kd)  
  
  for (i in 1:nTarget) {
    ###################### change to factor  #################
    curY <- factor(Yfold[, i])

    model <- ksvm(trainK, curY, kernel = "matrix", C = cc, cross = 0)

    testK <- as.kernelMatrix(Kd[, SVindex(model), drop = F])
    
    #cat(length(SVindex(model)), " ")
    
    YpredByTarget[, i] <- predict(model, testK, type = "decision")
  }
  
  ## prediction by drug. Voi moi Drug thi du doan score cho toan bo Target.
  ## Do do, voi moi Drug, se co nTarget sample (duoc bieu dien bang kernel matrix Kt)
  YfoldT <- t(Yfold)
  nTarget <- nrow(YfoldT)
  nDrug <- ncol(YfoldT)
  YpredByDrug <- matrix(NA, nrow = nTarget, ncol = nDrug)
  
  ## pre-computing
  trainK <- as.kernelMatrix(Kt)
  
  for (i in 1:nDrug) {
    ##############################################
    curY <- factor(YfoldT[, i]) ## must factor()
    ##############################################
    model <- ksvm(trainK, curY, kernel = "matrix", C = cc, cross = 0)
    testK <- as.kernelMatrix(Kt[, SVindex(model), drop = F])
    YpredByDrug[, i] <- predict(model, testK, type = "decision")
  }
  
  ## get maximum
  Ypred <- pmax(YpredByTarget, t(YpredByDrug))
  
  return(Ypred)
}








