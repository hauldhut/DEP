setwd("d:/Setup/BioInformatic/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM_Copy/Code")
#setwd("C:/Users/Mrs Giang/OneDrive - Hanoi University of Science and Technology/R/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM/")
rm(list = ls())
library(caret)# 
library(Metrics)#use for rmse calculation
library(mltools)# use for auc calculation
library("ggpubr")#use for correlative caculation
#######################################
# Load packages and list of files

## load required packages
pkgs <- c(
  "matrixcalc",
  "data.table",
  "Rcpp",
  "ROCR",
  "Bolstad2",
  "MESS",
  "nloptr",
  "cluster",
  "kernlab"
)
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "doCrossValidation.R",
  "calAUPR.R",
  "evalMetrics.R",
  "getCvIndex.R",
  "blm.R",
  "SoftImpute.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp",
                    "fastKgipMat.cpp",
                    "log1pexp.cpp",
                    "sigmoid.cpp",
                    "fastSolve.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)


#######################################


#db <- "CCLE_IC50"
db<-"GDSC_AUC"
#db <- "GDSC_IC50"
switch (db,
        CCLE_IC50 = {
         flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_CCLE_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y) 
          #Y <- t(Y)         
        },
        GDSC_AUC = {
          #cat("ic data\n")
          flush.console()
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_AUC.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
          
        },
        GDSC_IC50 = {
          flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr}\n")
        )


  #### Error: wrong split: 
  ### if idxZeroCol and idxZeroRow = 0;  Y ko thay doi
  ### if idxZeroCol or idxZeroRow = 0 Y thay doi 1 lan
  # Y duoc gan voi matran moi loai cac cot co tong =0 hoac cac hang =0
  # 
  
  # 
  # idxZeroCol <- which(colSums(Y) == 0)
  # if (length(idxZeroCol)!=0){
  #   Y <- Y[, -idxZeroCol]
  #   
  # }
  # 
  #  ## which(colSums(Y) == 0)
  # 
  # idxZeroRow <- which(rowSums(Y) == 0)
  # if (length(idxZeroRow)!=0){
  #   Y <- Y[-idxZeroRow, ]
  #   
  #  }
 
 

## convert to kernel

## do cross-validation
kfold <- 10
numSplit <- 1

## caculate savedFolds

savedFolds <- doCrossValidation(inMat=Y, kfold = kfold, numSplit = numSplit)
  
## for saving results
 AUPRVec <- vector(length = kfold)
# AUCVec <- vector(length = kfold)
# finalResult <- matrix(NA, nrow = numSplit, ncol = 2)
# colnames(finalResult) <- c("AUPR", "AUC")
 CorSplit= vector(length = numSplit)
 RMSESplit= vector(length = numSplit)
 AUCSplit = vector(length = numSplit)
  AUCFold <- vector(length = kfold)
  CorFold <- vector(length = kfold)
  RMSEFold <- vector(length = kfold)
  
  ## parameters
  lambda <- 10#10
  rank <-50
  
  ## main loop
  for (i in 1:numSplit) {
    
    for (j in 1:kfold) {
      cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, 
          "/", kfold, "\n")
      flush.console()
      
      ## training set with the test set links removed
      Yfold <- savedFolds[[i]][[j]][[7]]
      
      ## extract test set
      testSet <- savedFolds[[i]][[j]][[2]]
          
      Ypred <- function.Impute_kfold(x = Yfold, lambdaImp = lambda, rankImp=rank) 
      
      testLabel <- Y[testSet]
      score <- Ypred[testSet]
      #test = unique(testLabel)
      #print(test)
      #result <- calAUPR(testLabel, score)
      #auc= auc_roc(score, testLabel)
      #print(auc)
      #AUPR(score, testLabel, AUPRVec, "KSeeds")
      
      result_cor <- cor(testLabel, score, method="pearson")
      result_rmse <- RMSE(testLabel, score)
      #AUPRVec[j] <- result[1, "aupr"]
      #AUCFold[j] <-auc
      CorFold[j] <- result_cor
      RMSEFold[j] <- result_rmse
      #print(AUCFold[j])
      print(CorFold[j])
      print(RMSEFold[j])
    }
    print("split thu")
    print(i)
    #print(ii)
    #AUCSplit[i]<-mean(AUCFold)
    CorSplit[i] <- mean(CorFold)
    RMSESplit[i] <- mean(RMSEFold)
    
    #print( AUCSplit[i] )
    print(CorSplit[i])
    print( RMSESplit[i] )
   
  }
  #AUC <- mean(AUCSplit)
  COR <- mean(CorSplit)
  RMSE <- mean(RMSESplit)
  
  print(" Means of Colleration ")
  #print(AUC)
  print(COR)
  print(RMSE)
  
  # save to file
  curDate <- format(Sys.time(), format = "%Y-%m-%d")
  curTime <- format(Sys.time(), format =  "%H.%M.%S")
  savedFileName <- paste0(db, "_full", curDate, "_", curTime, "_cor", COR, "+-", "_rmse", RMSE,".RData")
  cat("\n\n")
  print(savedFileName)
 
  save.image(file = paste0("../R.data/",savedFileName))
  
  
  
  




