setwd("d:/Setup/BioInformatic/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM_Copy/Code/")
pkgs <- c(
  "matrixcalc",
  "data.table",
  "Rcpp",
  "ROCR",
  "Bolstad2",
  "MESS",
  "nloptr",
  "cluster",
  "kernlab",
  "plyr",
  "doSNOW",
  "parallel",
  "snow",
  "foreach",
  "interators",
  "gtools",
  "DTHybrid",
  "data.table",
  "mltools",
  "Metrics"
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
  "nbi_ProfS_function.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp",
                    "fastKgipMat.cpp",
                    "log1pexp.cpp",
                    "sigmoid.cpp",
                    "fastSolve.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)

#Read data.frame
db <- "ProbS_CCLE_IC50"
#db<-"ProbS_GDSC_AUC"
#db <- "ProbS_GDSC_IC50"
switch (db,
        ProbS_CCLE_IC50 = {
          flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_CCLE_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y) 
          #Y <- t(Y)         
        },
        ProbS_GDSC_AUC = {
          #cat("ic data\n")
          flush.console()
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_AUC.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
          
        },
        ProbS_GDSC_IC50 = {
          flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
        },
        stop("db should be one of the follows: 
               {,GDSC_IC50,GDSC_AUC,CCLE_IC50}\n")
)
  kfold <- 10
  numSplit <- 1#5
  
  savedFolds <- doCrossValidation(inMat=Y, kfold = 10, numSplit = 1)
  
  ## for saving results
  CorSplit= vector(length = numSplit)
  RMSESplit= vector(length = numSplit)
  CorSDSplit= vector(length = numSplit)
  RMSESDSplit= vector(length = numSplit)
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
      
      Ypred <- nbi_ProbS(Y = Yfold, lambdaY=0.5) 
      
      
      testLabel <- Y[testSet]
      score <- Ypred[testSet]
      
      result_cor <- cor(testLabel, score, method="pearson")
      result_rmse <- rmse(testLabel, score)
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
    CorSDSplit[i] <- sd(CorFold)
    RMSESDSplit[i] <- sd(RMSEFold)
    
    #print( AUCSplit[i] )
    print(CorSplit[i])
    print( RMSESplit[i] )
  }
  
  COR <- mean(CorSplit)
  RMSE <- mean(RMSESplit)
  CORSD <- mean(CorSDSplit)
  RMSESD <- mean(RMSESDSplit)
  
  print(" Means of Colleration ")
  #print(AUC)
  print(COR)
  print(RMSE)
  
  print(CORSD)
  print(RMSESD)
  
  
  # save to file
  curDate <- format(Sys.time(), format = "%Y-%m-%d")
  curTime <- format(Sys.time(), format =  "%H.%M.%S")
  savedFileName <- paste0(db, "_full", curDate, "_", curTime, "_cor", COR, "+-", "_rmse", RMSE,".RData")
  cat("\n\n")
  print(savedFileName)
  
  save.image(file = paste0("../R.data/",savedFileName))
