#setwd("d:/Setup/BioInformatic/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM_Copy/")
#library(netpredictor)
library(igraph)
library(netpredictor)
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
    "iterators",
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
  
  #Read data.frame
  #db <- "nbi_CCLE_IC50"
  db<-"nbi_GDSC_AUC"
  #db <- "nbi_GDSC_IC50"
  switch (db,
          nbi_CCLE_IC50 = {
            flush.console()
            
            Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_CCLE_IC50.txt_AdjMat.txt')
            Y <- as.matrix(Y) 
            #Y <- t(Y)         
          },
          nbi_GDSC_AUC = {
            #cat("ic data\n")
            flush.console()
            Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_AUC.txt_AdjMat.txt')
            Y <- as.matrix(Y)
            #Y <- t(Y)
            
          },
          nbi_GDSC_IC50 = {
            flush.console()
            
            Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_IC50.txt_AdjMat.txt')
            Y <- as.matrix(Y)
            #Y <- t(Y)
          },
          stop("db should be one of the follows: 
               {,GDSC_IC50,GDSC_AUC,CCLE_IC50}\n")
          )
  
  m=as.matrix(Y)
  
  nrow=dim(m)[1]
  
  CORVec <- vector(length=nrow)
  RMSEVec <- vector(length = nrow)
  CORSDVec <- vector(length = nrow)
  RMSESDVec <- vector(length = nrow)
  mt<-c()
  for (ii in 1:nrow){
    Y <- t(m[ii,])
    x=ii-1
    if(x>1){
      Xii <- as.matrix(m[1:x,])
      
    }else{
      Xii <- t(m[1,])
      
    }
    
    
    z=ii+1
    if (z<nrow){
      
      Zii <-as.matrix(m[z:nrow,])
      
    } else {
      
      Zii=t(m[nrow,])
    }
    
  ## do cross-validation
  kfold <- 10
  numSplit <- 1
  
  ## k-fold method
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
      #testLabel <- savedFolds$split_1$fold_1$testLabel
      testLabel <- savedFolds[[i]][[j]][[1]]
      testSet <- savedFolds[[i]][[j]][[2]]
      
      #print(length(testSet))
      #merge matrix Yfold and submatrix
      if(ii==1){
        
        Yfold <-rbind(Yfold, Zii)
      } else if (ii==nrow) {
        
        Yfold <-rbind(Xii, Yfold)
      } else {
        
        Yfold <- rbind(Xii,Yfold, Zii)				
      }
     
	    Ypred <- nbiNet(Yfold, alpha=0.5, lamda=0.5, s1=NA, s2=NA,format = "matrix")
      
	    n<-length(savedFolds[[i]][[j]][[3]])
	    testSet<-vector(length=n)
	    for(k in 1:n){
	      r<-savedFolds[[i]][[j]][[3]][k]
	      c<-savedFolds[[i]][[j]][[4]][k]
	      testSet[k]<-Ypred[ii,c]
	    }
      
	    result_cor <- cor(testLabel, testSet, method="pearson")
	    result_rmse <- rmse(testLabel, testSet)
	    
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
  print("hang thu")
  print(ii)
  CORVec[ii]=mean(CorSplit)
  RMSEVec[ii] =mean(RMSESplit)
  CORSDVec[ii]=mean(CorSDSplit)
  RMSESDVec[ii] =mean(RMSESDSplit)
  print(CORVec[i])
  print( RMSEVec[i] )
  }
  COR <- mean(CORVec)
  RMSE <- mean(RMSEVec)
  CORSD <- mean(CORSDVec)
  RMSESD <- mean(RMSESDVec)
  
  print(" Means of Colleration ")
  #print(AUC)
  print(COR)
  print(RMSE)
  
  print(CORSD)
  print(RMSESD)
  
  
  # save to file
  curDate <- format(Sys.time(), format = "%Y-%m-%d")
  curTime <- format(Sys.time(), format =  "%H.%M.%S")
  savedFileName <- paste0(db, "_Row_", curDate, "_", curTime, "_cor", COR, "+-", "_rmse", RMSE,".RData")
  cat("\n\n")
  print(savedFileName)
  
  save.image(file = paste0("../R.data/",savedFileName))


