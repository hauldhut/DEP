setwd("/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Rcode/dthybrid/")
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
  "doCrossValidation.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


## sourceCPP required C++ files
# cppSourceNames <- c("fastKF.cpp",
#                     "fastKgipMat.cpp",
#                     "log1pexp.cpp",
#                     "sigmoid.cpp",
#                     "fastSolve.cpp")
# cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)

#Read data.frame
#db <- "DTH_CCLE_IC50"
#db<-"DTH_GDSC_AUC"
#db <- "DTH_GDSC_IC50"

#db <- "Mat_SeqBased_Matched_w_HomoNets"
db <- "Mat_SeqBased" #Matched_w_HeterNets
#db <- "Mat_SharedGeneBased" #Matched_w_HeterNets
switch (db,
        Mat_SeqBased_Matched_w_HomoNets = {
          flush.console()
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_w_HomoNets/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_w_HomoNets/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_w_HomoNets/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)  
        },
        Mat_SeqBased = {
          flush.console()
          
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)  
        },
        Mat_SharedGeneBased = {
          flush.console()
          
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SharedGeneBased/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SharedGeneBased/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SharedGeneBased/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)    
        },
        DTH_CCLE_IC50 = {
          flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_CCLE_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y) 
          #Y <- t(Y)         
        },
        DTH_GDSC_AUC = {
          #cat("ic data\n")
          flush.console()
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_AUC.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
          
        },
        DTH_GDSC_IC50 = {
          flush.console()
          
          Y <- read.table('../Dataset/PreviousData/DrugCelllineResponse_GDSC_IC50.txt_AdjMat.txt')
          Y <- as.matrix(Y)
          #Y <- t(Y)
        },
        
        stop("db should be one of the follows: 
               {,GDSC_IC50,GDSC_AUC,CCLE_IC50}\n")
)
  kfold <- 5
  numSplit <- 1#5
  
  savedFolds <- doCrossValidation(inMat=Y, kfold = kfold, numSplit = numSplit)
  
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
      
      #Ypred <- computeRecommendation(Yfold, lambda=0.5, alpha=0.5, S=NA, S1=NA, cl=NA)
      Ypred <- computeRecommendation(Yfold, lambda=0.5, alpha=0.5, S=sd, S1=st, cl=NA)
      
      testLabel <- Y[testSet]
      score <- Ypred[testSet]
      
      ########## Start Evaluation
      pred <- prediction(score,testLabel)
      #AUCROC
      
      perfROC <- performance(pred,"tpr","fpr")
      plot(perfROC,colorize=FALSE)
      
      perfAUC <- performance(pred,"auc")
      auc<-perfAUC@y.values[[1]]
      
      #AUCPR
      perfPR <- performance(pred,"prec","rec")
      plot(perfPR,colorize=FALSE)
      Recall<-perfPR@x.values[[1]]#Recall
      Precision<-perfPR@y.values[[1]]#Precision
      #aupr_spline <- try(MESS::auc(Recall, Precision, type = 'spline'), silent = TRUE)
      aupr_simpson <- Bolstad2::sintegral(Recall, Precision)$int
      #print(aupr_simpson)
      
      result_cor <- auc
      result_rmse <- aupr_simpson
      ########## End Evaluation
      
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
  CORSD <- sd(CorSDSplit)
  RMSESD <- sd(RMSESDSplit)
  
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
  
