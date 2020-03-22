#setwd("d:/Setup/BioInformatic/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM_Copy/Code")
#setwd("C:/Users/Mrs Giang/OneDrive - Hanoi University of Science and Technology/R/Dataset/chemogenomicAlg4DTIpred-master/chemogenomicAlg4DTIpred-master/aucAUPR/BLM/")
rm(list = ls())
library(caret)
library(Metrics)
library(mltools)
library(data.table)
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


db <- "CCLE_IC50"
#db<-"GDSC_AUC"
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
             {,GDSC_IC50,GDSC_AUC,CCLE_IC50}\n")
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


    
    ## parameters
    lambda <- 10#10
    rank <-7
    m=as.matrix(Y)
    
    nrow=dim(m)[1]
    
    CORVec <- vector(length=nrow)
    RMSEVec <- vector(length = nrow)
    AUCVec <- vector(length = nrow)
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
    
    ## convert to kernel
    
    ## do cross-validation
		kfold <- 10
		numSplit <- 1#5
		
		## DT-Hybrid method
		#savedFolds <- doCrossVal(Y, nfold = kfold, nsplit = numSplit)
		
		savedFolds <- doCrossValidation(inMat=Y, kfold = 10, numSplit = 1)
		
		CorSplit= vector(length = numSplit)
		RMSESplit= vector(length = numSplit)
		
		AUCSplit = vector(length = numSplit)
		
		## for saving results
		## main loop
		for (i in 1:numSplit) {
		  AUCFold <- vector(length = kfold)
		  CorFold <- vector(length = kfold)
		  RMSEFold <- vector(length = kfold)
		  
		  
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
			
			  # use SoftImpute to caculate for kfold method
			 
			
			Ypred <- function.Impute_kfold(x=Yfold, rankImp = rank,lambdaImp = lambda) 
		
			# set testSet from Ypredition
			
				n<-length(savedFolds[[i]][[j]][[3]])
				testSet<-vector(length=n)
				for(k in 1:n){
				  r<-savedFolds[[i]][[j]][[3]][k]
				  c<-savedFolds[[i]][[j]][[4]][k]
				  testSet[k]<-Ypred[ii,c]
				}
			
				#auc= auc_roc(testSet, testLabel)
				result_cor <- cor(testLabel, testSet, method="pearson")
				result_rmse <- rmse(testLabel, testSet)
				
  			#result <- calAUPR(testLabel, score)
  			#result <- cor(testLabel, testSet , method = "pearson")
  			
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
		print("hang thu")
		print(ii)
		CORVec[ii]=mean(CorSplit)
		RMSEVec[ii] =mean(RMSESplit)
		#AUCVec[ii] =mean(AUCSplit)
		#print( AUCVec[i] )
		print(CORVec[i])
		print( RMSEVec[i] )
		
    }
    
  COR <- mean(CORVec)
  RMSE <- mean(RMSEVec)
 
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

