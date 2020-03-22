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
  "parallel",
  "snow",
  "foreach",
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


# ## sourceCPP required C++ files
# cppSourceNames <- c("fastKF.cpp",
#                     "fastKgipMat.cpp",
#                     "log1pexp.cpp",
#                     "sigmoid.cpp",
#                     "fastSolve.cpp")
# cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)


sd <- read.table("/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/DOIDSimMat.txt",check.names = FALSE)
sd <- as.matrix(sd)
#st <- read.table("/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/EnhSimMat_SharedGene.txt",check.names = FALSE)
st <- read.table("/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/EnhSimMat_Seq.txt",check.names = FALSE)
st <- as.matrix(st)
Y <- read.table("/Users/admin/Manuscripts/75 Mat4DiseaseEnhancerPrediction (DEP)/Data/UnMatched/DOIDEnhMat.txt",check.names = FALSE)
Y <- as.matrix(Y) 

db <- "MatchedHeterNet"
#db <- "MatchedHomoNet" 

switch (db,
        MatchedHeterNet = {
          flush.console()

          sharedRows <- intersect(rownames(sd), rownames(Y))
          sharedCols<-intersect(colnames(st), colnames(Y))
          sdm<-sd[sharedRows, sharedRows]
          stm<-st[sharedCols, sharedCols]
          Ym<-Y[sharedRows,sharedCols]
        },
        MatchedHomoNet = {
          flush.console()
          
          sharedRows <- intersect(rownames(sd), rownames(Y))
          sharedCols<-intersect(colnames(st), colnames(Y))
          
          sdm<-sd
          stm<-st
          Ym<-matrix(0, nrow = dim(sd)[1], ncol = dim(st)[1])
          rownames(Ym) <-rownames(sd)
          colnames(Ym) <-colnames(st)
          
          Ym[sharedRows,sharedCols]<-Y[sharedRows,sharedCols]
        },
        stop("db should be one of the follows: 
               {,GDSC_IC50,GDSC_AUC,CCLE_IC50}\n")
)

NuRows=dim(Ym)[1]


AUROCAll <- vector()
AUPRAll <- vector()
AUROCAll.SD <- vector()
AUPRAll.SD <- vector()

FPRAll <-vector()
TPRAll <-vector()

RecAll <-vector()
PreAll <-vector()

for (ii in 1:NuRows){
    Oneidxs <-which(Ym[ii,]>0)
    kfold <-length(Oneidxs)
    
    if(kfold<2) next()
    
    AUROCDrug <- vector(length = kfold)
    AUPRDrug <- vector(length = kfold)
    
    FPRDrug <-vector()
    TPRDrug <-vector()
    
    RecDrug <-vector()
    PreDrug <-vector()
    
    for (j in 1:kfold) {
      print(paste("Drug",ii,":", rownames(Ym)[ii], "kfold:", j, "/", kfold))
      flush.console()
      
      testLabel <-Ym[ii,]
      
      Yfold <- Ym
      Yfold[ii,Oneidxs[j]] <-0
      
      #Ypred <- computeRecommendation(Yfold, lambda=0.5, alpha=0.5, S=NA, S1=NA, cl=NA)
      Ypred <- computeRecommendation(Yfold, lambda=0.5, alpha=0.5, S=sdm, S1=stm, cl=NA)
          
      score <- Ypred[ii,]
      
      ########## Start Evaluation
      pred <- prediction(score,testLabel)
      #AUCROC
      
      perfROC <- performance(pred,"tpr","fpr")
      
      
      fpr<-perfROC@x.values[[1]]
      tpr<-perfROC@y.values[[1]]
      
      # print(length(fpr))
      # print(length(tpr))
      

      FPRDrug<-c(FPRDrug,fpr)
      TPRDrug<-c(TPRDrug,tpr)
      
      # plot(perfROC,colorize=FALSE)
      
      perfAUC <- performance(pred,"auc")
      auc<-perfAUC@y.values[[1]]
      
      #AUCPR
      perfPR <- performance(pred,"prec","rec")
      # plot(perfPR,colorize=FALSE)
      rec<-perfPR@x.values[[1]]#Recall
      pre<-perfPR@y.values[[1]]#Precision
      # aupr <- try(MESS::auc(rec, pre, type = 'spline'), silent = TRUE)
      aupr <- Bolstad2::sintegral(rec, pre)$int
      
      RecDrug <-c(RecDrug, rec)
      PreDrug <-c(PreDrug, pre)
      
      
      
    
      ########## End Evaluation
      print(auc)
      AUROCDrug[j] <- auc
      
      print(aupr)
      AUPRDrug[j] <- aupr
      
    }
    print(paste("Drug", ii, ":", mean(AUROCDrug)))
    plot(FPRDrug, TPRDrug)
    
    AUROCAll <- c(AUROCAll, mean(AUROCDrug))
    AUROCAll.SD <- c(AUROCAll.SD, sd(AUROCDrug))
    
    AUPRAll <- c(AUPRAll, mean(AUPRDrug))
    AUPRAll.SD <- c(AUPRAll.SD, sd(AUPRDrug))
    
    FPRAll<-c(FPRAll,FPRDrug)
    TPRAll<-c(TPRAll,TPRDrug)
    
    RecAll<-c(RecAll,RecDrug)
    PreAll<-c(PreAll,PreDrug)
}

AUROCavg <- mean(AUROCAll)
AUROCavg.sd <- sd(AUROCAll)

AUPRavg <- mean(AUPRAll)
AUPRavg.sd <- sd(AUPRAll)

print("Final results")

print(paste("AUROCavg:",AUROCavg,"AUROCavg.sd:",AUROCavg.sd))

print(paste("AUPRavg:", AUPRavg, "AUPRavg.sd:",AUPRavg.sd))

plot(FPRAll, TPRAll)
plot(RecAll, PreAll)

# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "_", curDate, "_", curTime, "_AUROC_", AUROCavg, "+-", "_SD_", AUROCavg.sd,".RData")
cat("\n\n")
print(savedFileName)
