##Neu loi lien quan den gfortran khi bien dich file *.cpp.
##Download va cai dat gfortran tai http://gcc.gnu.org/wiki/GFortranBinaries#MacOS

setwd("/Users/admin/Manuscripts/42 Disease-associated Enhancers/Rcode/chemogenomicAlg4DTIpred/aucAUPR/dthybrid")

rm(list = ls())

## current data set name
db <- "enh"#en/ic/GPCR/ic/nr/kd/miRWalk/TargetScan

switch (db,
        miRNA = {
          cat("miRNA data\n")
          flush.console()
          sd <- read.table("DisSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("miRNASimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("miRNADisMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        miRWalk = {
          cat("miRWalk data\n")
          flush.console()
          sd <- read.table("miRWalk_DisSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("miRWalk_miRNASimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("miRWalk_miRNADisMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        TargetScan = {
          cat("TargetScan data\n")
          flush.console()
          sd <- read.table("TargetScan_DisSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("TargetScan_miRNASimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("TargetScan_miRNADisMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        en = {
          cat("en data\n")
          flush.console()
          sd <- read.table("e_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("e_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("e_admat_dgc.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        ic = {
          cat("ic data\n")
          flush.console()
          sd <- read.table("ic_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("ic_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("ic_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        gpcr = {
          cat("gpcr data\n")
          flush.console()
          sd <- read.table("gpcr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("gpcr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("gpcr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        nr = {
          cat("nr data\n")
          flush.console()
          sd <- read.table("nr_simmat_dc.txt")
          sd <- as.matrix(sd)
          st <- read.table("nr_simmat_dg.txt")
          st <- as.matrix(st)
          Y <- read.table("nr_admat_dgc.txt")
          Y <- as.matrix(Y)
          Y <- t(Y)
        },
        kd = {
          cat("kd data\n")
          Y <- read.table("drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt")
          Y[Y <= 30] <- 1
          Y[Y > 30] <- 0
          Y <- as.matrix(Y)
          sd <- read.table("drug-drug_similarities_2D.txt")
          sd <- as.matrix(sd)
          st <- read.table("target-target_similarities_WS_normalized.txt")
          st <- as.matrix(st)
        },
        dr = {
          cat("dr data\n")
          flush.console()
          sd <- read.table("DrugSim-Fdataset.txt")
          sd <- as.matrix(sd)
          st <- read.table("DiseaseSim-Fdataset.txt")
          st <- as.matrix(st)
          Y <- read.table("DiDrA-Fdataset.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        enh = {
          cat("enh data\n")
          flush.console()
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr}\n")
        )
dim(sd)
dim(st)
dim(Y)

if (db == "kd") {
  dim(Y) ## 68 * 442
  dim(sd) ## 68 * 68
  dim(st) ##  442 * 442
  
  
  idxZeroCol <- which(colSums(Y) == 0)
  Y <- Y[, -idxZeroCol]
  st <- st[-idxZeroCol, -idxZeroCol]
  
  ## which(colSums(Y) == 0)
  
  idxZeroRow <- which(rowSums(Y) == 0)
  Y <- Y[-idxZeroRow, ]
  sd <- sd[-idxZeroRow, -idxZeroRow]
  
  which(rowSums(Y) == 0)
  which(colSums(Y) == 0)
  
  dim(Y)  ## 65 373
  dim(sd) ## 65 65
  dim(st) ## 373 373
  
  sd[1:3, 1:3]
  st[1:3, 1:3]
  Y[1:3, 1:3]
}

## load required packages
pkgs <- c("matrixcalc", "data.table", "Rcpp", "ROCR", "Bolstad2", "MESS")
rPkgs <- lapply(pkgs, require, character.only = TRUE)

## source required R files
rSourceNames <- c(
  "doCVPositiveOnly.R",
  "doCVPositiveOnly3.R",
  "doCrossVal.R",
  "calAUPR.R",
  "evalMetrics.R",
  "Recommendation.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

## do cross-validation
kfold <- 5
numSplit <- 1

## DT-Hybrid method
#savedFolds <- doCrossVal(Y, nfold = kfold, nsplit = numSplit)
nlinks=2
rowOK<-which(rowSums(Y)>=nlinks)
colOK<-which(colSums(Y)>=nlinks)

matOK<-matrix(0,nrow = nrow(Y),ncol=ncol(Y))
matOK[rowOK, colOK]<-Y[rowOK, colOK]

OneIdx<-which(Y!=0)
onePos <- which(Y!= 0, arr.ind = TRUE)#Tra ve index dang row, col
ZeroIdx<-which(Y==0)

onePos <- onePos[onePos[, "row"] %in% rowOK, ]
onePos <- onePos[onePos[, "col"] %in% rowOK, ]

# NuRevOneLinks<-round(length(OneIdx)/kfold)
# NuRevZeroLinks<-round(length(ZeroIdx)/kfold)

dim(onePos)
onePos
NuOKElement<-nrow(onePos)
NuRevElement<-ceiling(NuOKElement/kfold)

# main loop
par(mfrow=c(3,3))# 4 figures arranged in 1 rows and 2 columns
library(ROCR)

  testLabelAll<-c()
  scoreAll<-c()
  for (j in 1:kfold) {
    cat("kfold:", j, "/", kfold, "\n")
    
    RevIdx<-sample(NuOKElement, NuRevElement)
    rowRevIdx<-onePos[RevIdx,1]
    colRevIdx<-onePos[RevIdx,2]
    #Convert array index to linear index
    RevOneIdx<-(colRevIdx-1)*nrow(Y)+rowRevIdx
    
    #RevOneIdx<-sample(OneIdx, NuRevOneLinks)
    #RevZeroIdx<-sample(ZeroIdx, NuRevZeroLinks)
    RevZeroIdx<-ZeroIdx
    
    RevIdx<-c(RevOneIdx,RevZeroIdx)
    
    Yfold<-Y
    Yfold[rowRevIdx,colRevIdx]<-0
    
    ## row is drug, col is target
    if (db == "en") {
      theAl <- 0.4
    } else if (db == "ic") {
      theAl <- 0.3
    } else if (db == "gpcr") {
      theAl <- 0.2
    } else {
      theAl <- 0.4
    }
    
    Ypred <- computeRecommendation(A = Yfold, lambda = 0.5, 
                                   alpha = theAl, S = sd, S1 = st)
    
    testLabel <- Y[RevIdx]
    score <- Ypred[RevIdx]
    print(length(testLabel))
    
    
    pred <- prediction(score,testLabel)
    perf <- performance(pred,"tpr","fpr")
    plot(perf,colorize=FALSE)
    perf <- performance(pred,"auc")
    auc<-perf@y.values[[1]]
    print(auc)
    
    
    testLabelAll<-c(testLabelAll,testLabel)
    scoreAll<-c(scoreAll,score)
    
  }
  print(length(testLabelAll))
  
  pred <- prediction(scoreAll,testLabelAll)
  
  perf <- performance(pred,"tpr","fpr")
  plot(perf,colorize=FALSE)
  FPR<-perf@x.values
  TPR<-perf@y.values
  perf <- performance(pred,"auc")
  auc<-perf@y.values[[1]]
  print(auc)
  
  perf <- performance(pred,"prec","rec")
  plot(perf,colorize=FALSE)
  Recall<-perf@x.values#Recall
  Precision<-perf@y.values#Precision
  
  aupr_spline <- try(MESS::auc(Recall, Precision, type = 'spline'), silent = TRUE)
  aupr_simpson <- Bolstad2::sintegral(Recall, Precision)$int
  print(aupr_simpson)
