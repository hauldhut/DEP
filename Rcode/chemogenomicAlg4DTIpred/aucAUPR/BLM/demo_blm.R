##Neu loi lien quan den gfortran khi bien dich file *.cpp.
##Download va cai dat gfortran tai http://gcc.gnu.org/wiki/GFortranBinaries#MacOS

setwd("/Users/admin/Manuscripts/42 Disease-associated Enhancers/Rcode/chemogenomicAlg4DTIpred/aucAUPR/BLM")

rm(list = ls())

## current data set name
db <- "enh"#en/ic/gpcr/ic/nr/kd/dr
print(db)
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
          sd <- read.table("DrugSim.txt")
          sd <- as.matrix(sd)
          st <- read.table("DiseaseSim.txt")
          st <- as.matrix(st)
          Y <- read.table("DiDrA.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        gwas = {
          cat("gwas data\n")
          flush.console()
          sd <- read.table("/Users/admin/Data/GWAS/CAUSALdb/MeSHIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/GWAS/CAUSALdb/SNPSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/GWAS/CAUSALdb/SNPMeSHIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        enh = {
          cat("enh data\n")
          flush.console()
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_2_HomoNets/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_2_HomoNets/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased_Matched_2_HomoNets/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr, kd, dr}\n")
        )

dim(Y)

if (db == "miRNA") {
  dim(Y) ## 68 * 442
  dim(sd) ## 68 * 68
  dim(st) ##  442 * 442
  
  
  idxZeroCol <- which(colSums(Y) == 0)
  Y <- Y[, -idxZeroCol]#Loai bo cac Target chua biet bat ky associated Drug nao trong matrix Y
  st <- st[-idxZeroCol, -idxZeroCol]#Loai bo tat ca cac Target do trong matrix st
  
  ## which(colSums(Y) == 0)
  
  idxZeroRow <- which(rowSums(Y) == 0)
  Y <- Y[-idxZeroRow, ]#Loai bo cacs Drug chua biet bat ky associated Target nao torng Y
  sd <- sd[-idxZeroRow, -idxZeroRow]##Loai bo tat ca cac Drug do trong matrix sd
  
  which(rowSums(Y) == 0)#Check lai xem 2 thao tac tren da work chua
  which(colSums(Y) == 0)#Check lai xem 2 thao tac tren da work chua
  
  #Kiem tra kich thuoc matrix sau khi loai bo cac row/col toan 0
  dim(Y)  ## 65 373 
  dim(sd) ## 65 65
  dim(st) ## 373 373
  
  sd[1:3, 1:3]
  st[1:3, 1:3]
  Y[1:3, 1:3]
}

dim(Y)

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
  "doCrossVal.R",
  "calAUPR.R",
  "evalMetrics.R",
  "getCvIndex.R",
  "blm.R"
)
rSN <- lapply(rSourceNames, source, verbose = FALSE)


## sourceCPP required C++ files
cppSourceNames <- c("fastKF.cpp",
                    "fastKgipMat.cpp",
                    "log1pexp.cpp",
                    "sigmoid.cpp",
                    "fastSolve.cpp")
cppSN <- lapply(cppSourceNames, sourceCpp, verbose = FALSE)


## convert to kernel
# Preprocessing
isKernel <- TRUE
if (isKernel) {
  #Xu ly voi matrix sd
  #Bien thanh matrix doi xung neu chua doi xung
  if (!isSymmetric(sd)) {
    sd <- (sd + t(sd)) / 2
  }
  
  #Bien thanh ma tran duong, ban xac dinh
  epsilon <- 0.1
  while (!is.positive.semi.definite(sd)) {
    sd <- sd + epsilon * diag(nrow(sd))
  }
  
  #Xu ly voi matrix st
  if (!isSymmetric(st)) {
    st <- (st + t(st)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(st)) {
    st <- st + epsilon * diag(nrow(st))
  }
}

#Chuyen vi vai tro cua Drug va Target. Do vai tro cua Drug va Target la tuong duong ve mat thuat toan
Y <- t(Y)
#swap sd <-> st
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

## do cross-validation
kfold <- 5
numSplit <- 5#number of trials. Repeat the kFold process numSplit times

## DT-Hybrid method
print("Hello World")
print(dim(Y))
#print(Y)

#Phan tach du lieu Y thanh 2 phan (1: Yfold (training set); 2: testSet
#Lap lai numSplit trials, moi trial thi lam kfold.

savedFolds <- doCrossVal(Y, nfold = kfold, nsplit = numSplit)

## for saving results
#Moi fold se output ra AUPR va AUC --> kfold se tao ra 1 vector gom kfold gia tri nhu vay
AUPRVec <- vector(length = kfold)
AUCVec <- vector(length = kfold)

#Tao matrix chua file result. So dong bang so trials (numSplit), so cot chinh la AUPR va AUC
finalResult <- matrix(NA, nrow = numSplit, ncol = 2)
colnames(finalResult) <- c("AUPR", "AUC")


## parameters
cc <- 1

## main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, 
        "/", kfold, "\n")
    flush.console()
    
    ## training set with the test set links removed
    Yfold <- savedFolds[[i]][[j]][[1]]
    cat(dim(Yfold),"\n")
    
    Kd <- sd
    Kt <- st
    
    ## extract test set
    testSet <- savedFolds[[i]][[j]][[2]]
    dim(testSet)    
    
    Ypred <- blm(Yfold = Yfold, Kd = Kd, Kt = Kt, cc = cc) 
    
    testLabel <- Y[testSet]#observed label: 1 (co interaction) va 0 (unknown)
    score <- Ypred[testSet]#Ypred duoc xay dung dua tren Yfold (Real value: -1 -> +1)
    
    #Tinh AUPR va AUC
    result <- calAUPR(testLabel, score)
    
    AUPRVec[j] <- result[1, "aupr"]
    AUCVec[j] <- result[1, "auc"]
  }
  AUPR <- mean(AUPRVec)
  AUC <- mean(AUCVec)
  finalResult[i, "AUPR"] <- AUPR
  finalResult[i, "AUC"] <- AUC
}

#Lam tron den 3 so sau dau . thap phan
auc <- round(mean(finalResult[, "AUC"]), 3)
aucSD <- round(sd(finalResult[, "AUC"]), 3)

aupr <- round(mean(finalResult[, "AUPR"]), 3)
auprSD <- round(sd(finalResult[, "AUPR"]), 3)



# save to file
curDate <- format(Sys.time(), format = "%Y-%m-%d")
curTime <- format(Sys.time(), format =  "%H.%M.%S")
savedFileName <- paste0(db, "_", curDate, "_", curTime, "_auc", auc, "+-", aucSD, "_aupr", aupr, "+-", auprSD, ".RData")
cat("\n\n")
print(savedFileName)
# save.image(file = savedFileName)









