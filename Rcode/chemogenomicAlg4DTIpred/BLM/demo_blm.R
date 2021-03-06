##Neu loi lien quan den gfortran khi bien dich file *.cpp.
##Download va cai dat gfortran tai http://gcc.gnu.org/wiki/GFortranBinaries#MacOS

#setwd("/Users/admin/Tools/DTI/chemogenomicAlg4DTIpred/BLM/")

rm(list = ls())

## current data set name
db <- "enh"#en/ic/gpcr/ic/nr/dr

switch (db,
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
        gwas = {
          cat("gwas data\n")
          flush.console()
          sd <- read.table("/Users/admin/Data/GWAS/CAUSALdb/LD=0.2/MeSHIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/GWAS/CAUSALdb/LD=0.2/SNPSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/GWAS/CAUSALdb/LD=0.2/SNPMeSHIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        enh = {
          cat("enh data\n")
          flush.console()
          sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat/DOIDSimMat.txt")
          sd <- as.matrix(sd)
          st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat/EnhSimMat.txt")
          st <- as.matrix(st)
          Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat/EnhDOIDMat.txt")
          Y <- as.matrix(Y) 
          Y <- t(Y)         
        },
        stop("db should be one of the follows: 
             {en, ic, gpcr, nr}\n")
        )

getwd()
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
isKernel <- TRUE
if (isKernel) {
  if (!isSymmetric(sd)) {
    sd <- (sd + t(sd)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(sd)) {
    sd <- sd + epsilon * diag(nrow(sd))
  }
  if (!isSymmetric(st)) {
    st <- (st + t(st)) / 2
  }
  epsilon <- 0.1
  while (!is.positive.semi.definite(st)) {
    st <- st + epsilon * diag(nrow(st))
  }
}

Y <- t(Y)
tmp <- sd
sd <- st
st <- tmp

Y[1:3, 1:3]

## do cross-validation
kfold <- 5
numSplit <- 5

## DT-Hybrid method
savedFolds <- doCVPositiveOnly3(Y, kfold = kfold, numSplit = numSplit)

## saving results
resMetrics <- matrix(NA, nrow = kfold, ncol = 1)
colnames(resMetrics) <- c("MPR")
resMetrics <- as.data.frame(resMetrics)
finalResult <- vector("list", length = numSplit)

## parameters
cc <- 1

## main loop
for (i in 1:numSplit) {
  for (j in 1:kfold) {
    cat("numSplit:", i, "/", numSplit, ";", "kfold:", j, 
        "/", kfold, "\n")
    flush.console()
    
    ## training set with the test set links removed
    Yfold <- savedFolds[[i]][[j]][[6]]

    Kd <- sd
    Kt <- st
    
    ## extract test set
    testSet <- savedFolds[[i]][[j]][[1]]
    knownDrugIndex <- savedFolds[[i]][[j]][[4]]
    knownTargetIndex <- savedFolds[[i]][[j]][[5]] 
    testIndexRow <- savedFolds[[i]][[j]][[2]]      
    testIndexCol <- savedFolds[[i]][[j]][[3]]      
    
    Ypred <- blm(Yfold = Yfold, Kd = Kd, Kt = Kt, cc = cc) 
    
    ## result
    result2 <- evalMetrics(Ypred = Ypred,testSet = testSet)
    resMetrics[j, ] <- result2
  }
  finalResult[[i]] <- resMetrics
}

# combine result
resCom <- as.data.frame(data.table::rbindlist(finalResult))
resMean <- colMeans(resCom)
se <- sqrt(var(resCom[, 1]) / length(resCom[, 1]))
cat("BLM:", "MPR =", round(resMean, 3), "+\\-", round(se, 3),  "\n")
flush.console()

savedName <- paste0(db, "_BLM.RData")
save.image(savedName)









