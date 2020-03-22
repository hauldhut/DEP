#library(netpredictor)
library(igraph)
#net.perf<- function(A,S1,S2,restart=0.8,alpha=0.5,lamda=0.5,relinks=100,numT=2,norm="laplace",Calgo = c("rwr","nbi","netcombo","all")){
sd <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/DOIDSimMat.txt")
sd <- as.matrix(sd)
st <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhSimMat.txt")
st <- as.matrix(st)
Y <- read.table("/Users/admin/Data/Enhancer/DiseaseEnhancer/Mat_SeqBased/EnhDOIDMat.txt")
Y <- as.matrix(Y) 
#Y <- t(Y)    

A=Y
S1=st
S2=sd
alpha=0.5
restart=0.8
lamda=0.5
kfold=5
relinks=length(which(Y>0))/kfold
numT=0
norm="chen"#normalise=c("laplace","none","chen")
Calgo = c("rwr")


  auctop = numeric()
  aucc = numeric()
  bdr  = numeric()
  efc   = numeric()
  ranks = numeric()
  totallinks = sum(A)
  
  m = dim(A)[1] ## rows for targets 
  n = dim(A)[2] ## columns for drugs
  
  if (!exists('S1') || !exists('S2')){
    stop("You must submit s1 and s2 matrices.\n")
  }
  
  if (nrow(S1)!=m | ncol(S1) != m){
    stop("Your number of targets does not match with target similarity matrix.\n")
  }
  
  if (nrow(S2)!=n | ncol(S2) != n){
    stop("Your number of targets does not match with target similarity matrix.\n")
  }
  
  
  ## Get the name of the algorithm.
  algo <- "netcombo"#match.arg(Calgo)
  g1 <- graph.incidence(A)
  eg <- get.edgelist(g1)
  c <- data.frame(table(eg[,2]))#Check number of associated elements
  c <- c[c$Freq>numT,]#Chi lay cac drug co number of associated elements >numT
  
  drugnames <- unique(as.character(c$Var1))
  
  ids <- which(eg[,2] %in% drugnames)#Chi lay cac edge (=1) co drug thoa man (...>numT)
  re <- eg[sample(ids,size = relinks,replace=FALSE),]#sampling so edge = relinks to be removed (re)
  
  
  if (totallinks <= relinks){
    stop("Total links removed is less than equal given links to be removes. Give a sensible value.")
  }
  
  SampledGraph <- g1
  for (i in 1:dim(re)[1])
  {
    if (are.connected(SampledGraph, re[i,1], re[i,2])) 
      SampledGraph <- delete.edges(SampledGraph, E(SampledGraph, P=c(re[i,1], re[i,2])))
  }
  g1 = SampledGraph#Remove la chuyen 1 thanh 0 --> size (num of drugs vaf num of targets khong doi)
  #Check
  #length(get.vertex.attribute(g1)$name)
  #dim(S1)[1]+ dim(S2)[1]
  
  Sg_t <- get.incidence(SampledGraph)
  
  #Sg_t <- randomizeMatrix(Sg_t,null.model = "frequency",iterations = 1000)
  
  #mat<-tMat(Sg_t,as.matrix(S1),as.matrix(S2),normalise="laplace")
  
  drugs <- re[,2]
  
  message(sprintf("Detected (%s) drugs & (%s) proteins with (%s) interactions...",n,m,totallinks))
  message(sprintf("Running prediction for (%s) links removed using (%s) .. ",as.character(relinks),as.character(algo)))
  
  # performances <- function(predictR,m,re){
  #   
  #   s1<-predictR[1:m,]
  #   s1<- scale(s1, center=FALSE, scale=colSums(s1,na.rm=TRUE))
  #   s1[is.na(s1)] <- 0
  #   test <- data.frame(re)
  #   for (dis in 1:dim(s1)[2]) {
  #     
  #     drugname = colnames(s1)[dis]
  #     subfr <- test[test$X2==drugname,]
  #     p1name<-as.character(subfr$X1)
  #     id = which(rownames(s1) %in% p1name)
  #     clabel <- rep(0,m)
  #     clabel[id] <- 1
  #     res = cbind(s1[,dis],clabel)
  #     colnames(res)[1] <- "score"
  #     
  #     d <- res[order(-res[,1]),]
  #     ac <- auac(d[,1], d[,2])
  #     au <- auc(d[,1], d[,2])
  #     at <-  auc(d[,1], d[,2],top=0.1)
  #     bd <- bedroc(d[,1], d[,2])
  #     ef <- enrichment_factor(d[,1], d[,2],top=0.1)
  #     aucc <- c(aucc, ac)
  #     bdr <- c(bdr,bd)
  #     efc <- c(efc,ef) 
  #     auctop <- c(auctop,at)
  #     
  #   }
  #   
  #   scores = c(list(auac = mean(aucc),auc= mean(au),auctop = mean(auctop),bdr = mean(bdr),efc = mean(efc)))
  #   return (scores)
  # }
  
  if (algo == "rwr"){
    #par="True"
    message(sprintf("Running RWR Algorithm"))
    mat = biNetwalk(g1,s1=S1,s2=S2,restart=restart,normalise=norm,verbose=T)
    predictR <- mat[,colnames(mat) %in% drugs]#Tra ve gia tri predict cho cac drug ung voi edge da xoa
    
    dim(predictR)#|St|x|drugs|
    #scores <- performances(predictR,m,re)
    
    s1<-predictR[1:m,]#m = dim(A)[1] ## rows for targets 
    s1<- scale(s1, center=FALSE, scale=colSums(s1,na.rm=TRUE))
    s1[is.na(s1)] <- 0
    test <- data.frame(re)#re chua cac edge ban dau truoc khi bi remove
    library(ROCR)
    testLabelAll<-c()
    scoreAll<-c()
    
    for (dis in 1:dim(s1)[2]) {#Tinh cho tung drug trong removed list
      
      drugname = colnames(s1)[dis]
      subfr <- test[test$X2==drugname,]#Chon cac edge tuong ung voi drug do (test chinh la re)
      p1name<-as.character(subfr$X1)#Lay ten cac target co edge
      id = which(rownames(s1) %in% p1name)#Lay id cua cac target do
      clabel <- rep(0,m)#Tao lai vector label voi tat ca cac target deu =0
      clabel[id] <- 1#Tru vi tri cua cac target noi voi drug hien ta (dis)
      res = cbind(s1[,dis],clabel)
      colnames(res)[1] <- "score"
      auc_unordered <- auc(s1[,dis], clabel)
      
      d <- res[order(-res[,1]),]
      
      testLabelAll<-c(testLabelAll,d[,2])
      scoreAll<-c(scoreAll,d[,1])
      
      ac <- auac(d[,1], d[,2])
      au <- auc(d[,1], d[,2])
      at <-  auc(d[,1], d[,2],top=0.1)
      bd <- bedroc(d[,1], d[,2])
      ef <- enrichment_factor(d[,1], d[,2],top=0.1)
      aucc <- c(aucc, ac)
      bdr <- c(bdr,bd)
      efc <- c(efc,ef) 
      auctop <- c(auctop,at)
      
      cat(dis,"\t", drugname,"\auc_unordered\t",auc_unordered,"\tau\t",au,"\n")
    }
    par(mfrow=c(1,2))# 2 figures arranged in 1 rows and 2 columns
    pred <- prediction(scoreAll,testLabelAll)
    
    perf <- performance(pred,"tpr","fpr")
    plot(perf,colorize=FALSE)
    perf@x.values
    perf@y.values
    perf <- performance(pred,"auc")
    auc<-perf@y.values[[1]]
    print(auc)
    
    perf <- performance(pred,"prec","rec")
    plot(perf,colorize=FALSE)
    Precision <- perf@x.values
    Recall<-perf@y.values
    
    aupr_spline <- try(MESS::auc(Recall[[1]], Precision[[1]], type = 'spline'), silent = TRUE)
    aupr_simpson <- Bolstad2::sintegral(Recall[[1]], Precision[[1]])$int
    
    perf <- performance(pred,"auc")
    auc<-perf@y.values[[1]]
    
    scores = c(list(auac = mean(aucc),auc= mean(au),auctop = mean(auctop),bdr = mean(bdr),efc = mean(efc)))
    
    #return (scores)
  }else if (algo == "nbi"){
    message(sprintf("Running NBI Algorithm"))
    #S1 = S1[rownames(S1) %in% rownames(N_M),colnames(S1) %in% rownames(N_M)]
    #S2 = S2[rownames(S2) %in% colnames(N_M),colnames(S2) %in% colnames(N_M)]   
    mat <- nbiNet(Sg_t, lamda=lamda, alpha=alpha, s1=as.matrix(S1), s2=as.matrix(S2),format = "matrix")
    predictR <- mat[,colnames(mat) %in% drugs]
    scores <- performances(predictR,m,re)
    #return (scores)
  }else if(algo == "netcombo"){
    message(sprintf("Running NetCombo Algorithm"))
    #par="True"
    mat1 = biNetwalk(g1,s1=S1,s2=S2,normalise=norm,verbose=T,restart=restart)
    mat2 <- nbiNet(Sg_t,lamda=lamda, alpha=alpha, s1=as.matrix(S1), s2=as.matrix(S2),format = "matrix")
    mat = (mat1+mat2)/2
    predictR <- mat[,colnames(mat) %in% drugs]
    scores <- performances(predictR,m,re)
    #return (scores)
  } else if (algo == "all"){
    
    message(sprintf("Running all the algorithms ..."))
    #par="True"
    mat1 <- biNetwalk(g1,s1=S1,s2=S2,normalise=norm,verbose=T)
    mat2 <- nbiNet(Sg_t, lamda=0.5, alpha=0.5, s1=as.matrix(S1), s2=as.matrix(S2),format = "matrix")
    mat3 <- (mat1+mat2)/2
    predictR1 <- mat1[,colnames(mat1) %in% drugs]
    predictR2 <- mat2[,colnames(mat2) %in% drugs]
    predictR3 <- mat3[,colnames(mat3) %in% drugs]
    
    scores1 <- performances(predictR1,m,re)
    scores2 <- performances(predictR2,m,re)
    scores3 <- performances(predictR3,m,re)        
    
    list1 = list(type = 'rwr',score=scores1)
    list2 = list(type = 'nbi',score=scores2)
    list3 = list(type = 'netcombo',score=scores3)
    scoreList = list(list1,list2,list3)
    #return (scoreList)
    
  }
  
  ## Get the transition matrix for a Bipartite Graph. Input a sequence similarity matrix (s1),
  ## chemical similarity matrix (s2) and drug target adjacency matrix (g1) where rows are 
  ## protein targets and columns as Drugs.
  ## @export
  
  tMat <- function(g1,s1,s2,normalise="chen"){
    
    cl <- makeCluster(detectCores())
    
    g1 <- t(g1)
    seq<- as.matrix(s1)          ## sequence similairty matrix normalised between 0 and 1 
    drugProt <- as.matrix(g1)           ## drug target matrix   
    csim <- as.matrix(s2)         ## drug similarity matrix normalised between 0 and 1
    
    new.drug_drug <- suppressWarnings(snow::parMM(cl,drugProt,t(drugProt)))
    
    new.prot_prot <- suppressWarnings(snow::parMM(cl,t(drugProt),drugProt))
    
    #calculate drug-drug similarity based on shared proteins based on jaccard similarity            
    norm_drug <- jaccard.sim(drugProt)
    
    # Jaccard similarity of two proteins based on shared compounds
    norm_prot <- jaccard.sim(t(drugProt))
    
    # Normalizing the matrices with equal weights
    drug.similarity.final <- 0.5*(csim)+0.5*(norm_drug)
    prot.similarity.final <- 0.5*(seq)+0.5*(norm_prot)
    
    
    
    if(normalise == "laplace"){
      
      D1  <- diag(x=(rowSums(drugProt))^(-0.5))
      D2  <- diag(x=(colSums(drugProt))^(-0.5))   
      
      MTD1 <- suppressWarnings(snow::parMM(cl,D1,g1))
      MTD <- suppressWarnings(snow::parMM(cl,MTD1,D2))
      
      D3  <- diag(x=(rowSums(prot.similarity.final))^(-0.5))
      
      MTT1 <- suppressWarnings(snow::parMM(cl,D3,seq))
      MTT  <- suppressWarnings(snow::parMM(cl,MTT1,D3))
      
      D4  <- diag(x=(rowSums(drug.similarity.final))^(-0.5))
      MDD1  <- suppressWarnings(snow::parMM(cl,D4,csim))
      MDD  <- suppressWarnings(snow::parMM(cl,MDD1,D4))
      
      M1<-cbind(MTT,t(MTD))
      M2<-cbind(MTD,MDD)
      M <- rbind(M1,M2)
      M <- as.matrix(M)
      M[is.na(M)]<-0
      n =c(colnames(g1),rownames(g1))
      rownames(M) <- n
      colnames(M) <- n
      # Returning the final matrix 
      return(as.matrix(M))
    }    
    if(normalise == "chen"){
      
      ADT<-colSums(drugProt)
      ATD<-rowSums(drugProt)
      Sd<-rowSums(drug.similarity.final)
      St<-rowSums(prot.similarity.final)
      MDT<-mat.or.vec(nrow(drugProt),nrow(seq))
      MTD<-mat.or.vec(nrow(seq),nrow(drugProt))
      print (dim(MDT))
      print (dim(MTD))
      
      A<-t(drugProt)
      
      D1  <- diag(x=(rowSums(drugProt))^(-0.5))
      D2  <- diag(x=(colSums(drugProt))^(-0.5))   
      
      #MTD1 <- suppressWarnings(snow::parMM(cl,D1,g1))
      ##MTD <- suppressWarnings(snow::parMM(cl,MTD1,D2))
      #MDT <- t(MTD)
      D3  <- diag(x=(rowSums(prot.similarity.final))^(-0.5))
      
      MTT1 <- suppressWarnings(snow::parMM(cl,D3,seq))
      MTT  <- suppressWarnings(snow::parMM(cl,MTT1,D3))
      
      D4  <- diag(x=(rowSums(drug.similarity.final))^(-0.5))
      MDD1  <- suppressWarnings(snow::parMM(cl,D4,csim))
      MDD  <- suppressWarnings(snow::parMM(cl,MDD1,D4))
      
      for (i in 1:nrow(seq)){
        for (j in 1:nrow(seq)){
          if (ATD[i]==0){
            
            MTT[i,j]<-prot.similarity.final[i,j]/St[i]
          }
          else{
            MTT[i,j]<-(0.8*(prot.similarity.final[i,j])/St[i])
          }
        }
      }
      
      for (i in 1:ncol(drugProt)){
        for (j in 1:ncol(drugProt)){
          if (ADT[i]==0){
            MDD[i,j]<-drug.similarity.final[i,j]/Sd[i]
          }
          else{
            MDD[i,j]<-(0.8*(drug.similarity.final[i,j])/Sd[i])
          }
        }
      }
      
      for (i in 1:ncol(drugProt)){
        for (j in 1:nrow(seq)){
          if (ADT[i]!=0){
            MDT[i,j]<- (0.2*A[i,j])/ADT[i]
          }
          else{
            MDT[i,j]<-0
          }
        }
      }
      for (i in 1:nrow(seq)){
        for (j in 1:ncol(drugProt)){
          if (ATD[i]!=0){
            MTD[i,j]<- (0.2*A[j,i])/ATD[i]
          }
          else{
            MTD[i,j]<-0
          }
        }
      }
      
      M1<-cbind(MTT,MTD)
      M2<-cbind(MDT,MDD)
      M <- rbind(M1,M2)
      M <- as.matrix(M)
      M[is.na(M)]<-0
      n =c(colnames(g1),rownames(g1))
      rownames(M) <- n
      colnames(M) <- n
      # Returning the final matrix 
      return(as.matrix(M))
      
    }
    
    
    if(normalise == "none"){
      MDD <- drug.similarity.final
      MTT <- prot.similarity.final
      MTD <- drugProt
      M1<-cbind(MTT,t(MTD))
      M2<-cbind(MTD,MDD)
      M <- rbind(M1,M2)
      M <- as.matrix(M)
      M[is.na(M)]<-0
      n =c(colnames(g1),rownames(g1))
      rownames(M) <- n
      colnames(M) <- n
      M <- colNorm(M)
      
      # Returning the final matrix 
      return(as.matrix(M))
    }
    on.exit(stopCluster(cl))
  }
  
  biNetwalk <- function(g1,s1,s2,normalise=c("laplace","none","chen"), dataSeed=NULL,restart=0.8,verbose=T,weight=FALSE) {
    
    startT <- Sys.time()
    
    if (!exists('s1') || !exists('s2')){
      stop("You must submit s1 and s2 matrices.\n")
    }
    
    if (class(g1) != "igraph"){
      stop("The function applies to 'igraph' object.\n")
    }
    
    if (!bipartite.mapping(g1)$res){
      stop("The function applies to bipartite graphs.\n")
    }
    
    if(verbose){
      now <- Sys.time()
      message(sprintf("First, get the adjacency matrix of the input graph (%s) ...", as.character(now)), appendLF=T)
    }
    if(is.null(restart) || is.na(restart) || restart<0 || restart>100){
      c <- 0.8
    }
    else{
      c <- restart
    }
    normalise <- match.arg(normalise)
    if (weight){
      if ("weight" %in% list.edge.attributes(g1)){
        adjM <- get.incidence(g1, attr="weight", names=T)
        if(verbose){
          message(sprintf("Notes: using weighted graph!"), appendLF=T)
        }
      }
    }else{
      adjM <- get.incidence(g1, attr=NULL, names=T)
      if(verbose){
        message(sprintf("Note: using unweighted graph!"), appendLF=T)
      }
    }
    adjM <- as.matrix(adjM)
    # get the transition matrix
    W = tMat(adjM,s1,s2,normalise=normalise)
    message(sprintf("got the transition matrix for RWR"))
    if(is.null(dataSeed)){
      
      M<-Matrix(adjM)
      M2<-0.99*M
      d<-Matrix(0.01*diag(nrow(s2)))
      P0matrix<-rBind(M2,d)
      
    }else{
      
      # part of the section for input file name
      drug.names <- as.character(unique(dataSeed$V2))
      P0matrix <- matrix(0,nrow=nrow(W),ncol=length(drug.names))
      
      for (i in 1:length(drug.names)){
        sub.fr <- dataSeed[dataSeed$V2==drug.names[i],]
        proteins <- as.character(sub.fr$V1)
        ind1 <- match(proteins, rownames(W))
        ind2 <- match(drug.names[i],rownames(W))
        ind <- append(ind1,ind2)
        nodes_mapped <- rownames(W)[ind[!is.na(ind)]]
        if(length(nodes_mapped)!=length(ind)){
          warning("The row names of input dataSeed do not contain all those in the input graph.\n")
        }
        
        P0matrix[ind[!is.na(ind)],i] <- 1 
      }
      P0matrix <- colNorm(P0matrix)
      
    }
    
    
    if (exists("W")){
      rmat <- rwr(W,P0matrix,r=c)
    } else{
      stop("Transition matrix couldnt be generated..")
    }
    
    if (!exists("rmat")){
      stop("Couldn't return the RWR matrix. \n")
    }else{
      if(verbose){
        now <- Sys.time()
        message(sprintf("Rescaling steady probability vector (%s) ...", as.character(now)), appendLF=T)
      }
      rmat[rmat < 1e-06] <- 0
      rmat <- rmat[1:nrow(adjM),]
      
      rmat <- colNorm(as.matrix(rmat))
      rownames(rmat)<- rownames(adjM)
      if(!is.null(dataSeed)){
        colnames(rmat)<- drug.names
        invisible(rmat)
      } else {
        colnames(rmat)<- colnames(adjM)
        invisible(rmat)
      }
      endT <- Sys.time()
      runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
      message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=T)    
      invisible(rmat)
      
    }
    
  }
#}