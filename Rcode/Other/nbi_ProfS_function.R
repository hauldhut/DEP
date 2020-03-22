nbi_ProbS <- function(Y, lambdaY){
  n <- nrow(Y)
  m <- ncol(Y)
  # remove col, row which have sum ==0
    idxZeroCol <- which(colSums(Y) == 0)
    if (length(idxZeroCol)!=0){
      Y <- Y[, -idxZeroCol] 
      
    }
    i=length(idxZeroCol)
    idxZeroRow <- which(rowSums(Y) == 0)
    if (length(idxZeroRow)!=0){
      Y <- Y[-idxZeroRow, ]
      
    }
  
  Ky <- diag(1/colSums(Y)) # tao matran duong cheo = length(sumcol)
  
  Ky[is.infinite(Ky) | is.na(Ky)] <- 0 # is.infinite: tra ve 1 vector chieu dai bang chieu cua Ky
  #  cho biet phan tu nao la finite - huu han (not missing) hoac vo han
  
  kx <- rowSums(Y)
  Nx <- 1/(matrix(kx, nrow=n, ncol=n, byrow=TRUE)^(lambda) * 
             
             matrix(kx, nrow=n, ncol=n, byrow=FALSE)^(1-lambda))
  
  Nx[is.infinite(Nx) | is.na(Nx)] <- 0 
  kx[is.infinite(kx) | is.na(kx)] <- 0 
  
  # this is the first resource pass from X to Y nodes where the degree information is passed . Look for dimension of adjacency matrix
  
  W <- t(Y %*% Ky) 
  #print("first resource pass from X to Y nodes")
  #print.table(W)
  # Again here it returns back the resource to the X nodes
  
  W <- Y %*% W 
  #print("Again here it returns back the resource to the X nodes")
  #print.table(W)
  # This is the final scaling which needs to be done you cannot do mulliplication here you are just scaling
  
  W <- Nx * W 
  #print("W <- Nx * W ")
  #print.table(W)
  rownames(W) <- rownames(Y)
  colnames(W) <- rownames(Y)
  R <- W %*% Y
  #print("ket qua R")
  
  #print.table(R)
  return(R)
 }
