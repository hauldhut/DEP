
function.Impute_kfold <- function(x,rankImp,lambdaImp){
  library(softImpute)
  # lay so hang so cot
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  np=nrow*ncol
  
  # count zero/nonzero in the matrix
  zero <-length(which(x==0)) 
  nonzero <-length(which(x!=0))
  
  xna<-x
  # xac dinh vetor chi so (thu tu) phan tu bang Zero trong ma tran
  ix<-which(xna==0)
  #lay cac chi so phan tu  =0
  #gan cac phan tu 0 => NA
  xna[ix]=NA
  
  ###uses sparse matrix method for matrices of class "Incomplete"
  xnaC<-as(xna,"Incomplete")
  fit<-softImpute(xnaC,rank=rankImp,lambda=lambdaImp,type="svd")
  xc=biScale(fit,col.scale=FALSE,row.scale=FALSE) 
  ximp<-complete(xc,fit)
  return(ximp)
  
}