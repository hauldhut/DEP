nmf <- function (nsclc2){
  
  system.time(nsclc2.nmf <- nnmf(nsclc2, 2));
  nsclc2.hat.nmf <- with(nsclc2.nmf, W %*% H);
  return <- nsclc2.hat.nmf
  
}