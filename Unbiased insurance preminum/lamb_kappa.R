
combine_kappa_level <- function(data.train){
  M <- c()
  for(i in 1:3){
    M <- rbind(M, cbind( diag(1, nrow = 6, ncol = 6), i-1))
  }  
  M<- data.frame(M)
  colnames(M) <- c("TypeMisc","TypeCity","TypeCounty", "TypeSchool","TypeTown","TypeVillage","col.Cov.Idx")
  lamb_kappa <- rep(0,dim(data.train)[1])
  
  for(level in 1:nrow(M)){
    lamb_broad <- data.frame((matrix(rep(M[level,], each = nrow(data.train)),nrow = nrow(data.train))))
    lamb_level.id = apply((data.train[,3:9] == lamb_broad),1,sum) == 7
    lamb_kappa[lamb_level.id] = level  
  }
  data.train$lamb_kappa <- lamb_kappa
  return(data.train)
}




