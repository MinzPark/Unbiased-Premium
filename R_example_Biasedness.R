# import function to calculate each model

path <- "C:/Users/user/Downloads/Unbiased-Premium-main/Unbiased-Premium-main" # Download path 
path_ftn <- path
source(paste0(path_ftn,"/Pois_Gamma_RE.R"))

total_lam <- matrix(c(3,1,4,6),nrow = 2, ncol =2)

nsim= 1000

tau=c(3,10)
a = c(2)



## data
ls <- list()

ls[[1]] <- tau 
ls[[2]] <- a

a_and_tau <- data.frame(as.matrix(expand.grid(ls[[2]],ls[[1]]))); 
colnames(a_and_tau) <-c('a','tau')
total_df <- data.frame(c())
result_df <- data.frame(c())

for(lam in 1:dim(total_lam)[1]){
  for( scen in 1:dim(a_and_tau)[1]){
    
    lambs = total_lam[lam,]
    lam_vec = c(rep(lambs[1], floor(nsim/2)), rep(lambs[2], nsim - floor(nsim/2)))
    a   = a_and_tau[scen, ]$a
    tau = a_and_tau[scen, ]$tau
    
    l <- list(a = a, lam_ = paste0('(',paste0(lambs, collapse = ','),')'), tau = tau, 
              CPost.bayes.lamb1_mean = NA,
              CPost.bayes.lamb1_std  = NA,
              CPost.bayes.lamb2_mean = NA,
              CPost.bayes.lamb2_std  = NA,
              
              CPost.cred.lamb1_mean  = NA,
              CPost.cred.lamb1_std   = NA,
              CPost.cred.lamb2_mean  = NA,
              CPost.cred.lamb2_std   = NA)
    result <- list(a = a, lam_ = paste0('(',paste0(lambs, collapse = ','),')'), tau = tau, 
                   CPost.bayes.dix        = NA,
                   CPost.bayes.dix_std    = NA,
                   CPost.bayes.mse        = NA,
                   CPost.bayes.mse_std    = NA,
                   
                   CPost.cred.dix        = NA,
                   CPost.cred.dix_std    = NA,
                   CPost.cred.mse        = NA,
                   CPost.cred.mse_std    = NA)
    
    CPost.bayes_mean = data.frame(lamb1 = double(),lamb2 = double())
    CPost.cred_mean  = data.frame(lamb1 = double(),lamb2 = double())
    CPost.bayes_dix  = data.frame(dix = double())
    CPost.cred_dix   = data.frame(dix = double())
    CPost.bayes_mse  = data.frame(mse = double())
    CPost.cred_mse   = data.frame(mse = double())
    for(iter in 1:5){
      theta = rep(rgamma(nsim, shape=a, rate=a), rep((tau+1), nsim)) # rgamma의 vector가 3번씩 nsim만큼 생성됨 0.1,0.1,0.1,0.2,0.2,0.2,...
      lambda <- rep(lam_vec, rep((tau+1),nsim))
      Ns=matrix(rpois(nsim*(tau+1), lambda=theta* lambda ), ncol=(tau+1), byrow=TRUE)
      sum_n = rowSums(Ns[,1:tau])
      
      id_lam.1 <- lam_vec == lambs[1]
      id_lam.2 <- lam_vec == lambs[2]
      
      c1 = lambs[1]^(2+sum_n)/((a+lambs[1]*tau)^(a+sum_n)) 
      c2 = lambs[2]^(2+sum_n)/((a+lambs[2]*tau)^(a+sum_n)) 
      CPost.bayes = (a+sum_n)*( c1/(a+lambs[1]*tau)/(c1+c2) + c2/(a+lambs[2]*tau)/(c1+c2) )

      CPost.bayes_mean[iter, "lamb1"] <- mean(CPost.bayes[id_lam.1], na.rm = TRUE)
      CPost.bayes_mean[iter, "lamb2"] <- mean(CPost.bayes[id_lam.2], na.rm = TRUE)
      
      # dix
      CPost.bayes_dix[iter, 'dix'] <- apply(CPost.bayes_mean[iter,],1,function(x){var(x)}) /var(CPost.bayes,na.rm = TRUE)
      
      # mse
      threshold <- quantile(CPost.bayes, na.rm = TRUE, 0.98) # remove outlier due to param given the assumped distn 
      threshold.idx <-which(CPost.bayes < threshold)
      CPost.bayes_mse[iter, 'mse'] <- mean((Ns[threshold.idx,tau+1] -  lam_vec[threshold.idx] * CPost.bayes[threshold.idx])^2, na.rm = TRUE) 
      
      #CPost cred
      P <- prop_coeff(lam_vec, a, tau)
      
      m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lam_vec^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
      m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
      
      CPost.cred <- (m2_alpha0 + m2_alpha1 * sum_n/tau)
      CPost.cred_mean[iter, "lamb1"] <- mean(CPost.cred[id_lam.1])
      CPost.cred_mean[iter, "lamb2"] <- mean(CPost.cred[id_lam.2])
      # dix
      CPost.cred_dix[iter, 'dix'] <- apply(CPost.cred_mean[iter,],1,function(x){var(x)}) /var(CPost.cred,na.rm = TRUE)

      
      # mse
      CPost.cred_mse[iter, 'mse'] <- mean((Ns[,tau+1] -  lam_vec * CPost.cred)^2, na.rm = TRUE) 
    }
    # save simulation for each scenario
    
    l['CPost.bayes.lamb1_mean'] <- mean(CPost.bayes_mean$lamb1)
    l['CPost.bayes.lamb1_std']  <- sqrt(1/100*var(CPost.bayes_mean$lamb1))
    l['CPost.bayes.lamb2_mean'] <- mean(CPost.bayes_mean$lamb2)
    l['CPost.bayes.lamb2_std']  <- sqrt(1/100*var(CPost.bayes_mean$lamb2))
    
    l['CPost.cred.lamb1_mean'] <- mean(CPost.cred_mean$lamb1)
    l['CPost.cred.lamb1_std']  <- sqrt(1/100*var(CPost.cred_mean$lamb1))
    l['CPost.cred.lamb2_mean'] <- mean(CPost.cred_mean$lamb2)
    l['CPost.cred.lamb2_std']  <- sqrt(1/100*var(CPost.cred_mean$lamb2))
    
    result['CPost.bayes.dix']        <- mean(CPost.bayes_dix$dix)
    result['CPost.bayes.dix_std']    <- sqrt(1/100*var(CPost.bayes_dix$dix))
    result['CPost.bayes.mse']        <- mean(CPost.bayes_mse$mse)
    result['CPost.bayes.mse_std']    <- sqrt(1/100*var(CPost.bayes_mse$mse))
    
    result['CPost.cred.dix']        <- mean(CPost.cred_dix$dix)
    result['CPost.cred.dix_std']    <- sqrt(1/100*var(CPost.cred_dix$dix))
    result['CPost.cred.mse']        <- mean(CPost.cred_mse$mse)
    result['CPost.cred.mse_std']    <- sqrt(1/100*var(CPost.cred_mse$mse))
    
    total_df <- rbind(total_df, data.frame(l))  
    result_df <- rbind(result_df, data.frame(result))  
  }
}

# save the results of ex.9

ifelse(dir.exists(paste0(path_ftn,"/result")), F, dir.create(paste0(path_ftn,"/result")))
new_folder_ex9 <- paste0(path_ftn,"/result/", format(Sys.Date(),"%m%d"), "_ex6")
ifelse(dir.exists(new_folder_ex9), F, dir.create(new_folder_ex9))

write.csv(total_df, file = paste0(new_folder_ex9, "/result_ex6.csv"))
write.csv(result_df, file = paste0(new_folder_ex9, "/result_perf_ex6.csv"))




