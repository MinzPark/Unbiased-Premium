# import function to calculate each model
path <- "C:/Users/none/Desktop/대학원/논문" # Download path 
setwd(path)

path_ftn <- paste0(path,'/Unbiased insurance preminum')
source(paste0(path_ftn,"/Pois_Gamma_RE.R"))

total_lam <- matrix(c(3,1,4,6),nrow = 2, ncol =2)


nsim= 5

tau=c(3,10)
a = c(2,10)



## data
ls <- list()

ls[[1]] <- tau 
ls[[2]] <- a

a_and_tau <- data.frame(as.matrix(expand.grid(ls[[2]],ls[[1]]))); 
colnames(a_and_tau) <-c('a','tau')
total_df <- data.frame(c())
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
    
    CPost.bayes_mean = data.frame(lamb1 = double(),lamb2 = double())
    CPost.cred_mean  = data.frame(lamb1 = double(),lamb2 = double())
    
    for(iter in 1:100){
      theta = rep(rgamma(nsim, shape=a, rate=a), rep(tau, nsim)) # rgamma의 vector가 3번씩 nsim만큼 생성됨 0.1,0.1,0.1,0.2,0.2,0.2,...
      lambda <- rep(lam_vec, rep(tau,nsim))
      Ns=matrix(rpois(nsim*tau, lambda=theta* lambda ), ncol=tau, byrow=TRUE)
      sum_n = rowSums(Ns)
      
      id_lam.1 <- lam_vec == lambs[1]
      id_lam.2 <- lam_vec == lambs[2]
      
      c1 = lambs[1]^(2+sum_n)/((a+lambs[1]*tau)^(a+sum_n))
      c2 = lambs[2]^(2+sum_n)/((a+lambs[2]*tau)^(a+sum_n))
      CPost.bayes = (a+sum_n)/(a+lambs[1]*tau)*c1/(c1+c2)+
        (a+sum_n)/(a+lambs[2]*tau)*c2/(c1+c2)
      CPost.bayes_mean[iter, "lamb1"] <- mean(CPost.bayes[id_lam.1])
      CPost.bayes_mean[iter, "lamb2"] <- mean(CPost.bayes[id_lam.2])
      
      #CPost cred
      P <- prop_coeff(lam_vec, a, tau)
      
      m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lam_vec^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
      m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
      
      CPost.cred <- lam_vec * (m2_alpha0 + m2_alpha1 * sum_n/tau)
      CPost.cred_mean[iter, "lamb1"] <- mean(CPost.cred[id_lam.1])
      CPost.cred_mean[iter, "lamb2"] <- mean(CPost.cred[id_lam.2])
      
    }
    # save simulation for each scenario
    
    l['CPost.bayes.lamb1_mean'] <- round(mean(CPost.bayes_mean$lamb1),3)
    l['CPost.bayes.lamb1_std']  <- round(sqrt(1/100*var(CPost.bayes_mean$lamb1)),3)
    l['CPost.bayes.lamb2_mean'] <- round(mean(CPost.bayes_mean$lamb2),3)
    l['CPost.bayes.lamb2_std']  <- round(sqrt(1/100*var(CPost.bayes_mean$lamb2)),3)
    
    l['CPost.cred.lamb1_mean'] <- round(mean(CPost.cred_mean$lamb1),3)
    l['CPost.cred.lamb1_std']  <- round(sqrt(1/100*var(CPost.cred_mean$lamb1)),3)
    l['CPost.cred.lamb2_mean'] <- round(mean(CPost.cred_mean$lamb2),3)
    l['CPost.cred.lamb2_std']  <- round(sqrt(1/100*var(CPost.cred_mean$lamb2)),3)
    
    total_df <- rbind(total_df, data.frame(l))  
  }
}

# save the results of ex.9

ifelse(dir.exists(paste0(path_ftn,"/result")), F, dir.create(paste0(path_ftn,"/result")))
new_folder_ex9 <- paste0(path_ftn,"/result/", format(Sys.Date(),"%m%d"), "_ex9")
ifelse(dir.exists(new_folder_ex9), F, dir.create(new_folder_ex9))

write.csv(total_df, file = paste0(new_folder_ex9, "/result_ex9.csv"))




