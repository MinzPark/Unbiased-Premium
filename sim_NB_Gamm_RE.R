# -------------------------------------------------------------------------------
# Simulation for NB-Gamma random effect Model
# -------------------------------------------------------------------------------

## Calculate premium, DIX, HMSE of Builhmann, commercial, Generalized commercial

# create alpha0 plot folder
set.seed(1234)
path <- "C:/Users/user/Downloads/Unbiased_Premium" # Download path 

setwd(path)

#path_ftn <- paste0(path, "/Unbiased insurance preminum")
path_ftn <- path

ifelse(dir.exists(paste0(path_ftn,"\\m3_alpha0_plot")),F, dir.create(paste0(path_ftn,"\\m3_alpha0_plot") ) )


# import function to calculate each model

source(paste0(path_ftn,"/NB_Gamma_RE.R"))


#standard from real data lamb_kappa

lamb.c1 <- c( 0.05,0.33,0.61)

# generate each lambda case

lamb.c2 <- c(0.25,0.33,0.42)
lamb.c3 <- c(0.39,0.67,0.95)
lamb.c4 <- c(0.03,0.22,0.75)
total_lam <- matrix(rbind(lamb.c1, lamb.c2,lamb.c3, lamb.c4), ncol = length(lamb.c1))

# save results
df_tau = df <- data.frame()


for(scen_num in 1:(dim(total_lam)[1])){
  startTime <- Sys.time() 
  for( tau in c(5,10)){
    
    HMSE <- data.frame(sprem = double(),cprem = double(),gprem = double())
    
  for(iter in 1:100){

    
  
    # save param
    l <- list( N_sim = NA, tau = tau, lam_ = NA, a = NA, k = NA, 
               s_hmse.mean  = NA, s_hmse.std  = NA, 
               c_hmse.mean  = NA, c_hmse.std  = NA, 
               gc_hmse.mean = NA, gc_hmse.std = NA,
               c_dix.mean = NA)
    
    lambs = total_lam[scen_num,];
  
    a = 0.8; N_sim = 5000; n_model = 3; k = 3;
    
    # create lamb.vec with each prob
    lam.vec_wo_last <- rep(lambs[-length(lambs)], rep(N_sim/length(lambs), length(lambs[-length(lambs)])))
    lam_vec <- c(rep(lambs[length(lambs)], rep(N_sim - length(lam.vec_wo_last),1)), lam.vec_wo_last)

    R = rgamma(N_sim, shape = a, rate = a)
    
    # y1,...,ytau,ytau+1 빈공간 생성
    Y = data.frame(matrix(nrow = N_sim, ncol = tau+1))
    for(i in 1:N_sim){
      for(j in 1: (tau+1) ){
        Y[i,j] = rnbinom(1, size = k, prob = k/(k+lam_vec[i]*R[i]) )  
      }
    }
    
    # calculate (alpha0, alpha1) for CPrem
    P <- prop_coeff(lam_vec, a, tau, k)
    
    m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lam_vec^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
    m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
    
    
    
    # optimize alpha0 for GPrem
    result <- optim(par = c(0), 
                    fn = optim_alpha , 
                    lam_vec = lam_vec, Y = Y,
                    control = list(fnscale = 1), method="BFGS")
    
    # optimal alpha0 for GPrem
    
    m3_alpha0 <- result$par

      
    # plot alpha0, notice) draw just one plot whether it repeats or not
    start_time <- Sys.time()[3]
    if(iter == 1){
      alpha0_plot(path_ftn, m3_alpha0, m2_alpha0, m2_alpha1, lam_vec, Y, R)
      print('break?')
    }

    end_time <- Sys.time()[3]
    execution_time <- end_time - start_time
    if(is.na(execution_time) ==FALSE ) {cat("Execution time:", execution_time, "seconds\n");}


    ##########################################################
    # final calculated premium, DIX, HMSE for each models
    ##########################################################
    
    # estimation of yhat for each model
    Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)
    

    
    for(i in (1:n_model)){
      # HMSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,i])
      HMSE[iter, i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,i])
    }  
    
  }
    
    # take round each value
    HMSE <- round(HMSE,4)
    HMSE;
    
    m2_DIX <- cal_UB(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 )
    m2_DIX <- round(m2_DIX,4)
    
    
    # to save each case 
    l['N_sim'] <- N_sim; l['tau'] <- tau; l['lam_'] <- paste(lambs, collapse = ","); l['a'] <- a; l['k'] <- k;
    
    l['s_hmse.mean']  <- mean(HMSE[,1]); l['s_hmse.std']  <- sqrt(1/100*var(HMSE[,1]))
    l['c_hmse.mean']  <- mean(HMSE[,2]); l['c_hmse.std']  <- sqrt(1/100*var(HMSE[,2]))
    l['gc_hmse.mean'] <- mean(HMSE[,3]); l['gc_hmse.std'] <- sqrt(1/100*var(HMSE[,3]))
    l['c_dix.mean']   <- m2_DIX[[1]]
    
    df_tau <- rbind(df_tau, data.frame(l))
    print(df_tau)
  }
  df <- rbind(df, df_tau)
  endTime <- Sys.time()
  cat('running time each scenario:',endTime - startTime)
}


ifelse(dir.exists(paste0(path_ftn,"\\result")),F, dir.create(paste0(path_ftn,"\\result")) )

new_result_sim_nb <- paste0(path_ftn,"\\result\\", format(Sys.Date(),"%m%d"),"_Sim")

if (!file.exists(new_result_sim_nb)){
  dir.create(new_result_sim_nb)
}


filename <- paste0("/",format(Sys.Date(),"%m%d"), "_sim_NB_gamm.csv")
write.csv(df_tau, file = paste0(new_result_sim_nb,"/sim_NB_gamm.csv"))

