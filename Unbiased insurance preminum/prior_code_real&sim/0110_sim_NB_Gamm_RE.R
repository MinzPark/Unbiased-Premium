# -------------------------------------------------------------------------------
# Simulation for NB-Gamma random effect Model
# -------------------------------------------------------------------------------

## Calculate premium, DIX, HMSE of Builhmann, commercial, Generalized commercial

# create alpha0 plot folder
set.seed(123)
path <- "C:/Users/none/Desktop/대학원/논문/Unbiased insurance preminum"

ifelse(dir.exists(paste0(path,"\\m3_alpha0_plot")),F, dir.create(paste0(path,"\\m3_alpha0_plot") ) )


# import function to calculate each model

source(paste0(path,"/NB_Gamma_RE.R"))

#source("C:/Users/parkminji/Downloads/민지/Unbiased insurance preminum/NB_Gamma_RE.R")


# lam_vec case
total_lam <- matrix(c(0.1,0.4,0.6,0.1,0.5,0.5,1.0,0.2,0.9,0.6,1.4,1.2),nrow = 4, ncol =3)

# save results
df_tau = df <- data.frame()
for( tau in c(5,10)){
  for(i in 1:(dim(total_lam)[1])){
    # save param
    l <- list( N_sim = NA, tau = tau, lam_ = NA, a = NA, k = NA, s_hmse = NA, c_hmse = NA, gc_hmse = NA, c_dix = NA)

    
    # #lambs = c(0.1, 0.2, 1.2); 
    lambs = total_lam[i,];
    a = 0.8; N_sim = 10000; n_model = 3; k = 8;
    
    lam_vec = c(rep(lambs[1], floor(N_sim/3)),rep(lambs[2], floor(N_sim/3)),rep(lambs[3], N_sim-2*floor(N_sim/3)))
    R = rgamma(N_sim, shape = a, rate = a)
    
    # y1,...,ytau,ytau+1 빈공간 생성
    Y = data.frame(matrix(nrow = N_sim, ncol = tau+1))
    for(i in 1:N_sim){
      for(j in 1: (tau+1) ){
        Y[i,j] = rnbinom(1, size = k, prob = k/(k+lam_vec[i]*R[i]) )  
      }
    }
    
    # 각 모델에 대한 mod2의 알파0,알파1과 mod3의 알파0 계산
    P <- prop_coeff(lam_vec, a, tau, k)
    
    m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lam_vec^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
    m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
    
    
    
    # m3에 대한 alpha0 
    result <- optim(par = c(0), 
                    fn = optim_alpha , 
                    lam_vec = lam_vec, Y = Y,
                    control = list(fnscale = 1), method="BFGS")
    
    # alpha 최적일때 err값
    # (optim_alpha(lam_vec, Y, par=c(result$par)))
    
    m3_alpha0 <- result$par
    
    
    ##########################################################
    # plotting estimation of alpha0 for generalized comm
    ##########################################################
    m3_alpha0_seq <- seq(0,1,length = 30)
    HMSE_alpha0 <- rep(0, length(m3_alpha0_seq))
    
    for(i in 1: length(m3_alpha0_seq)){
      # 각 모델별 예측값 추출
      Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0 = m3_alpha0_seq[i])
      # calculate HMSE of Gprem to check which is the best alpha0
      HMSE_alpha0[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,3])
    }
    
    # best HMSE for alpha0
    Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0 = m3_alpha0)
    best_HMSE <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,3])
    
    
    # save png for each alpha0
    
    new_folder_name <- paste0(path,"\\m3_alpha0_plot\\", format(Sys.Date(),"%m%d"),"_NB_Gamm")
    if (!file.exists(new_folder_name)){
      dir.create(new_folder_name)
    }
    
    setwd(new_folder_name)
    png(filename=paste0(paste(c('tau',tau,'lambs',lambs), collapse = "_"),'.png'),
        width=600,height=600)
    x_axis_tick=seq(0,1,length=5)
    y_axis_tick=round(seq(min(HMSE_alpha0),max(HMSE_alpha0),length=5),2)
    plot(m3_alpha0_seq, HMSE_alpha0,type = 'l', lwd = 1.5, xlab = expression(alpha[0]), ylab = "HMSE",xlim = c(0,1), ylim = c(min(HMSE_alpha0),max(HMSE_alpha0)),axes=FALSE, cex.lab=1.3)
    axis(side=1,at=x_axis_tick)
    axis(side=2,at=y_axis_tick)
    points(m3_alpha0, best_HMSE,pch = 17, cex = 1.5)
    
    
    ##########################################################
    # final calculated premium, DIX, HMSE for each models
    ##########################################################
    
    # 각 모델별 예측값 추출 
    Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)
    
    # HMSE 계산
    HMSE <- rep(0, n_model)
    
    for(i in (1:n_model)){
      HMSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,i])
    }  
    
    
    # 각 모델별 unbiasedness 확인
    m2_DIX <- cal_UB(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 )
    
    HMSE <- round(HMSE,4)
    HMSE;
    m2_DIX
    
    
    # to save each case 
    l['N_sim'] <- N_sim; l['tau'] <- tau; l['lam_'] <- paste(lambs, collapse = ","); l['a'] <- a; l['k'] <- k;
    
    l['s_hmse'] <- HMSE[1]
    l['c_hmse'] <- HMSE[2]
    l['gc_hmse'] <- HMSE[3]
    l['c_dix'] <- m2_DIX[[1]]
    
    df_tau <- rbind(df_tau, data.frame(l))
  }
  df <- rbind(df, df_tau)
  print(df)
}

filename <- paste0("/",format(Sys.Date(),"%m%d"), "_sim_NB_gamm.csv")
write.csv(df_tau, file = paste0(path,"/sim_NB_gamm.csv"))

