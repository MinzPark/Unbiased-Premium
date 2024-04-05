# -------------------------------------------------------------------------------
# 필요 함수 정의 (proposition coeff, Yhat)
# -------------------------------------------------------------------------------

# 각 모델의 optm값( 즉 model3의 경우 alpha0, alpha1 // model2의 경우 alpha0) 계산에 이용되는 값 계산
# E_ysq, E_lam_ysq, E_lamsq_ysq, E_lamsq_ybar, E_lamsq_ybarsq, E_lam_y, E_lamsq_y, E_lam_y_ybar, E_lamsq_y_ybar
prop_coeff <- function(lam_vec, a, tau, k){
  #E_ysq = var(u_lam_R) + mean(v_lam_R) + mean(u_lam_R)^2
  E_lam_vec    = mean(lam_vec)
  E_lam_vec_sq = mean(lam_vec^2)
  E_lam_vec_th = mean(lam_vec^3)
  E_lam_vec_for = mean(lam_vec^4)
  
  
  E_lamsq_ybar   = E_lam_vec_th
  E_lamsq_ybarsq = 1/tau * E_lam_vec_th + 1/(tau*k) * E_lam_vec_for * (1+1/a) + E_lam_vec_for * (1+1/a)
  E_lam_y        = E_lam_vec_sq
  E_lam_y_ybar   = E_lam_vec_th * (1/a + 1)
  
  return (c(E_lamsq_ybar = E_lamsq_ybar, 
            E_lamsq_ybarsq = E_lamsq_ybarsq, E_lam_y = E_lam_y, E_lam_y_ybar = E_lam_y_ybar))
}


Yhat_ftn <- function( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0){
  #NB-Gamma
  Yhat = data.frame(matrix(nrow = N_sim, ncol = n_model))
  
  # model 1,2,3 y_hat값 생성
  for(n in 1: N_sim){
    u_ = lam_vec[n]  
    v_ = lam_vec[n]+lam_vec[n]^2/k*(1/a+1)
    a_ = lam_vec[n]^2/a
    z_ = tau*a_/(v_+tau*a_)
    Yhat[n,1] = z_ * apply(Y[n,1:tau],1,mean) + (1-z_)*u_
    Yhat[n,2] = lam_vec[n] * (m2_alpha0 + m2_alpha1 * sum(Y[n,1:tau])/tau)
    Yhat[n,3] = lam_vec[n]/(m3_alpha0 + (1- m3_alpha0) * lam_vec[n]) * ( m3_alpha0 + (1- m3_alpha0) * sum(Y[n,1:tau])/tau)
  }
  return( Yhat )
}



cal_pred <- function(lam_vec, Y = Y[,tau+1], Theta = R, type,Prem = Yhat){
  
  if(type == "HMSE"){
    prediction = mean( (lam_vec * Theta - Prem)^2)
  }
  else if(type == "MSE"){
    prediction = mean( (     Y          - Prem)^2)
  }
  else if(type == "MAE"){
    prediction = mean(abs(     Y          - Prem))
  }
  
  return( prediction )
}


# 각 모델별 unbiasedness 계산

cal_UB <- function(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 ){
  var_exp_prem = mean( ( alpha0 + alpha1 * lam_vec)^2) - (mean( alpha0 + alpha1 * lam_vec))^2
  var_prem     = alpha1^2 * mean( 1/tau * (lam_vec + lam_vec^2/k *(1+1/a)) + lam_vec^2* 1/a) + var_exp_prem
  return(var_exp_prem / var_prem)
}


optim_alpha <- function(lam_vec, Y, par){
  alpha = par[1]
  err <-  Y[,tau+1] - lam_vec/(alpha +(1-alpha) * lam_vec)*(alpha + (1-alpha)*apply(Y[,1:tau],1,mean))
  return( mean(err^2) )
}






alpha0_plot <- function(path_ftn, m3_alpha0, m2_alpha0, m2_alpha1, lam_vec, Y, R){
  
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
  
  new_folder_name <- paste0(path_ftn,"\\m3_alpha0_plot\\", format(Sys.Date(),"%m%d"),"_NB_Gamm")
  
  # version with triangle
  
  new_folder_name_tri <- paste0(new_folder_name,"_NB_Gamm_tri")
  if (!file.exists(new_folder_name_tri)){
    dir.create(new_folder_name_tri)
  }
  
  setwd(new_folder_name_tri)
  
  png(filename=paste(c('tau', tau,'alpha0.sim.nb',scen_num,'png'),collapse = '.'),width=2200,height=2000)
  
  #dev.new()
  
  x_axis_tick=seq(0,1,length=3)
  y_axis_tick=round(seq(min(HMSE_alpha0),max(HMSE_alpha0),length=3),2)
  par(oma = c(1,1,1,1))
  par(mar=c(25,20,20,5)) # mar=c(아래,왼쪽,위,오른쪽)
  #par(las = 1)
  oldp <- par(mgp=c(3, 6.5, 0.5)) # mgp = c(메인 타이틀, 눈금 표시 레이블, 눈금 표시)
  
  
  plot(m3_alpha0_seq,HMSE_alpha0,type = 'l', xlab =''
       , ylab = "",xlim = c(0,1),
       ylim = c(min(HMSE_alpha0),max(HMSE_alpha0))
       ,axes=FALSE, xaxt = 'n',yaxt = 'n', lwd = 5)
  
  axis(side=1,at=x_axis_tick, cex.axis = 10, lwd = 4) #, mgp=c(0, 0,-100)
  axis(side=2,at=y_axis_tick, cex.axis = 10, lwd = 4,mgp=c(0, 2,0))
  
  mtext(expression(alpha[0]), side = 1, line= 20, cex = 12)
  mtext("HMSE", side = 2, line= 12, cex = 12)
  
  points(m3_alpha0, best_HMSE,pch = 17, cex = 8)
  
  
  dev.off()
}