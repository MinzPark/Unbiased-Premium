# -------------------------------------------------------------------------------
# 필요 함수 정의 (proposition coeff, Yhat)
# -------------------------------------------------------------------------------

# 각 모델의 optm값( 즉 model3의 경우 alpha0, alpha1 // model2의 경우 alpha0) 계산에 이용되는 값 계산
# E_ysq, E_lam_ysq, E_lamsq_ysq, E_lamsq_ybar, E_lamsq_ybarsq, E_lam_y, E_lamsq_y, E_lam_y_ybar, E_lamsq_y_ybar
prop_coeff <- function(lam_vec, a, tau){
  #E_ysq = var(u_lam_R) + mean(v_lam_R) + mean(u_lam_R)^2
  E_lam_vec    = mean(lam_vec)
  E_lam_vec_sq = mean(lam_vec^2)
  E_lam_vec_th = mean(lam_vec^3)
  E_lam_vec_for = mean(lam_vec^4)
  
  
  E_lamsq_ybar   = E_lam_vec_th
  E_lamsq_ybarsq = 1/tau * E_lam_vec_th +E_lam_vec_for *(1+1/a)
  E_lam_y        = E_lam_vec_sq
  E_lam_y_ybar   = E_lam_vec_th * (1/a + 1)
  
  return (c(E_lamsq_ybar = E_lamsq_ybar, 
            E_lamsq_ybarsq = E_lamsq_ybarsq, E_lam_y = E_lam_y, E_lam_y_ybar = E_lam_y_ybar))
}


Yhat_ftn <- function( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0){
  #Poisson-Gamma
  Yhat = data.frame(matrix(nrow = N_sim, ncol = n_model))
  
  # model 1,2,3 y_hat값 생성
  for(n in 1: N_sim){
    Yhat[n,1] = lam_vec[n] * (a+sum(Y[n,1:tau]))/(a+tau*lam_vec[n]) ### 수정필요 builmann으로
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
  
  return( prediction )
}


# 각 모델별 unbiasedness 계산

cal_UB <- function(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 ){
  var_exp_prem = mean( ( alpha0 + alpha1 * lam_vec)^2) - (mean( alpha0 + alpha1 * lam_vec))^2
  var_prem     = alpha1^2 * mean( 1/tau * lam_vec + lam_vec^2 *(1+1/a) ) + var_exp_prem
  return(var_exp_prem / var_prem)
}


optim_alpha <- function(lam_vec, Y, par){
  alpha = par[1]
  err <-  Y[,tau+1] - lam_vec/(alpha +(1-alpha) * lam_vec)*(alpha + (1-alpha)*apply(Y[,1:tau],1,mean))
  return( mean(err^2) )
}