# -------------------------------------------------------------------------------
# 필요 함수 정의 (proposition coeff, Yhat)
# -------------------------------------------------------------------------------

# 각 모델의 optm값( 즉 model3의 경우 alpha0, alpha1 // model2의 경우 alpha0) 계산에 이용되는 값 계산
# E_ysq, E_lam_ysq, E_lamsq_ysq, E_lamsq_ybar, E_lamsq_ybarsq, E_lam_y, E_lamsq_y, E_lam_y_ybar, E_lamsq_y_ybar
prop_coeff <- function(lam_vec, a, tau, k){
  #NB-Gamma의 경우 아래 proposition 사용
  #E_ysq = var(u_lam_R) + mean(v_lam_R) + mean(u_lam_R)^2
  E_lam_vec    = mean(lam_vec)
  E_lam_vec_sq = mean(lam_vec^2)
  E_lam_vec_th = mean(lam_vec^3)
  E_lam_vec_for = mean(lam_vec^4)
  
  
  E_lamsq_ybar   = E_lam_vec_th
  E_lamsq_ybarsq = 1/tau * E_lam_vec_th + 1/(tau*k) * E_lam_vec_for * (1+1/a) + E_lam_vec_for * (1+1/a)
  E_lam_y        = E_lam_vec_sq
  E_lamsq_y      = E_lam_vec_th
  E_lam_y_ybar   = E_lam_vec_th * (1/a + 1)
  
  return (c(E_lamsq_ybar = E_lamsq_ybar, 
            E_lamsq_ybarsq = E_lamsq_ybarsq, E_lam_y = E_lam_y, E_lamsq_y = E_lamsq_y, E_lam_y_ybar = E_lam_y_ybar))
}


Yhat_ftn2 <- function( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0){
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
# -------------------------------------------------------------------------------
# 시뮬레이션aaa
# -------------------------------------------------------------------------------



# 1. HMSE의 차이를 기준으로 k선택
a = 0.8; N_sim = 1000; lambs = c(0.1, 0.2, 1.2); tau = 5; n_model = 3

k_seq = seq(1, 3, length = 100)
diff_HMSE <- rep(0, length(k_seq))


for(ii in (1:length(k_seq))){
  k = k_seq[ii]
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
  (optim_alpha(lam_vec, Y, par=c(result$par)))
  
  m3_alpha0 <- result$par
  
  
  # 각 모델별 예측값 추출 
  Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)
  
  # HMSE, MSE, MAE 계산
  HMSE <- rep(0, n_model)
  for(i in (1:n_model)){
    HMSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,i])
  }
  
  
  diff_HMSE[ii] = abs(HMSE[1] - HMSE[2])
  
}

plot(k_seq, diff_HMSE)
print(k_seq[which.max(diff_HMSE)])

# 결정된 parameter k를 이용해, 최종 HMSE, DIX 산출
a = 0.8; N_sim = 10000; lambs = c(0.1, 0.2, 1.2); tau = 5; n_model = 3
k = k_seq[which.max(diff_HMSE)] 



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
(optim_alpha(lam_vec, Y, par=c(result$par)))

m3_alpha0 <- result$par


# 각 모델별 예측값 추출 
Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)

# HMSE, MSE, MAE 계산
HMSE <- rep(0, n_model)
MSE <- rep(0, n_model)
for(i in (1:n_model)){
  HMSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'HMSE', Prem = Yhat[,i])
  MSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'MSE', Prem = Yhat[,i])
}




# 각 모델별 unbiasedness 확인
m2_DIX <- cal_UB(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 )


HMSE;
MSE;
m2_DIX


