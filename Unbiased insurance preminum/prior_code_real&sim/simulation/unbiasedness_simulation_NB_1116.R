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
  
  E_ysq          = E_lam_vec_sq + E_lam_vec + 1/k * E_lam_vec_sq * (1/a+1)
  E_lam_ysq      = E_lam_vec_sq + 1/k * E_lam_vec_th * (1+1/a) + E_lam_vec_sq
  E_lamsq_ysq    = E_lam_vec_th + 1/k * E_lam_vec_for * (1+1/a) + E_lam_vec_th
  E_lamsq_ybar   = E_lam_vec_th
  E_lamsq_ybarsq = 1/tau * E_lam_vec_th + 1/(2*k) * E_lam_vec_for * (1+1/a) + E_lam_vec_for * (1+1/a)
  E_lam_y        = E_lam_vec_sq
  E_lamsq_y      = E_lam_vec_th
  E_lam_y_ybar   = E_lam_vec_th * (1/a + 1)
  E_lamsq_y_ybar = E_lam_vec_for * (1+1/a)
  E_lam_ybarsq = 1/tau * E_lam_vec_sq + 1/(2*k) * E_lam_vec_th * (1+1/a) + E_lam_vec_th * (1+1/a)
  E_lamth_ybar = E_lam_vec_for
  E_lamth_y    = E_lam_vec_for
  
  return (c(E_ysq = E_ysq, E_lam_ysq = E_lam_ysq, E_lamsq_ysq = E_lamsq_ysq, E_lamsq_ybar = E_lamsq_ybar, 
            E_lamsq_ybarsq = E_lamsq_ybarsq, E_lam_y = E_lam_y, E_lamsq_y = E_lamsq_y, E_lam_y_ybar = E_lam_y_ybar,
            E_lamsq_y_ybar = E_lamsq_y_ybar, E_lam_ybarsq= E_lam_ybarsq, E_lamth_ybar = E_lamth_ybar, E_lamth_y = E_lamth_y))
}


Yhat_ftn <- function( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0){
  Yhat = data.frame(matrix(nrow = N_sim, ncol = n_model))
    
  # model 1,2,3 y_hat값 생성
  for(n in 1: N_sim){
    Yhat[n,1] = lam_vec[n] * (a+sum(Y[n,1:tau]))/(a+tau*lam_vec[n]) ### 수정필요 builmann으로
    Yhat[n,2] = lam_vec[n] * (m2_alpha0 + m2_alpha1 * sum(Y[n,1:tau])/tau)
    Yhat[n,3] = lam_vec[n]/(m3_alpha0 + (1- m3_alpha0) * lam_vec[n]) * ( m3_alpha0 + (1- m3_alpha0) * sum(Y[n,1:tau])/tau)
  }
  return( Yhat )
}


cal_mse <- function(n_model, y_true = Y[,6]){
  mse = data.frame(matrix(nrow = 1, ncol = n_model))
  for(j in 1: n_model){
    mse[1,j] = mean((Yhat[,j]-y_true)^2)
  }
  return( mse )
}


# 각 모델별 unbiasedness 계산

cal_UB <- function(lambs, n_model,Yhat, lam_vec){
  E_y_given_lam <- data.frame(matrix(nrow = length(lambs), ncol = n_model+1))
  
  total <- Yhat
  total$lam_vec <- lam_vec 
  for( i in 1: length(lambs)){
    E_y_given_lam[i,] = apply(total[total$lam_vec == lambs[i],],2,mean) # 각 모델의 lamb별 평균
  }
  
  total_unbiasedness_ch = apply((E_y_given_lam[,1:n_model]-lambs)^2/lambs,2,sum)
  return(total_unbiasedness_ch)
}


optim_alpha <- function(lam_vec, Y, par){
  alpha = par[1]
  tau <- ncol(Y[,-1])
  err <-  Y[,tau+1] - lam_vec/(alpha +(1-alpha) * lam_vec)*(alpha + (1-alpha)*apply(Y[,1:tau],1,mean))
  return( mean(err^2) )
}
# -------------------------------------------------------------------------------
# 시뮬레이션
# -------------------------------------------------------------------------------

a =0.1; N_sim = 10000; lambs = c(0.3, 0.6, 0.9); tau = 5; n_model = 3
k =3 # NB의 param

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

m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lambs^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))
m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))

# m3에 대한 alpha0 

result <- optim(par = c(0), 
      fn = optim_alpha , 
      lam_vec = lam_vec, Y = Y,
      control = list(fnscale = 1), method="L-BFGS-B")

# alpha 최적일때 err값
(optim_alpha(lam_vec, Y, par=c(result$par)))

m3_alpha0 <- result$par



######################
# model1 bayesian을 불만으로 수정해야함. 
######################

# 각 모델별 예측값 추출 
Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)

# MSE 계산
mse_model <- cal_mse(n_model, y_true = Y[,6])
print(mse_model)
# 각 모델별 unbiasedness 확인
total_unbiasedness_ch <- cal_UB(lambs,n_model,Yhat, lam_vec)

total_unbiasedness_ch

