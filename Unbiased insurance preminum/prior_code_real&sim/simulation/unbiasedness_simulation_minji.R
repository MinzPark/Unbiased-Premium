options(digits=22)


a =0.1; N_sim = 1000; lambs = c(0.8,1, 1.2); tau = 5; n_model = 3

lam_vec = c(rep(lambs[1], floor(N_sim/3)),rep(lambs[2], floor(N_sim/3)),rep(lambs[3], N_sim-2*floor(N_sim/3)))
R = rgamma(N_sim, shape = a, rate = a)

# y1,...,ytau,ytau+1 빈공간 생성
Y = data.frame(matrix(nrow = N_sim, ncol = tau+1))
for(i in 1:N_sim){
  for(j in 1: (tau+1) ){
    Y[i,j] = rpois(1, lam_vec[i]*R[i])  
  }
}


# 각 모델의 optm값( 즉 model3의 경우 alpha0, alpha1 // model2의 경우 alpha0) 계산에 이용되는 값 계산
# E_ysq, E_lam_ysq, E_lamsq_ysq, E_lamsq_ybar, E_lamsq_ybarsq, E_lam_y, E_lamsq_y, E_lam_y_ybar, E_lamsq_y_ybar
prop_coeff <- function(lam_vec, a, tau){
  #E_ysq = var(u_lam_R) + mean(v_lam_R) + mean(u_lam_R)^2
  E_lam_vec    = mean(lam_vec)
  E_lam_vec_sq = mean(lam_vec^2)
  E_lam_vec_th = mean(lam_vec^3)
  E_lam_vec_for = mean(lam_vec^4)
  
  E_ysq          = E_lam_vec_sq*(1/a+1) + E_lam_vec
  E_lam_ysq      = 2*E_lam_vec_sq
  E_lamsq_ysq    = 2*E_lam_vec_th
  E_lamsq_ybar   = E_lam_vec_th
  E_lamsq_ybarsq = 1/tau * E_lam_vec_th + E_lam_vec_for*(1+1/a)
  E_lam_y        = E_lam_vec_sq
  E_lamsq_y      = E_lam_vec_th
  E_lam_y_ybar   = E_lam_vec_th * (1/a + 1)
  E_lamsq_y_ybar = E_lam_vec_for * (1+1/a)
  E_lam_ybarsq = 1/tau * E_lam_vec_sq + E_lam_vec_th * (1+1/a)
  E_lamth_ybar = E_lam_vec_for
  E_lamth_y    = E_lam_vec_for
  
  return (c(E_ysq = E_ysq, E_lam_ysq = E_lam_ysq, E_lamsq_ysq = E_lamsq_ysq, E_lamsq_ybar = E_lamsq_ybar, 
            E_lamsq_ybarsq = E_lamsq_ybarsq, E_lam_y = E_lam_y, E_lamsq_y = E_lamsq_y, E_lam_y_ybar = E_lam_y_ybar,
            E_lamsq_y_ybar = E_lamsq_y_ybar, E_lam_ybarsq= E_lam_ybarsq, E_lamth_ybar = E_lamth_ybar, E_lamth_y = E_lamth_y))
}


# 각 모델에 대한 mod2의 알파0,알파1과 mod3의 알파0 계산
P <- prop_coeff(lam_vec, a, tau)
m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lambs^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))
m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))

# m3에 대한 alpha0 
m3_alpha0_n = P["E_lamsq_y_ybar"] - P["E_lamth_y"] - P["E_lamsq_ybarsq"] + P["E_lamth_ybar"]
m3_alpha0_d = P["E_lam_y_ybar"] - P["E_lamsq_y_ybar"] - P["E_lamsq_ybar"] + P["E_lamsq_ybarsq"] - P["E_lamsq_y"]  + P["E_lamth_y"] + mean(lam_vec^3) - P["E_lamth_ybar"]

m3_alpha0 = - m3_alpha0_n/m3_alpha0_d

Yhat = data.frame(matrix(nrow = N_sim, ncol = n_model))

# model 1,2,3 y_hat값 생성
for(n in 1: N_sim){
  Yhat[n,1] = lam_vec[n] * (a+sum(Y[n,1:tau]))/(a+tau*lam_vec[n])
  Yhat[n,2] = lam_vec[n] * (m2_alpha0 + m2_alpha1 * sum(Y[n,1:tau])/tau)
  Yhat[n,3] = lam_vec[n]/(m3_alpha0 + (1- m3_alpha0) * lam_vec[n]) * ( m3_alpha0 + (1- m3_alpha0) * sum(Y[n,1:tau])/tau)
}


lam_vec[1] * (a+sum(Y[1,1:tau]))/(a+tau*lam_vec[1])
lam_vec[1]/(m3_alpha0 + (1- m3_alpha0) * lam_vec[1]) * ( m3_alpha0 + (1- m3_alpha0) * sum(Y[1,1:tau])/tau)


# MSE 계산
mean((Yhat$X1-Y[,6])^2) # 베이지안
mean((Yhat$X2-Y[,6])^2) # 기존 불만
mean((Yhat$X3-Y[,6])^2) # 개선안
#GLM MSE -> 가장 안좋은 결과 산출
(sum((Y[1: (floor(N_sim/3)),6]-lambs[1])^2)+sum((Y[floor(N_sim/3+1): (2*floor(N_sim/3)),6]-lambs[2])^2)+sum((Y[(2*floor(N_sim/3)+1): N_sim,6]-lambs[3])^2))/N_sim


# 각 모델별 unbiasedness 계산
total <- Yhat
total$lam_vec <- lam_vec

# 모델별 lam[i]와 비교
yhat_mean_lam1 = apply(total[total$lam_vec == lambs[1],],2,mean)
yhat_mean_lam2 = apply(total[total$lam_vec == lambs[2],],2,mean)
yhat_mean_lam3 = apply(total[total$lam_vec == lambs[3],],2,mean)

# 모델별 ubiasedness
total_unbiasedness = data.frame(matrix(nrow = 1, ncol = n_model))
for(i in 1: n_model){
  total_unbiasedness[1,i] = (yhat_mean_lam1[i]- lambs[1])^2/lambs[1]  +
                            (yhat_mean_lam2[i]- lambs[2])^2/lambs[2]  +
                            (yhat_mean_lam3[i]- lambs[3])^2/lambs[3]
}

total_unbiasedness

