# -------------------------------------------------------------------------------
# Simulation for NB-Gamma random effect Model
# -------------------------------------------------------------------------------

## Calculate premium, DIX, HMSE of Builhmann, commercial, Generalized commercial

# import function to calculate each model

#source("C:/Users/none/Desktop/대학원/논문/Unbiased insurance preminum/NB_Gamma_RE.R") 
source("C:/Users/parkminji/Downloads/민지/Unbiased insurance preminum/NB_Gamma_RE.R")
total_lam <- matrix(c(0.1,0.4,0.6,0.1,0.5,0.5,1.0,0.2,0.9,0.6,1.4,1.2),nrow = 4, ncol =3)

a = 0.8; N_sim = 10000; lambs = c(0.1, 0.2, 1.2); tau = 5; n_model = 3; k = 2


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
m3_alpha0_seq <- seq(0,1,length = 20)
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


plot(m3_alpha0_seq, HMSE_alpha0)
title(main= paste0("lambs =", lambs))
points(m3_alpha0, best_HMSE,pch = 19, col ='blue')



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


HMSE;
MSE;
m2_DIX


