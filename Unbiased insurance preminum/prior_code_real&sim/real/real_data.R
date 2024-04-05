
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
    Yhat[n,1] = lam_vec[n] * (a+sum(Y[n,1:tau]))/(a+tau*lam_vec[n])
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
  total$lam_vec <- lam_vec # tmp <- cbind(Yhat, lam_vec)이건 왜 안됨..?
  for( i in 1: length(lambs)){
    E_y_given_lam[i,] = apply(total[total$lam_vec == lambs[i],],2,mean) # 각 모델의 lamb별 평균
  }
  
  total_unbiasedness_ch = apply((E_y_given_lam[,1:n_model]-lambs)^2/lambs,2,sum)
  return(total_unbiasedness_ch)
}

# -------------------------------------------------------------------------------
# 실제 데이터 
# -------------------------------------------------------------------------------

# 데이터 불러오기
setwd("C:/Users/none/Desktop/대학원/논문/data")
getwd()
source("repeated_data_management.R")
head(data.train)
sum(is.na(data.train$n)) # 미싱데이터 없음

# beta_hat, a_hat, lam_hat, p_hat 추정
uniq_pol = unique(data.train$PolicyNum)

person = c(NULL)
for(i in 1: length(uniq_pol)){
  ind = uniq_pol[i] == data.train$PolicyNum
  person = c(person, rep(i, sum(ind)))
}

library(nimble)
library(MCMCvis)


# -------------------------------------------------------------------------------
# X maxtrix 만들기 : City, County, Village, School, Town, col.Cov.IDx
# -------------------------------------------------------------------------------
inc <- function(x)
{
  eval.parent(substitute(x <- x + 1))
}
col2 = col3 = rep(0, times = nrow(data.train))
j = 1
for(i in data.train$col.Cov.Idx){
  if ( as.numeric(i) > 1 ){
      col2[j] = 1
    if(as.numeric(i) > 2){
      col3[j] =1} 
  }
  j = inc(j)
}

c(data.train$col.Cov.Idx,col2, col3)

data.train$col2 <-col2
data.train$col3 <-col3

x_df <- data.frame(cbind(rep(1,nrow(data.train)), data.train[,c("TypeCity","TypeCounty","TypeVillage","TypeSchool","TypeTown","col2","col3")])) # [상수항, params]
X = data.matrix(x_df)
head(X)
X_ori = model.matrix(n~ TypeCity+TypeCounty+ TypeVillage+TypeSchool+TypeTown+factor(col.Cov.Idx), data = data.train) #(col.Cov.Idx가 [(0,0),(0,1),(1,1)]조합)
head(X_ori)

model <- nimbleCode({
  for(j in 1:JJ){
    lambda[j] <- exp(inprod(X[j,1:KK], beta.hat[1:KK]))
    y[j] ~ dnegbin(size = k, prob = k/(k+lambda[j]*R[Person[j]]))
  }
  #Prior of R
  for(n in 1: N.person){
    R[n] ~ dgamma(mean = 1, sd = 1/a.hat)
  }
  
  # prior of Betas
  for(j in 1:KK){
    beta.hat[j] ~ dnorm(mean =0, sd = 10)
  }
  
  # prior of a
  a.hat ~ dunif(min = 0 , max = 10)
})

my.data <- list( X = X, y = data.train$n)
my.constants <- list(JJ = dim(X)[1], KK = dim(X)[2], Person = person, N.person = max(person))
parameters.to.save <- c("beta.hat", "a.hat", "lambda")

pois_model <- nimbleModel(code = model, name = "pois_model", constants = my.constants, data = my.data)
showCompilerOutput = TRUE
Cpois_model <- compileNimble(pois_model)

printErrors()

pois_modelConf <- configureMCMC(pois_model, print = TRUE)

pois_modelConf$addMonitors('beta.hat', 'a.hat','lambda')

MCMC <- buildMCMC(pois_modelConf)
CMCMC <- compileNimble(MCMC, project = pois_model)

n.iter <- 10000
n.burnin <- 2000
n.chains <- 2

set.seed(123)
samples <- runMCMC(CMCMC, niter = n.iter, nburnin = n.burnin, nchains = n.chains)

MCMCtrace(object = samples, pdf = FALSE, ind = TRUE, params = "beta.hat")
MCMCtrace(object = samples, pdf = FALSE, ind = TRUE, params = "a.hat")

str(samples[[1]])

result <- MCMCsummary(object = samples, round = 2)
dim(result)

# save estimated params

a.hat <- result[1,1]
beta.hat <- result[2:9,1]; # beta.hat
lambda.hat <- result[10:1476,1]

# ---------------------------------------------------------------------------
N_sim = length(uniq_pol); 
a = a.hat;lambs =unique(lambda.hat); tau = length(unique(data.train$Year)); n_model = 3
k =2 # NB의 param
lam_vec = lambda.hat



# Y 구성

Y = data.frame(matrix(nrow = N_sim, ncol = tau))
ID = unique(data.train$PolicyNum); 
for(i in 1:length(ID)){
  Y[i,tau+1] = ID[i]
  ys = data.train[data.train$PolicyNum == ID[i], ]$n
  for(j in 1: length(ys)){
    Y[i,j] = ys[j]
  }
}
#Y$n <- data.valid$n # target은 이름이 n인 열로 구성

Y

# 각 모델에 대한 mod2의 알파0,알파1과 mod3의 알파0 계산
P <- prop_coeff(lam_vec, a, tau, k)


m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lambs^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))
m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lambs^2))

# m3에 대한 alpha0 
m3_alpha0_n = P["E_lamsq_y_ybar"] - P["E_lamth_y"] - P["E_lamsq_ybarsq"] + P["E_lamth_ybar"]
m3_alpha0_d = P["E_lam_y_ybar"] - P["E_lamsq_y_ybar"] - P["E_lamsq_ybar"] + P["E_lamsq_ybarsq"] - P["E_lamsq_y"]  + P["E_lamth_y"] + mean(lam_vec^3) - P["E_lamth_ybar"]

m3_alpha0 = - m3_alpha0_n/m3_alpha0_d


# 각 모델별 예측값 추출 
Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0)

# MSE 계산
mse_model <- cal_mse(n_model, y_true = Y[,tau+1])

# 각 모델별 unbiasedness 확인
total_unbiasedness_ch <- cal_UB(lambs,n_model,Yhat, lam_vec)

total_unbiasedness_ch

