
#####################################################
# Load data
#####################################################


set.seed(1234)
path <- "C:/Users/none/Desktop/대학원/논문" # Download path 

setwd(path)

path_ftn <- paste0(path, "/Unbiased insurance preminum")

ifelse(dir.exists(paste0(path_ftn,"\\m3_alpha0_plot")),F, dir.create(paste0(path_ftn,"\\m3_alpha0_plot")) ) 

source(paste0(path_ftn,"/lamb_kappa.R"))

data.train <- read.csv('./data/data_train.csv')
data.valid <- read.csv('./data/data_valid.csv')

dim(data.train)
dim(data.valid)

#####################################################
# Import package
#####################################################
library(tidyverse)
library(nimble)
library(MCMCvis)
library(plotrix)
#####################################################
# data preprocessing
#####################################################
# (1) uniq_id for PolicyNum

id_uniq = unique(data.train$PolicyNum)
length(id_uniq) # 497

# (2) idx per id per PolicyNum's Year since each ID has different year of frequency (i.e not the same observations, 5,5,4,...)
id_idx=c(NULL)

for(i in 1:length(id_uniq)){
  ind = id_uniq[i] == data.train$PolicyNum 
  id_idx=c(id_idx, rep(i, sum(ind)))
}

head(id_idx, 20)    
tail(id_idx, 20) 
length(id_idx)


data.train <- combine_kappa_level(data.train)
data.valid <- combine_kappa_level(data.valid)

#####################################################
# estimate a.hat, lambda.hat, R, beta.hat using nimble
#####################################################

round(apply(data.train, 2, mean)*100,2) 

# Make design matrix 
X = model.matrix(n~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+C2+C3, data = data.train)


model <- nimbleCode({
  
  # likelihood
  for(j in 1:JJ){
    lambda[j] <-exp(inprod(X[j,1:KK], beta.hat[1:KK]))
    y[j] ~ dpois( lambda[j] * R[id_idx[j]] )  
  }
  
  # prior of R
  for(n in 1:N.id_idx){
    R[n] ~ dgamma(mean=1, sd=1/sqrt(a.hat))
  }
  
  # prior of Betas 
  for(j in 1:KK){
    beta.hat[j] ~ dnorm(mean=0, sd=5)
  } 
  
  # prior of alpha
  a.hat ~ dunif(min=0, max=10)
})


my.data <- list(X = X, y = data.train$n)
my.constants <- list(JJ=dim(X)[1], KK=dim(X)[2],id_idx=id_idx, N.id_idx=max(id_idx))
parameters.to.save <- c("beta.hat", "a.hat", "lambda", "R")


my_model <- nimbleModel(code = model, name = "my_model", constants = my.constants, data = my.data)
Cmy_model  <- compileNimble(my_model)
my_modelConf <- configureMCMC(my_model, print = TRUE)
my_modelConf$addMonitors('beta.hat', "a.hat", "lambda", "R")
MCMC <- buildMCMC(my_modelConf)
CMCMC <- compileNimble(MCMC, project = my_model)


n.iter <- 12000
n.burnin <- 2000
n.chains <- 2

samples <- runMCMC(CMCMC, niter = n.iter,
                   nburnin = n.burnin,
                   nchains = n.chains)

MCMCtrace(object = samples,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "beta.hat")
dev.off()
MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "a.hat")



# save each estimation ( ex, R.hat, a.hat, beta.hat, k.hat, lambda.hat)
# total num of estimation = 497 + 1 + 8 + 1 + 2057

str(samples) # property of MCMC parameters

# round 2 point for each parameters of 10000 samples 
result <- MCMCsummary(object = samples, round = 2)

# check betas estimation
result[499:506,]

head(result); dim(result)

# save estimated params
R <- result[1:497,1] ;R #497
a.hat <- result[498,1] ;a.hat #1
beta.hat <-  result[499:506,1] ;beta.hat #8
lambda.hat <- result[507:2581,1] # 2075




# risk ratio per PolcyNum
summary(lambda.hat)
table(lambda.hat) # num = 17 since the one combination is missing

#####################################################
# estimation of beta.hat[1]-beta.hat[8]
# combine parameter all estimation( beta.hat,  a.hat)
#####################################################
beta.est.mean <- result[499:506,1];
beta.est.sd <- result[499:506,2];
mu <- rep(0,8)
pval <- pnorm(-abs((beta.est.mean - mu)/beta.est.sd))
beta.pval <-round(pval*2, 4)

beta.df <- data.frame(cbind(beta.est.mean= beta.est.mean, beta.est.sd = beta.est.sd,beta.pval = beta.pval))

a.hat.result = k.hat.result <- rep(NA, dim(beta.df)[2])
a.hat.result[1] <- a.hat;
param_result <- rbind(beta.df, a.hat)
rownames(param_result) <- c('beta0','TypeCity','TypeCounty','TypeSchool','TypeTown','TypeVillage',"c2","C3",'a.hat')

param_result

#####################################################
# Fill missing data for Y's using Y's distn assumption(ex rpois(mean= lambda * R, phi) with sampled R.hat , lambda.hat
#####################################################

# Y 
tau = 5 ; Y = data.frame(matrix(nrow = length(id_uniq), ncol = tau))


for(i in 1:length(id_uniq)){
  Y[i,tau+1] = id_uniq[i]
  ys = data.train[data.train$PolicyNum == id_uniq[i], ]$n
  for(j in 1: length(ys)){
    Y[i,j] = ys[j]
  }
}

colnames(Y) <- c("y_1", "y_2", "y_3",  "y_4", "y_5", "ID")

# Sample R.hat and select randomly R.hat[1] for each PolicyNum


R.hat <- rep(0, 497)
for(i in 1:497){
  R.hat[i] <- sample(samples$chain1[,paste0("R[", i, "]")],1)
}

lambda_id = data.frame(unique(cbind(data.train[,c("PolicyNum","lamb_kappa")], lambda.hat))); colnames(lambda_id) <- c('ID',"lambda_kappa_level", 'lambda.hat_per_id')

############################################
# chk1) Ratio of Risk group for data.train
###########################################
table(lambda_id[,'lambda.hat_per_id']);
round(table(lambda_id[,'lambda.hat_per_id'])/497,3);

for (i in 1:nrow(Y)){
  for (j in 1:tau){
    if (is.na(Y[i,j])==TRUE)
      Y[i,j] = rpois(1,lambda_id[,'lambda.hat_per_id'][i]*R.hat[i])}
}
head(Y)

# total Y : columns = ID, y_1,y_2,y_3,y_4, y_5, lambda.hat_per_ID
total <- merge(Y, lambda_id, by = 'ID')
total <- merge(total, data.valid[,c('PolicyNum','n')], by.x = "ID", by.y = 'PolicyNum') 

head(total, 20)


############################################
# chk2) E[Y|lamb = lambda_kappa] , lamb.hat, lamb_level
###########################################
# import function to combine lamb_kappa with data.valid
source(paste0(path_ftn,"/lamb_kappa.R"))

head(total)

lambda_kappa_level <- unique(total[,"lambda_kappa_level"])
i = 1

lamb.hat.csv <- data.frame(lambda_kappa_level = double(),
                           lambda_hat = double(),
                           E_y = double()
)

for( i in 1:length(lambda_kappa_level)){
  total_id.per.kappa <- total[total$lambda_kappa_level == lambda_kappa_level[i],]
  size_y <- dim(total_id.per.kappa[, c('y_1','y_2','y_3','y_4','y_5')])
  sum_y  <- sum(total_id.per.kappa[, c('y_1','y_2','y_3','y_4','y_5')])
  E_y    <- sum_y/(size_y[1] * size_y[2])
  
  lamb.hat.csv[i,1] <- total_id.per.kappa$lambda_kappa_level[1]
  lamb.hat.csv[i,2] <- total_id.per.kappa$lambda.hat_per_id[1]
  lamb.hat.csv[i,3] <- E_y
}


# number of each lambda_kappa 
table(total[,'lambda_kappa_level'])
lamb.hat.csv <- lamb.hat.csv[order(lamb.hat.csv$lambda_kappa_level),]


############################################
# chk2) Ratio of Risk group for data.valid
###########################################
idx = unique(data.train$PolicyNum) %in% data.valid$PolicyNum

table(lambda_id[idx, ]$lambda);
round(table(lambda_id[idx, ]$lambda.hat_per_id)/dim(data.valid)[1], 3);


#####################################################
# Import function to measure Cred for Bulhmann, Comm, Generalized Comm with Pois-Gamma Random Effect Model
#####################################################

# fill the Y given estimated R.hat

Y <- total[,c('y_1','y_2','y_3','y_4','y_5','n')]

a = a.hat; # param of Gamma
N_sim = dim(Y)[1]; # 379
tau = 5; 
n_model = 3; 
lam_vec <- lambda_id[idx,'lambda.hat_per_id']; #379



## Builhmann, commercial, Generalized commercial의 예측값, MSE, DIX 계산

# import function to calculate each model
source(paste0(path_ftn,"/Pois_Gamma_RE.R"))

P <- prop_coeff(lam_vec, a, tau)
m2_alpha1 = (P["E_lam_y"] * P["E_lamsq_ybar"] - P["E_lam_y_ybar"] * mean(lam_vec^2))       /(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))
m2_alpha0 = (P["E_lamsq_ybar"] * P["E_lam_y_ybar"] - P["E_lam_y"] * P["E_lamsq_ybarsq"])/(P["E_lamsq_ybar"]^2-P["E_lamsq_ybarsq"]*mean(lam_vec^2))

result_alpha0 <- optim(par = c(0), 
                fn = optim_alpha , 
                lam_vec = lam_vec, Y = Y,
                control = list(fnscale = 1), method="BFGS")


(optim_alpha(lam_vec, Y, par=c(result_alpha0$par)))

m3_alpha0 <- result_alpha0$par

##########################################################
# plotting estimation of alpha0 for generalized comm
##########################################################
m3_alpha0_seq <- seq(0,1,length = 30)
MSE_alpha0 <- rep(0, length(m3_alpha0_seq))

for(i in 1: length(m3_alpha0_seq)){

  Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0 = m3_alpha0_seq[i])
  # calculate HMSE of Gprem to check which is the best alpha0
  MSE_alpha0[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'MSE', Prem = Yhat[,3])
}

# best HMSE for alpha0
Yhat <- Yhat_ftn( N_sim, n_model, tau, lam_vec, Y, m2_alpha0, m2_alpha1, m3_alpha0 = m3_alpha0)
best_MSE <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'MSE', Prem = Yhat[,3])

# save png for each alpha0

new_folder_name <- paste0(path_ftn,"\\m3_alpha0_plot\\", format(Sys.Date(),"%m%d"),"_Real")

if (!file.exists(new_folder_name)){
  dir.create(new_folder_name)
}

setwd(new_folder_name)
png(filename='real.nb.png',width=2200,height=2200)

#dev.new()

x_axis_tick=seq(0,1,length=5)
y_axis_tick=round(seq(min(MSE_alpha0),max(MSE_alpha0),length=5),2)

png(filename='real.pois.png',width=2200,height=2000)
par(oma = c(1,1,1,1))
par(mar=c(25,25,25,5))
oldp <- par(mgp = c(3,6.5,0.5))


plot(m3_alpha0_seq,MSE_alpha0,type = 'l', xlab =''
     , ylab = "",xlim = c(0,1),
     ylim = c(min(MSE_alpha0),max(MSE_alpha0))
     ,axes=FALSE, xaxt = 'n',yaxt = 'n', lwd = 5)

axis(side=1,at=x_axis_tick, cex.axis = 10, lwd = 4)
axis(side=2,at=y_axis_tick, cex.axis = 10, lwd = 4, mgp = c(0,5,2))

mtext(expression(alpha[0]), side = 1, line= 22, cex = 12)
mtext("HMSE", side = 2, line= 15, cex = 12)

points(m3_alpha0, best_MSE,pch = 17, cex = 8)


dev.off()


MAE = HMSE = MSE <- rep(0, n_model)
for(i in (1:n_model)){
  MSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R, type = 'MSE', Prem = Yhat[,i])
  HMSE[i] <- cal_pred(lam_vec, Y = Y[,tau+1], Theta = R.hat[idx], type='HMSE', Prem = Yhat[,i])
}

# Unbiasedness for CPrem DIX
m2_DIX <- cal_UB(lam_vec, a, alpha0 = m2_alpha0,alpha1 = m2_alpha1 )

# to data.frame each criterion value(MSE, HMSE) for each model(Spre, Cprem, GCPrem)
criterion.result <- data.frame(rbind(MSE,HMSE)); 
colnames(criterion.result) <- c("Sprem","Cprem","GCprem")
criterion.result <- rbind(criterion.result,c(0,m2_DIX[1],0))
rownames(criterion.result)[4] <- 'DIX'


# save all results

ifelse(dir.exists(paste0(path_ftn,"\\result")),F, dir.create(paste0(path_ftn,"\\result")) ) 

new_result_pois <- paste0(path_ftn,"\\result\\", format(Sys.Date(),"%m%d"),"_Real")

if (!file.exists(new_result_pois)){
  dir.create(new_result_pois)
}

setwd(new_result_pois)

filename <- paste0("/",format(Sys.Date(),"%m%d"), "_real_Pois_gamm_result.csv")

write.csv(criterion.result, file = paste0(new_result_pois,"/real_Pois_gamm_result.csv"))


filename <- paste0("/",format(Sys.Date(),"%m%d"), "_pois_param_est_result.csv")

write.csv(param_result, file = paste0(new_result_pois,"/pois_param_est_result.csv"))


filename <- paste0("/",format(Sys.Date(),"%m%d"), "_pois_lamb.csv")
write.csv(lamb.hat.csv, file = paste0(new_result_pois,"/pois_lamb.csv"))
