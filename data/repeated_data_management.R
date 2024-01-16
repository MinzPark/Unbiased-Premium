######################################################################
# data prerpocessing for data.train
######################################################################

# Define function for preprocessing data.train and data.valid

############################################################
# divide into 3 category class for collision Coverage of data.train and data.valid
############################################################

Cov_idx <- function(data.train){
  mycov<-sort(data.train$col.Cov)
  length(mycov);
  
  data.train$col.Cov.Idx<- as.factor( (data.train$col.Cov>mycov[dim(data.train)[1] /3  ]) +
                                        (data.train$col.Cov>mycov[dim(data.train)[1] /3*2])  )
  unique(data.train$col.Cov.Idx)
  table(data.train$col.Cov.Idx);
  
  # set equally the same Coverage.Idx each PolicyNum 
  
  id_uniq = unique(data.train$PolicyNum)
  length(id_uniq) # 497
  
  not_equal_cov<-col.cov.idx <- c()
  for(i in 1:length(id_uniq)){
    tmp <- data.train[id_uniq[i] == data.train$PolicyNum, ] # select data.train per PolicyNum
    
    # check whose not the same Coverage.Idx
    table_cov_per_id <- table(tmp$col.Cov.Idx)
    if (sum(table(tmp$col.Cov.Idx)%%dim(tmp)[1]) != 0){ not_equal_cov <-c(not_equal_cov, id_uniq[i]) }
    
    # if each PolicyNum does not have the same Coverage.Idx, revert to the maximal frequency of Coverage.Idx
    col.cov.idx <- c(col.cov.idx, rep(as.numeric(names(which.max(table(tmp$col.Cov.Idx)))) ,dim(tmp)[1])) 
  }
  
  
  data.train$col.Cov.Idx <- col.cov.idx
  
  # one-hot encoding for categorical variable (Coverage.Idx)
  data.train$C2  <- ifelse(data.train$col.Cov.Idx==1, 1, 0)
  data.train$C3  <- ifelse(data.train$col.Cov.Idx==2, 1, 0)
  
  return(data.train)
}


refine_ex_var <- function(data.train){
  
  id_uniq = unique(data.train$PolicyNum)
  param <- c("TypeCity","TypeCounty", "TypeSchool","TypeTown","TypeVillage")
  not_equal_type <- c()
  
  for(i in (1:length(id_uniq))){
    id_idx <- id_uniq[i] == data.train$PolicyNum
    df_id <- data.train[id_idx, ]
    ref_ex_df <- df_id[,param]
    chk <- apply(ref_ex_df,2,sum)
    
    
    if( sum(chk %% dim(df_id)[1]) >0 ){
      not_equal_type <- c(not_equal_type, id_uniq[i])
      
      mode_col <- names(which.max(df_id))
      data.train[id_idx, mode_col] <- 1; data.train[id_idx, param[!param %in% mode_col]] <- 0
    }
  }
  
  return(data.train)
}



############################################################
# data load
############################################################
# setwd("C:/Users/none/Desktop/대학원/논문/data")
setwd(getwd())
load("./data/data.RData")
head(data)

# (1) data load and define data.train

data$col.freq<-data$FreqCN+data$FreqCO         #Collision old and new
data$col.Cov<-data$CoverageCN+data$CoverageCO  #Coverage old and new

index<-(data$CoverageCN>0 | data$CoverageCO>0)
data.train<-data[index,]



# (2) 1) refine the explantory variable
#     2) divide into 3 category class for collision Coverage of data.train


data.train <- Cov_idx(data.train)
data.train <- refine_ex_var(data.train)

# (3) Define frequency and severity of collision of old and new
data.train$n <-data.train$col.freq
data.train$s <- data.train$yAvgCN * data.train$FreqCN   +  data.train$yAvgCO * data.train$FreqCO

# # Define frequency and severity of collision of old and new
# index2<-(data.train$col.freq>0)
# data.train.sev<-data.train[index2,]
# data.train.sev$m<- (data.train.sev$yAvgCN * data.train.sev$FreqCN   +  data.train.sev$yAvgCO * data.train.sev$FreqCO)/(data.train.sev$FreqCN+data.train.sev$FreqCO)


# (5) Delete except explanatory parameter:  "PolicyNum", "Year", "TypeCity","TypeCounty", "TypeSchool","TypeTown","TypeVillage","col.Cov.Idx" and target: "n"

data.train <- data.train[,c("PolicyNum", "Year", "TypeMisc","TypeCity","TypeCounty", "TypeSchool","TypeTown","TypeVillage","col.Cov.Idx", "C2","C3","n")]

round(apply(data.train,2,mean)*100,2)

######################################################################
# data prerpocessing for data.test
######################################################################

load("./data/dataout.RData")
head(dataout)


# Out of test data, use index which is from the training data & coverage>0
out.idx1 <- dataout$PolicyNum %in% data.train$PolicyNum # filtering test data to fit data.train indicator
out.idx2 <- dataout$CoverageCO>0 | dataout$CoverageCN>0 # focus on collision coverage
out.idx <- out.idx1 & out.idx2
length( unique(data.train$PolicyNum) )
sum(out.idx)
data.valid<-dataout[out.idx,] #test data


# Now, define explanatory variable in test data. 
# (I am not sure the explanatory variable will be ever used or not.)
out.idx.sev<- data.valid$PolicyNum %in% data.train$PolicyNum
data.valid$col.freq<-data.valid$FreqCN+data.valid$FreqCO
data.valid$col.Cov<-data.valid$CoverageCN+data.valid$CoverageCO


# divide into 3 category class for collision Coverage of data.valid
data.valid <- Cov_idx(data.valid)


#data.valid$col.Cov.Idx<- as.factor((data.valid$col.Cov>0.635)  +1)
data.valid$s <- data.valid$yAvgCN * data.valid$FreqCN   +  data.valid$yAvgCO * data.valid$FreqCO
data.valid$n <- data.valid$FreqCN   +   data.valid$FreqCO

data.valid <- data.valid[,c("PolicyNum", "Year","TypeMisc", "TypeCity","TypeCounty", "TypeSchool","TypeTown","TypeVillage","col.Cov.Idx", "C2","C3","n")]

# ###For the Research with Rosy
# typeof(data.train$s)
# head(data.train, 1)
# idx<-data.train$Year==2010
# 
# # ??  why to do? 
# sum(data.train$s[idx] )/sum(data.train[idx,]$n )
# sum(data.valid$s )/sum(data.valid$n )



## save data file
#write.xlsx( x = data.train, file = "data_train.xlsx", rowNames = TRUE) # data.train
#write.xlsx( x = data.valid, file = "data_valid.xlsx", rowNames = TRUE) # data.valid

