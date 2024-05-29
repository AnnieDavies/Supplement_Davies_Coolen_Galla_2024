#function to fit the models & do cross validation for JM and LM
run_sim <- function(dat_split){
  dat <- dat_split[[1]]
  splits <- dat_split[[2]]
  
  dat.id <- dat[!duplicated(dat$id), ]
  n <- nrow(dat.id)
  
  ## STEP 2: fit models to full data set
  JM1.models <- fit_JM1(dat, dat.id)
  JM1 <- JM1.models[[1]]
  JM2 <- fit_JM2(JM1.models)
  
  paramsJM1 <- JM_params(JM1, model = 1) #alpha, gamma
  paramsJM2 <- JM_params(JM2, model = 2) #alpha, gamma
  
  ## 2a store results in dataframe
  df.alldata <- data.frame(matrix(nrow = 2, ncol = 5))
  colnames(df.alldata) <- c("Method", "alpha", "gamma", "tau")
  
  df.alldata[1, ] <- c("JM1", paramsJM1[[1]], paramsJM1[[2]], NA)
  df.alldata[2, ] <- c("JM2", paramsJM2[[1]], paramsJM2[[2]], NA)
  
  ## 2b create data frames to store group level results
  df.LM <- data.frame(matrix(nrow = 10, ncol = 6))
  colnames(df.LM) <- c("Group", "PE1", "PE2", "PE3", "PE4", "PE5")
  
  df.JM1 <- data.frame(matrix(nrow = 10, ncol = 8))
  colnames(df.JM1) <- c("Group", "alpha", "gamma", "PE1", "PE2", "PE3", "PE4", "PE5")
  df.JM2 <- df.JM1
  
  ## STEP 3: Cross validation
  for(g in 1:10){
    #create test and training data------------------
    test_id <- splits[[g]]
    test_data <- dat[dat$id %in% test_id, ]
    train_data <- dat[!dat$id %in% test_id, ]
    test_data.id <- test_data[!duplicated(test_data$id), ]
    train_data.id <- train_data[!duplicated(train_data$id), ]
    
    ## STEP 3a: Fit models
    train_JM1.models <- fit_JM1(train_data, train_data.id)
    train_JM1 <- train_JM1.models[[1]]
    train_JM2 <- fit_JM2(train_JM1.models)
    
    train_paramsJM1 <- JM_params(train_JM1, model = 1) #alpha, gamma
    train_paramsJM2 <- JM_params(train_JM2, model = 2) #alpha, gamma
    
    ## Step 3b: prediction error
    #order of vectors: vec.pe.LM, vec.pe.JM1, vec.pe.JM2
    PE.vecs <- run_PE(test_data, train_data, train_JM1, train_JM2)
    
    vec.pe.LM <- PE.vecs[[1]]
    vec.pe.JM1 <- PE.vecs[[2]]
    vec.pe.JM2 <- PE.vecs[[3]]
    
    ##assign values to dfs
    df.LM <- add_df_LM(df.LM, vec.pe.LM, g)
    
    df.JM1 <- add_df_JM(df.JM1, vec.pe.JM1, g, train_paramsJM1)
    df.JM2 <- add_df_JM(df.JM2, vec.pe.JM2, g, train_paramsJM2)
    
  }

  #define list of results
  JMLM_results <- list(df.alldata, df.LM, df.JM1, df.JM2)
  #return list
  JMLM_results
}


add_df_LM <- function(df, vec.pe, g){
  df$Group[g] <- g
  
  df$PE1[g] <- vec.pe[1]
  df$PE2[g] <- vec.pe[2]
  df$PE3[g] <- vec.pe[3]
  df$PE4[g] <- vec.pe[4]
  df$PE5[g] <- vec.pe[5]
  
  df
}

add_df_JM <- function(df, vec.pe, g, params){
  df$Group[g] <- g
  
  df$alpha[g] <- params[[1]]
  df$gamma[g] <- params[[2]]
  
  df$PE1[g] <- vec.pe[1]
  df$PE2[g] <- vec.pe[2]
  df$PE3[g] <- vec.pe[3]
  df$PE4[g] <- vec.pe[4]
  df$PE5[g] <- vec.pe[5]
  
  df
}
