export_JMdata <- function(data.list, scenario, JM, data_path){
  N <- length(data.list)
  
  alpha <- c()
  gamma <- c()
  PE1 <- c()
  PE2 <- c()
  PE3 <- c()
  PE4 <- c()
  PE5 <- c()
  
  df.results <- data.frame(
    alpha = numeric(N),
    gamma = numeric(N),
    avPE1 = numeric(N),
    avPE2 = numeric(N),
    avPE3 = numeric(N),
    avPE4 = numeric(N),
    avPE5 = numeric(N)
  )
  
  
  for(i in 1:length(data.list)){
    #parameter estimates
    df1 <- data.list[[i]][[1]]
    df.results$alpha[i] <- as.numeric(df1$alpha[JM]) #JM1 is row 1, JM2 is row2
    df.results$gamma[i] <- as.numeric(df1$gamma[JM]) #JM1 is row 1, JM2 is row2
    
    #Prediction error
    df2 <- data.list[[i]][[JM+2]] #JM1 is df[3], JM2 is df[4]
    #average over group
    df.results$avPE1[i] <- mean(as.numeric(df2$PE1), na.rm=TRUE)
    df.results$avPE2[i] <- mean(as.numeric(df2$PE2), na.rm=TRUE)
    df.results$avPE3[i] <- mean(as.numeric(df2$PE3), na.rm=TRUE)
    df.results$avPE4[i] <- mean(as.numeric(df2$PE4), na.rm=TRUE)
    df.results$avPE5[i] <- mean(as.numeric(df2$PE5), na.rm=TRUE)
    
    #per group
    alpha <- append(alpha, as.numeric(df2$alpha))
    gamma <- append(gamma, as.numeric(df2$gamma))
    PE1 <- append(PE1, as.numeric(df2$PE1))
    PE2 <- append(PE2, as.numeric(df2$PE2))
    PE3 <- append(PE3, as.numeric(df2$PE3))
    PE4 <- append(PE4, as.numeric(df2$PE4))
    PE5 <- append(PE5, as.numeric(df2$PE5))
  }
  df.group.results <- data.frame(
    alpha = alpha,
    gamma = gamma,
    PE1 = PE1,
    PE2 = PE2,
    PE3 = PE3,
    PE4 = PE4,
    PE5 = PE5
  )
  write.csv(df.results, file = paste(data_path, "S", scenario, "\\JM", JM, "\\results_JM", JM, ".csv", sep=""))
  write.csv(df.group.results, file = paste(data_path, "S", scenario, "\\JM", JM, "\\group_results_JM", JM, ".csv", sep=""))
  
}

export_LMdata <- function(data.list, scenario, data_path){
  N <- length(data.list)
  
  PE1 <- c()
  PE2 <- c()
  PE3 <- c()
  PE4 <- c()
  PE5 <- c()
  
  df.results <- data.frame(
    avPE1 = numeric(N),
    avPE2 = numeric(N),
    avPE3 = numeric(N),
    avPE4 = numeric(N),
    avPE5 = numeric(N)
  )
  
  
  for(i in 1:length(data.list)){
    
    #Prediction error
    df2 <- data.list[[i]][[2]] #LM is df[2]
    #average over group
    df.results$avPE1[i] <- mean(as.numeric(df2$PE1), na.rm=TRUE)
    df.results$avPE2[i] <- mean(as.numeric(df2$PE2), na.rm=TRUE)
    df.results$avPE3[i] <- mean(as.numeric(df2$PE3), na.rm=TRUE)
    df.results$avPE4[i] <- mean(as.numeric(df2$PE4), na.rm=TRUE)
    df.results$avPE5[i] <- mean(as.numeric(df2$PE5), na.rm=TRUE)
    
    #per group
    PE1 <- append(PE1, as.numeric(df2$PE1))
    PE2 <- append(PE2, as.numeric(df2$PE2))
    PE3 <- append(PE3, as.numeric(df2$PE3))
    PE4 <- append(PE4, as.numeric(df2$PE4))
    PE5 <- append(PE5, as.numeric(df2$PE5))
  }
  df.group.results <- data.frame(
    PE1 = PE1,
    PE2 = PE2,
    PE3 = PE3,
    PE4 = PE4,
    PE5 = PE5
  )
  write.csv(df.results, file = paste(data_path, "S", scenario, "\\LM\\results_LM.csv", sep=""))
  write.csv(df.group.results, file = paste(data_path, "S", scenario, "\\LM\\group_results_LM.csv", sep=""))
  
}


export_DKdata_python <- function(data.list, scenario, DK, data_path){
  
  if(DK=="A"){
    ind <- 1
  }else if(DK =="B"){
    ind <- 2
  }else{
    print("DK must be A or B")
    ind <- NA
  }
  
  N <- length(data.list)
  
  kappa <- c()
  tau <- c()
  gamma <- c()
  PE1 <- c()
  PE2 <- c()
  PE3 <- c()
  PE4 <- c()
  PE5 <- c()
  
  df.results <- data.frame(
    alpha = numeric(N),
    gamma = numeric(N),
    tau = numeric(N),
    avPE1 = numeric(N),
    avPE2 = numeric(N),
    avPE3 = numeric(N),
    avPE4 = numeric(N),
    avPE5 = numeric(N)
  )
  
  
  for(i in 1:length(data.list)){
    #parameter estimates
    df1 <- reticulate::py_to_r(data.list[[i]][[1]])
    df.results$alpha[i] <- as.numeric(df1$alpha[ind]) #Model A is row 1, Model B is row 2
    df.results$gamma[i] <- as.numeric(df1$gamma[ind]) #Model A is row 1, Model B is row 2
    df.results$tau[i] <- as.numeric(df1$tau[ind]) #Model A is row 1, Model B is row 2
    
    #Prediction error
    df2 <- reticulate::py_to_r(data.list[[i]][[ind+1]]) #A is df[2], B is df[3]
    #average over group
    df.results$avPE1[i] <- mean(as.numeric(df2$PE1), na.rm=TRUE)
    df.results$avPE2[i] <- mean(as.numeric(df2$PE2), na.rm=TRUE)
    df.results$avPE3[i] <- mean(as.numeric(df2$PE3), na.rm=TRUE)
    df.results$avPE4[i] <- mean(as.numeric(df2$PE4), na.rm=TRUE)
    df.results$avPE5[i] <- mean(as.numeric(df2$PE5), na.rm=TRUE)
    
    #per group
    kappa <- append(kappa, as.numeric(df2$kappa))
    tau <- append(tau, as.numeric(df2$tau))
    gamma <- append(gamma, as.numeric(df2$gamma))
    PE1 <- append(PE1, as.numeric(df2$PE1))
    PE2 <- append(PE2, as.numeric(df2$PE2))
    PE3 <- append(PE3, as.numeric(df2$PE3))
    PE4 <- append(PE4, as.numeric(df2$PE4))
    PE5 <- append(PE5, as.numeric(df2$PE5))
  }
  df.group.results <- data.frame(
    kappa = kappa,
    tau = tau,
    gamma = gamma,
    PE1 = PE1,
    PE2 = PE2,
    PE3 = PE3,
    PE4 = PE4,
    PE5 = PE5
  )
  write.csv(df.results, file = paste(data_path, "S", scenario, "\\DK", DK, "\\fromPython_results_DK", DK, ".csv", sep=""))
  write.csv(df.group.results, file = paste(data_path, "S", scenario, "\\DK", DK, "\\fromPython_group_results_DK", DK, ".csv", sep=""))
  
}


convert_py_res <- function(data.list){
  res.list <- list()
  
  for(i in 1:length(data.list)){
    #parameter estimates
    df1 <- reticulate::py_to_r(data.list[[i]][[1]])
    #Model A
    df2 <- reticulate::py_to_r(data.list[[i]][[2]])
    #Model B
    df3 <- reticulate::py_to_r(data.list[[i]][[3]])
    
    df.list <- list(df1, df2, df3)
    res.list[[i]] <- df.list
  }
  res.list
}