run_PE <- function(test_data, train_data, jointFit1, jointFit2){
  print("in run PE")
  vec.pe.LM <-c()
  vec.pe.JM1 <-c()
  vec.pe.JM2 <-c()
  
  win <- 2.0
  for(i in 0:4){
    
    timeBase <- 1.5+(i*2.0) #base time t*
    print(paste("timebase=", timeBase))
    ######################
    # LM model
    ######################
    
    #create LM data set at LM time = timeBase using training data
    train.LM <- dataLM(train_data, timeBase, respVar = "y", timeVar = "time", 
                       evTimeVar = "Time", idVar = "id", summary = "value")
    #Fit a standard Cox model to the landmark data set
    Cox.train.LM <- coxph(Surv(Time, event) ~ group + y, data = train.LM)
    
    
    #####################
    #PE using train_data
    #####################
    #LM
    print("PE LM")
    PE.LM <-PE.AD.coxph(Cox.train.LM, newdata = test_data, Tstart = timeBase, Thoriz = timeBase + win,
                      idVar = "id", timeVar = "time", respVar = "y", evTimeVar = "Time",
                      lossFun = "square", summary = "value")
    
    #JM1:
    print("PE JM1")
    PE.JM1<-PE.AD.JM(jointFit1, newdata = test_data, Tstart = timeBase, Thoriz = timeBase + win,
                      idVar = "id", timeVar = "time", respVar = "y", evTimeVar = "Time",
                      lossFun = "square", summary = "value", simulate = TRUE)

    #JM2:
    print("PE JM2")
    PE.JM2<-PE.AD.JM(jointFit2, newdata = test_data, Tstart = timeBase, Thoriz = timeBase + win,
                      idVar = "id", timeVar = "time", respVar = "y", evTimeVar = "Time",
                      lossFun = "square", summary = "value", simulate = TRUE)


    vec.pe.LM <-c(vec.pe.LM, PE.LM["prederr"])
    vec.pe.JM1 <-c(vec.pe.JM1, PE.JM1["prederr"])
    vec.pe.JM2 <-c(vec.pe.JM2, PE.JM2["prederr"])
    
  }
  vecs <- list(vec.pe.LM, vec.pe.JM1, vec.pe.JM2)
  vecs
}