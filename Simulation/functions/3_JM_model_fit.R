fit_JM1 <- function(data, data.id){
  longFit_lin <- lme(y ~ time, data = data, random = ~ time | id)
  survFit <- coxph(Surv(Time, event) ~ group, data = data.id, x = TRUE)
  jointFit1 <- jointModelBayes(longFit_lin, survFit, timeVar = "time", n.iter = 100000L)
  models <- list(jointFit1, longFit_lin, survFit)
  models
}

fit_JM2 <- function(models){
  jointFit1 <- models[[1]]
  longFit_lin <- models[[2]]
  survFit <- models[[3]]

  iForm <- list(fixed = ~ 0 + time + I(time^2/2),
                indFixed = 1:2, random = ~ 0 + time + I(time^2/2), 
                indRandom = 1:2)
  jointFit2 <- update(jointFit1, param = "td-extra", extraForm = iForm)
  jointFit2
}

JM_params <- function(jointFit, model){
  
  #alpha
  if(model==1){
    alpha <- unname(jointFit$postMeans$alphas['Assoct'])
  }else if(model==2){
    alpha <- unname(jointFit$postMeans$Dalphas['AssoctE'])
  }else{
    print("Model must be 1 or 2")
  }
  
  #gamma
  gamma <- unname(jointFit$postMeans$gammas['group1'])
  
  params <- list(alpha, gamma)
  params
  
}
