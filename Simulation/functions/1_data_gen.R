
#function to generate data for a given scenario
data_gen <- function(scenario, n){  
  
  K <- 15  # number of planned repeated measurements per subject, per outcome
  t.max <- 15 # maximum follow-up time
  
  ## 1a) longitudinal covariate measurements (both scenarios)-------------------
  
  # define parameters:
  betas <- c("Intercept" = 2.0, "Time" = 1.25)
  sigma.y <- 0.5 # measurement error standard deviation
  D <- matrix(c(0.55, 0.2, 0.15, 0.45), 2, 2)
  D <- (D + t(D)) / 2
  
  # at which time points longitudinal measurements are supposed to be taken
  times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) 
  
  # group indicator, i.e., '0' placebo, '1' active treatment
  group <- rep(0:1, each = n/2) 
  DF <- data.frame(year = times, group = factor(rep(group, each = K))) #replicates the group indiciators so it is repeated for each measurement
  
  # design matrices for the longitudinal measurement model
  X <- model.matrix(~ year, data = DF)
  Z <- model.matrix(~ year, data = DF)
  
  #simulate random effects
  b <- mvrnorm(n, rep(0, nrow(D)), D) #n = no. of observations sampled from the multivariate normal
  
  # simulate longitudinal responses
  id <- rep(1:n, each = K) #list of IDS to match the observations [i.e. repeat each ID K times]
  eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
  y <- rnorm(n * K, eta.y, sigma.y) #simulate n*K values from a normal dist with mean vals eta.y and variance sigma.y
  
  #1b) event times - depends on scenario----------------------------------------
  
  # parameters in both scenarios:
  phi <- 0.95 # shape for the Weibull baseline hazard
  meanCens0 <- 10 # mean of the uniform censoring distribution for group 0
  meanCens1 <- 14 # mean of the uniform censoring distribution for group 1
  
  # design matrix for the survival model
  W <- cbind("(Intercept)" = 1, "Group" = group)
  
  #simulate event times depending on scenario
  trueTimes <- sim_events(scenario, n, phi, W, betas, b)
  
  #For any failed attempt at finding root - delete the data for that individual from all relevant vectors/matrices/dataframes
  na.ind <- !is.na(trueTimes) #list of length n, =TRUE if not NA, = FALSE if NA
  trueTimes <- trueTimes[na.ind] #if an element of na.ind=FALSE, the corresponding trueTimes element is deleted
  W <- W[na.ind, , drop = FALSE] #deletes the row corresponding to na.ind=FALSE element
  group <- group[na.ind] #if an element of na.ind=FALSE, the corresponding group element is deleted
  long.na.ind <- rep(na.ind, each = K) #extend the na.ind vector for vectors with K elements per individual
  y <- y[long.na.ind] #delete any y when corresponding na.ind=FALSE
  X <- X[long.na.ind, , drop = FALSE] #deletes the row corresponding to na.ind=FALSE element ((long)Kxn version)
  Z <- Z[long.na.ind, , drop = FALSE] #deletes the row corresponding to na.ind=FALSE element ((long)Kxn version)
  DF <- DF[long.na.ind, ] #deletes the row (from df) corresponding to na.ind=FALSE element ((long)Kxn version)
  n <- length(trueTimes) #redefine number of individuals
  
  #1c) simulate censoring times-------------------------------------------------
  Ctimes <- numeric(n)
  Ctimes[group == 0] <- runif(sum(group == 0), 0, 2 * meanCens0)
  Ctimes[group == 1] <- runif(sum(group == 1), 0, 2 * meanCens1)
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  
  #1d) edit longitudinal measurements depending on events-----------------------
  
  # keep the non missing cases, i.e., drop the longitudinal measurements
  # that were taken after the observed event time for each subject.
  ind <- times[long.na.ind] <= rep(Time, each = K) #for each element of times, set ind=TRUE if the observation time<event time and FALSE otherwise
  y <- y[ind] #keep only the longitudinal measurements when ind=TRUE (i.e. obs time <= event time)
  X <- X[ind, , drop = FALSE] #same for the design matrices
  Z <- Z[ind, , drop = FALSE] #" "
  id <- id[long.na.ind][ind] #same for id
  id <- match(id, unique(id)) #renames the ID labels so that it goes 1,1,1,2,2,3.. without skipping any integers [in case an individual has been removed]
  
  dat <- DF[ind, ] #new data frame with only the observations to keep 
  #add variable columns to dat
  dat$id <- id 
  dat$y <- y 
  dat$Time <- Time[id] #event times (repeated for each observation time)
  dat$event <- event[id] #event indicator (" ")
  names(dat) <- c("time", "group", "id", "y", "Time", "event")
  
  #completed data frame for simulated data
  dat
}



sim_events <- function(scenario, n, phi, W, betas, b){
  
  if(scenario==1){
    trueTimes <- sim_events_s1(n, phi, W, betas, b)
  }else if(scenario==2){
    trueTimes <- sim_events_s2(n, phi, W, betas, b)
  }else{
    print("not a valid scenario")
  }
  trueTimes
}

## SCENARIO 1

sim_events_s1 <- function(n, phi, W, betas, b){
  #scenario specific parameters-----------------------------
  gammas <- c("(Intercept)" = -5.0, "Group" = 0.5) # coefficients for baseline covariates
  alpha <- 0.5 # association parameter 
  
  eta.t <- as.vector(W %*% gammas) # gamma0 + gamma1*groupi
  
  u <- runif(n) #u is the simulated event probability u ~ Unif(0,1)
  
  trueTimes <- numeric(n) #vector of 'doubles' - each item = 0
  for (i in 1:n) {
    #print(i)
    Up <- 50
    tries <- 5
    #uniroot(f, interval[low, up], ) - searches the interval low - upper for a root (f=0)
    #u and i are the parameters needed for the function invS
    Root <- try(uniroot(invS1, interval = c(1e-05, Up), u = u[i], i = i, betas = betas, b = b, phi = phi, eta.t = eta.t, alpha = alpha)$root, TRUE)

    #result of uniroot - $root = location of root, $f.root = value of function at root
    while(inherits(Root, "try-error") && tries > 0) {
      #if we can't find the root, keep increasing the interval (for a max no. of tries)
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS1, interval = c(1e-05, Up), u = u[i], i = i, betas = betas, b = b, phi = phi, eta.t = eta.t, alpha = alpha)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  trueTimes
}

h1 <- function (s, i, betas, b, phi, eta.t, alpha) {
  #coefficients of beta [fixed effects] (simple linear model)
  XX <- cbind(1, s)
  #coefficients of b [random effects]
  ZZ <- cbind(1, s)
  f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ])) #mi(t)
  hres <- exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha) #phi is sigma_t in the paper
  hres
}

#invS = log(u) + int_0^t h(s) ds = 0 [need to solve for = 0]
invS1 <- function (t, u, i, betas, b, phi, eta.t, alpha) {
  int_res <- integrate(h1, lower = 0, upper = t, i=i, betas=betas, b=b, phi=phi, eta.t=eta.t, alpha=alpha)
  res <- int_res$value + log(u)
  res
}


## SCENARIO 2

sim_events_s2 <- function(n, phi, W, betas, b){
  gammas <- c("(Intercept)" = -7.0, "Group" = 0.25) # coefficients for baseline covariates
  alpha <- 0.25 # association parameter (changed from 0.5 to 0.25)
  
  # simulate event times
  eta.t <- as.vector(W %*% gammas) # gamma0 + gamma1*groupi
  
  u <- runif(n) #generate n random numbers from a unif dist between 0 and 1
  #u is the simulated event probability
  
  trueTimes <- numeric(n) #vector of 'doubles' - each item = 0
  for (i in 1:n) {
    #print(i)
    Up <- 50
    tries <- 5
    #u and i are the parameters needed for the function invS
    Root <- try(uniroot(invS2, interval = c(1e-05, Up), u = u[i], i = i, betas = betas, b = b, phi = phi, eta.t = eta.t, alpha = alpha)$root, TRUE)
    #result of uniroot - $root = location of root, $f.root = value of function at root
    while(inherits(Root, "try-error") && tries > 0) {
      #if we can't find the root, keep increasing the interval (for a max no. of tries)
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS2, interval = c(1e-05, Up), u = u[i], i = i, betas = betas, b = b, phi = phi, eta.t = eta.t, alpha = alpha)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  trueTimes
}


h2 <- function (s, i, betas, b, phi, eta.t, alpha) {
  #integrate mi(t') [fixed effects] from 0 to t where mi(t')_FE=beta0 + beta2*t' 
  #gives: beta0*t + beta2*t^2/2 
  #Now enter the coefficients of beta into XX (t=s):
  XX <- cbind(s, (s^2)/2)
  #integrate mi(t') [random effects] from 0 to t where mi(t')_RE=bi0+bi1*t'
  ZZ <- cbind(s, (s^2)/2)
  f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ])) #mi(t)
  hres <- exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha) #phi is sigma_t in the paper
  hres
}



#invS = log(u) + int_0^t h(s) ds = 0 [need to solve for = 0]
invS2 <- function (t, u, i, betas, b, phi, eta.t, alpha){
  int_res <- integrate(h2, lower = 0, upper = t, i=i, betas=betas, b=b, phi=phi, eta.t=eta.t, alpha=alpha)
  res <- int_res$value + log(u)
  res
}

# function to split the data into groups
create_split <- function(dat){
  dat.id <- dat[!duplicated(dat$id), ]
  n <- nrow(dat.id)
  #split data into 10 groups---------------------------------------
  splits <- split(seq_len(n), sample(rep(seq_len(10), length.out = n))) #group ids
  
  vecs <- list()
  for(g in 1:10){
    test_id <- unname(unlist(splits[g], recursive = TRUE, use.names = TRUE))
    vecs <- c(vecs, list(test_id))
  }
  list(dat, vecs)
}

