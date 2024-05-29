# Edited version of prederrJM (from JMbayes package) for JMbayes objects
# Code copied from https://github.com/drizopoulos/JMbayes/tree/master/R/prederrJM.JMbayes.R
# Code edited by A Davies (2021)
# Edits are labelled by 'AD edit' throughout
# where practical I have shown the original code in a comment labelled #orig
# Edits are made to ensure calculation of prediction error (PE) exactly matches PE Eq 
# (see Eq (4.26) in Davies, Coolen and Galla (2021), or, equivalently the prediction 
# error equation (no number label) on pg 34 in Rizopoulos, D. (2016). The R package 
# JMbayes for fitting joint models for longitudinal and time-to-event data using MCMC.
# Journal of Statistical Software72(7), 1-46)
# NB: Ti = event time, t= base time, u=prediction time, delta_i = event indicator
# List of edits:
# (1) sum over i: change from Ti > t to Ti >= t 
# (2) 1st term of sum: change from I(Ti>u) to I(Ti>=u)
# (3) 2nd term of sum: change from delta_i*I(Ti<=u) to delta_i*I(Ti<u)
# (4) Allow calculation of PE if there are no events in the interval [t,u] (the 
#     second term in the sum is zero)
#
# NB: These edits were made and used for models where interval = "FALSE", and 
# censor type is right censoring (i.e. is_counting = FALSE)
# The validity of the function for interval = "TRUE" or is_counting=TRUE has 
# not been verified. 

PE.AD.JM <- function (object, newdata, Tstart, Thoriz, lossFun = c("square", "absolute"), 
                               interval = FALSE, idVar = "id", simulate = FALSE, M = 100, ...) {
  if (!inherits(object, "JMbayes"))
    stop("Use only with 'JMbayes' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(newdata[[idVar]]))
    stop("'idVar' not in 'newdata.\n'")
  lossFun <- if (is.function(lossFun)) {
    lf <- lossFun
    match.fun(lossFun)
  } else {
    lf <- match.arg(lossFun)
    if (lf == "absolute") function (x) abs(x) else function (x) x*x
  }
  id <- newdata[[idVar]]
  id <- match(id, unique(id))
  TermsT <- object$Terms$termsT
  SurvT <- model.response(model.frame(TermsT, newdata)) 
  is_counting <- attr(SurvT, "type") == "counting"
  Time <- if (is_counting) {
    ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
  } else {
    SurvT[, 1]
  }
  timeVar <- object$timeVar
  #AD edit (1)------------------------------------------------------------------
  ## keep individuals whose event time is Ti >= t (Tstart)
  #newdata2 <- newdata[Time > Tstart, ] #orig
  newdata2 <- newdata[Time >= Tstart, ] #AD edit
  #-----------------------------------------------------------------------------
  print("got up to newdata2")
  if(nrow(newdata2)>0){ #AD edit (20/02/24)
    print("newdata2 > 0")
    SurvT <- model.response(model.frame(TermsT, newdata2))
    if (is_counting) {
      id2 <- newdata2[[idVar]]
      f <- factor(id2, levels = unique(id2))
      Time <- ave(SurvT[, 2], f, FUN = function (x) tail(x, 1))
      delta <- ave(SurvT[, 3], f, FUN = function (x) tail(x, 1))
    } else {
      Time <- SurvT[, 1]
      delta <- SurvT[, 2]
    }
    #AD edit ---------------------------------------------------------------------
    ## keep observations <= t
    timesInd <- newdata2[[timeVar]] <= Tstart #orig
    
    ## Edit (2): alive individuals: Ti >= u 
    #aliveThoriz <- newdata2[Time > Thoriz & timesInd, ] #orig
    aliveThoriz <- newdata2[Time >= Thoriz & timesInd, ] #edit
    ## Edit (3): dead individuals: Ti < u and delta = 1
    #deadThoriz <- newdata2[Time <= Thoriz & delta == 1 & timesInd, ] #orig
    deadThoriz <- newdata2[Time < Thoriz & delta == 1 & timesInd, ] #edit
    
    ## censored individuals: Ti < u and delta = 0
    indCens <- Time < Thoriz & delta == 0 & timesInd #orig
    censThoriz <- newdata2[indCens, ] #orig
    #-----------------------------------------------------------------------------
    nr <- length(unique(newdata2[[idVar]]))
    idalive <- unique(aliveThoriz[[idVar]])
    iddead <- unique(deadThoriz[[idVar]])
    idcens <- unique(censThoriz[[idVar]])
    #AD edit (4)------------------------------------------------------------------
    #prederr <- if (length(unique(Time)) > 1 && nrow(aliveThoriz) > 1 &&  #orig
    #               nrow(deadThoriz) > 1) {                               #orig
    
    ## edit so that we don't get NA if no-one dies or is censored in the interval t->u
    
    ## work out sum if there are people alive after t (and t not equal to u)
    prederr <- if (length(unique(Time)) >= 1 && Tstart!=Thoriz) {
      ## work out term 1 (alive) if there are individuals alive after u
      if(nrow(aliveThoriz) >= 1){
        #orig:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Surv.aliveThoriz <- if (is_counting) {
          survfitJM(object, newdata = aliveThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = rep(Tstart, length(idalive)),
                    LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
          survfitJM(object, newdata = aliveThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = rep(Tstart, length(idalive)))
          
        }
        Surv.aliveThoriz <- sapply(Surv.aliveThoriz$summaries, "[", 2)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      } else{
        ## no-one alive after u
        Surv.aliveThoriz <- NA
      }
      
      ## work out term 2 (event) if there are individuals who died between t and u
      if(nrow(deadThoriz) >= 1){
        #orig:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Surv.deadThoriz <- if (is_counting) {
          survfitJM(object, newdata = deadThoriz, idVar = idVar, simulate = simulate,
                    survTimes = Thoriz, last.time = rep(Tstart, length(iddead)),
                    LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
          survfitJM(object, newdata = deadThoriz, idVar = idVar, simulate = simulate,
                    survTimes = Thoriz, last.time = rep(Tstart, length(iddead)))
        }
        Surv.deadThoriz <- sapply(Surv.deadThoriz$summaries, "[", 2)
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      } else{
        ## no-one died between t and u
        Surv.deadThoriz <- NA
      }
      ## work out term 3 (censor) if there are individuals who were censored between t and u
      if (nrow(censThoriz)>=1) {
    #---------------------------------------------------------------------------
        Surv.censThoriz <- if (is_counting) {
          survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = rep(Tstart, length(idcens)),
                    LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
          survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = rep(Tstart, length(idcens)))
        }
        tt <- Time[indCens]
        weights <- if (is_counting) {
          survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])],
                    LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
          survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, M = M,
                    survTimes = Thoriz, last.time = tt[!duplicated(censThoriz[[idVar]])])
        }
        Surv.censThoriz <- sapply(Surv.censThoriz$summaries, "[", 2)
        weights <- sapply(weights$summaries, "[", 2)
      } else {
        Surv.censThoriz <- weights <- NA
      }
      if (!interval) {
        (1/nr) * sum(lossFun(1 - Surv.aliveThoriz), lossFun(0 - Surv.deadThoriz),
                     weights * lossFun(1 - Surv.censThoriz) + (1 - weights) * lossFun(0 - Surv.censThoriz), 
                     na.rm = TRUE)
      } else {
        TimeCens <- object$y$Time
        deltaCens <- 1 - object$y$event
        KMcens <- survfit(Surv(TimeCens, deltaCens) ~ 1)
        times <- TimeCens[TimeCens > Tstart & TimeCens < Thoriz & !deltaCens]
        times <- sort(unique(times))
        k <- as.numeric(table(times))
        w <- summary(KMcens, times = Tstart)$surv / summary(KMcens, times = times)$surv
        prederr.times <- sapply(times, 
                                function (t) prederrJM(object, newdata, Tstart, t,
                                                       interval = FALSE, idVar = idVar, simulate = simulate)$prederr)
        num <- sum(prederr.times * w * k, na.rm = TRUE)
        den <- sum(w * k, na.rm = TRUE)
        num / den
      }
    } else {
      nr <- NA
      NA
    }
    out <- list(prederr = prederr, nr = nr, Tstart = Tstart, Thoriz = Thoriz, 
                interval = interval, classObject = class(object), 
                nameObject = deparse(substitute(object)), lossFun = lf)
  }else{
    print("newdata2 < 0")
    out <- list(prederr = 0, nr = 0, Tstart = Tstart, Thoriz = Thoriz, 
                interval = interval, classObject = class(object), 
                nameObject = deparse(substitute(object)), lossFun = lf)
  }
  class(out) <- "prederrJM"
  out
}