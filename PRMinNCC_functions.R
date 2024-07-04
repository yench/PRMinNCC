# Function to simulate competing risks data
simul.data = function(n, # cohort size
                      mctrl = NULL, # c(m1, m2): number of controls for cause 1 (m1) and cause 2 (m2) 
                      age.match = FALSE, # Does the NCC study match on continuous age (+/- 2 yr)?
                      max.t = 10, # maximum administrative censoring time
                      max.enter = 2, # maximum entry time 
                      Smax.censor = 0.9, # 1-probability of randomly censored at max.t
                      Smax1 = 0.90, Smax2 = 0.95, # 1-baseline probability of failing from cause 1 or 2 at max.t
                      rr11 = 2, rr12 = 2.5,  # hazard ratios associated with cause 1 (exp(beta11) and exp(beta12))
                      rr21 = 1.5, rr22 = 3,  # hazard ratios associated with cause 2 (exp(beta21) and exp(beta22))
                      mu = c(z1 = 0, z2 = 0), # mean of covariates z1 and z2
                      sig = c(sdz1 = 1, sdz2 = 1), # sd of covariates z1 and z2
                      alpha1 = 1, alpha2 = 1){ # alpha1=alpha2 => proportional risks
  
  Sig = matrix(0, nrow = 2, ncol = 2); diag(Sig) = sig^2
  Covmat = data.table(MASS::mvrnorm(n = n, mu = mu, Sigma = Sig))
  
  # Coefficients
  b10 = log(-log(Smax1)/(max.t^alpha1)); b11 = log(rr11); b12 = log(rr12)
  b20 = log(-log(Smax2)/(max.t^alpha2)); b21 = log(rr21); b22 = log(rr22)
  
    if (age.match==FALSE) {
    eta1 = with(Covmat, exp((b10 + b11*z1 + b12*z2)))
    eta2 = with(Covmat, exp((b20 + b21*z1 + b22*z2)))  
  } else if (age.match==TRUE){
    Covmat[ , age:= .(runif(n = n, min = 18, max = 65))]
    Covmat[ , z3 := scale(age)] # z3 = normalized age
    b13 = log(1.6); b23 = log(1.3) # beta13 and beta23
    eta1 = with(Covmat, exp((b10 + b11*z1 + b12*z2 + b13*z3)))
    eta2 = with(Covmat, exp((b20 + b21*z1 + b22*z2 + b23*z3))) 
  }
  
  # Survival data
  entrtime = runif(n, 0, max.enter)
  deadtime = rexp(n, rate = -log(Smax.censor)/max.t) # random censoring time
  censtime = pmin(max.t - entrtime, deadtime) # end of study
  failtime1 = rweibull(n, shape = alpha1, scale = eta1^(-1/alpha1)) # failure time due to cause 1
  failtime2 = rweibull(n, shape = alpha2, scale = eta2^(-1/alpha2)) # failure time due to cuase 2
  failtime = pmin(failtime1, failtime2)
  causetype = (failtime1<failtime2) + 2*(failtime1>=failtime2)
  eventime = pmin(failtime, censtime)
  ind.fail = ifelse(failtime <= censtime, causetype, 0)
  
  DATA = data.table(cbind(ID = 1:n, Covmat, eventime, ind.fail))
  
  if(!is.null(mctrl)){
    K = max(ind.fail)
    if(length(mctrl) == K){
      m = mctrl
    } else {
      m = rep(mctrl[1], K)
    }
    # At each failure time, get the (matched) risk set without the case
    if (age.match==FALSE){
      DATA[ , riskset := .(ifelse(ind.fail>0, lapply(ID, function(i){(which(eventime >= eventime[i] & ID != ID[i]))}), NA))]
    } else if (age.match==TRUE){
      DATA[ , riskset := .(ifelse(ind.fail>0, lapply(ID, function(i){(which(eventime >= eventime[i] & ID != ID[i] & abs(age[i] - age) <= 2))}), NA))]
    }
    DATA[ , nrisk := .(ifelse(ind.fail>0, lapply(ID, function(i){(sum(eventime >= eventime[i]))}), NA))]
    DATA[ , nrisk := as.integer(nrisk)]

    ph2data = list()
    for(k in 1:K){
      # Sample m[k] controls for type-k cases
      ph2idk = DATA[ind.fail == k, .(ID, riskset)]
      ph2idk[lengths(riskset)  > m[k], paste0("ctrl", 1:m[k]) := as.integer(lapply(.I, function(x){sample(riskset[[x]], m[k])}))]
      ph2idk[lengths(riskset) == m[k], paste0("ctrl", 1:m[k]) := as.integer(riskset)]
      setnames(ph2idk, "ID", "case")
      setorder(ph2idk, case)
      ph2idk[, matched.id := paste0(k, ",", 1:.N)]

      # Create NCC data
      uniidk = melt(ph2idk, id.vars = c("matched.id"), 
                    measure.vars = c("case", paste0("ctrl", 1:m[k])))
      uniidk[ , case := ifelse(variable == "case", 1, 0)]
      setnames(uniidk, "value", "ID")
      uniidk[ , variable := NULL]
      if (age.match==FALSE){
        ph2dak = uniidk[DATA, `:=`(z1 = z1, z2 = z2, eventime = eventime, ind.fail = ind.fail, nrisk = nrisk), on = .(ID)]
      } else if (age.match==TRUE){
        ph2dak = uniidk[DATA, `:=`(z1 = z1, z2 = z2, z3 = z3, age = age,
                                   eventime = eventime, ind.fail = ind.fail, nrisk = nrisk), on = .(ID)]
      }
      ph2dak[ , nrisk := max(nrisk, na.rm = TRUE), by = .(matched.id)]
      ph2dak[ , wgt := nrisk/(..m[..k] + 1)]
      ph2dak[ , nrisk := NULL]
      setorder(ph2dak, matched.id, -case)
      ph2data[[k]] = ph2dak
    }
    
    DATA[ , `:=` (riskset = NULL, nrisk = NULL)]
    setorder(DATA, ID)
    attr(DATA, "ph2data") = ph2data
  }
  return(DATA)
}

# Function to double-code data and set up cell-mean coding
repli.dat = function(dat){
  nrow = nrow(dat)
  n.cause = length(unique(dat$ind.fail))-1
  
  # Replicate data
  out = setDT(dat[rep(1:nrow, each = n.cause), ])
  out[ , cause:= as.factor(rep(1:n.cause, nrow))]
  
  if ("case" %in% colnames(dat)) {    # NCC data
    out[ , status := .(ifelse(ind.fail==cause & case==1, 1, 0))]
  } else {                            # Full cohort data
    out[ , status := .(ifelse(ind.fail==cause, 1, 0))]}

  # Cell-mean coding is optional, but saves some work for getting variances
  out[ , `:=` (z1.1 = z1, z1.2 = z1, z2.1 = z2, z2.2 = z2)]
  out[cause=="1", `:=` (z1.2 = 0, z2.2 = 0)] 
  out[cause=="2", `:=` (z1.1 = 0, z2.1 = 0)]
  
  return(out)
}

# Function to estimate baseline hazards
bhaz.wgt = function(coxph_object, cause = NULL, wgt = NULL, matched.id = NULL){ 
  coxdetail_object = coxph.detail(coxph_object, riskmat = TRUE)
  beta = coxph_object$coefficients
  Zbeta = coxdetail_object$x %*% beta
  failtime = coxdetail_object$y[ , ncol(coxdetail_object$y)-1]
  fail.ind = (coxdetail_object$y[ , ncol(coxdetail_object$y)] == 1)
  
  # Reorder objects that are not reordered by coxph
  sortorder = as.integer(rownames(Zbeta))
  cause = cause[sortorder]
  
  if(is.null(matched.id)){ # Full cohort Breslow estimator
    risk_mat = coxdetail_object$riskmat[sortorder, ]
    YexpbZ = sweep(risk_mat, 1, exp(Zbeta), `*`)
    
    dNt = matrix(0, nrow(risk_mat), ncol(risk_mat))
    evnt = apply(risk_mat, 1, function(x) { ifelse(all(x==0), 0, max(which(x == 1))) } )
    dNt[cbind(1:nrow(risk_mat), evnt)] = 1
    dNt = sweep(dNt, 1, coxdetail_object$y[,ncol(coxdetail_object$y)], `*`)
    
    Gt = colSums(YexpbZ)
    bhaz_vec = colSums(dNt)/Gt
    
  } else { # NCC Langholz-Borgan estimator
    matched.id = matched.id[sortorder]
    wgt = wgt[sortorder]
    if(!"cause2" %in% names(beta)){wgt = wgt*as.numeric(substr(matched.id, 1, 1)==cause)}
    bhaz_vec = 1/aggregate(exp(Zbeta)*wgt, list(matched.id), sum)$V1
  }
  
  table = data.table(failtime = failtime[fail.ind], bhaz = bhaz_vec, cause = cause[fail.ind])
  if ("cause2" %in% names(beta)){  # Model (2)
    table = table[order(failtime)]
    table[ , bHaz := cumsum(bhaz)] # Joint cumulative baseline
    if("t" %in% names(beta)){      # NCC time-dependent gamma
      table[ , bHaz.cause := ifelse(cause=="1", bHaz, 
                                    cumsum(bhaz*exp(beta["cause2"] + beta["t"]*failtime/10)))]
    } else {                       # Constant gamma
      table[ , bHaz.cause := ifelse(cause=="1", bHaz, bHaz*exp(beta["cause2"]))]
    }
  } else {
    table = table[order(cause, failtime)]
    table[ , bHaz := cumsum(bhaz), by = cause] # Cumulative baseline by cause
  }
  return(table)
}