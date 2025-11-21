## As above, but each age group has its own transmission parameter
general_sir <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[2:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=4)
  dS <- -(alphas*sir[1,]) * (beta*sir[2,] %*% t(C))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  dInc <- (alphas*sir[1,]) * (beta*sir[2,] %*% t(C))

  tmp <- as.vector(rbind(dS,dI,dR,dInc))
  
  return(list(c(tmp)))
}

## As above, but each age group has its own transmission parameter
general_sir_timevarying_C <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[2:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=4)
  dS <- -(alphas*sir[1,]) * (beta*sir[2,] %*% t(C[[floor(t)]]))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  dInc <- (alphas*sir[1,]) * (beta*sir[2,] %*% t(C[[floor(t)]]))
  
  tmp <- as.vector(rbind(dS,dI,dR,dInc))
  
  return(list(c(tmp)))
}


## As above, but each age group has its own transmission parameter
general_sir_timevarying_C_symp <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  kappa <- pars[2]
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[3:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=5)
  dS <- -(alphas*sir[1,]) * (beta*sir[2,] %*% t(C[[floor(t)]]))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  dD <-  kappa*(alphas*sir[1,]) * (beta*sir[2,] %*% t(C[[floor(t)]]))
  dInc <- (alphas*sir[1,]) * (beta*sir[2,] %*% t(C[[floor(t)]]))
  
  tmp <- as.vector(rbind(dS,dI,dR,dD,dInc))
  
  return(list(c(tmp)))
}




## As above, but each age group has its own transmission parameter
general_sir_explicit <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[2:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=4)
  N <- colSums(sir)
  M <- t(apply(C, 1, function(x) x/N))
  M[!is.finite(M)] <- 0
  dS <- -(alphas*sir[1,]) * apply(M, 1, function(x) sum(beta*sir[2,]*x))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  dInc <- (alphas*sir[1,]) * apply(M, 1, function(x) sum(beta*sir[2,]*x))
  
  tmp <- as.vector(rbind(dS,dI,dR,dInc))
  
  return(list(c(tmp)))
}

#' Epidemic Final Size Calculation ODE
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class using an SIR model
#' @param C1 the normalised contact matrix of contact frequencies between each age and immunity class
#' @param beta the disease specific beta (transmission rate). Note that another parameter will mediate the contact rate
#' @param Tg number of days spent infectious (ie. 1/gamma)
#' @param Ns the matrix of population sizes for each age/immunity combination (non-normalised) (ie. rows = ages, cols = immunity classes)
#' @param alphas a vector of values between 0 and 1 matching the number of immunity classes
#' @param ts vector of times to solve over
#' @param beta_scales vector giving relative transmissibility (0-1) for each age class
#' @param age_seed which age group to seed in
#' @param immunity_seed which immunity class to seed in
#' @return an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
#' @seealso \code{\link{epi_final_size}}
#' @export
epi_ode_size <- function(C1, beta, Tg, Ns, alphas, kappa = NULL,initial_immune_frac=0,
                                ts=seq(1,365,by=1),
                                age_seeds=c(1),immunity_seed=1,seed_size=1,
                                ver="fast",return_peak=FALSE,return_compartments=FALSE){
  C <- C1
  Nage <- nrow(Ns)
  Nimmunity <- ncol(Ns)
  
  ## Generate starting populations.
  ## This generates a vector with S, I and R for each age class and immunity class
  ## ie. S_11, I_11, R_11, incidence_11, S_12, I_12, R_12, incidence_12, S_22, ... etc.
  ## where S_ak gives the number susceptible in age class 1, immunity class k
  if(is.null(kappa)){
    long_Ns <- as.numeric(t(Ns))
    start <- NULL
    start[1] <- long_Ns[1]
    start[2] <- 0
    start[3] <- 0
    start[4] <- 0
    index <- 5
    for(i in 2:length(long_Ns)){
      start[index] <- (1-initial_immune_frac)*long_Ns[i] ## Susceptible
      index <- index + 1
      start[index] <- 0 ## Infected
      index <- index + 1
      start[index] <- initial_immune_frac*long_Ns[i] ## Recovered
      index <- index + 1
      start[index] <- 0
      index <- index + 1
    }
    start[start < 0] <- 0
    
    
    ## Move one susceptible to the infected classes for the seed population
    for(age_seed in age_seeds){
      start[(age_seed-1)*4*Nimmunity + (immunity_seed-1)*4+1] <- start[(age_seed-1)*4*Nimmunity + (immunity_seed-1)*4+1] - seed_size
      start[(age_seed-1)*4*Nimmunity + (immunity_seed-1)*4+2] <- start[(age_seed-1)*4*Nimmunity + (immunity_seed-1)*4+2] + seed_size
    }
  } else {
    long_Ns <- as.numeric(t(Ns))
    start <- NULL
    start[1] <- long_Ns[1]
    start[2] <- 0
    start[3] <- 0
    start[4] <- 0
    start[5] <- 0
    index <- 6
    for(i in 2:length(long_Ns)){
      start[index] <- long_Ns[i]
      index <- index + 1
      start[index] <- 0
      index <- index + 1
      start[index] <- 0
      index <- index + 1
      start[index] <- 0
      index <- index + 1
      start[index] <- 0
      index <- index + 1
    }
    start[start < 0] <- 0
    
    ## Move one susceptible to the infected classes for the seed population
    for(age_seed in age_seeds){
      start[(age_seed-1)*5*Nimmunity + (immunity_seed-1)*5+1] <- start[(age_seed-1)*5*Nimmunity + (immunity_seed-1)*5+1] - seed_size
      start[(age_seed-1)*5*Nimmunity + (immunity_seed-1)*5+2] <- start[(age_seed-1)*5*Nimmunity + (immunity_seed-1)*5+2] + seed_size
    }
  }
  #start[(age_seed + (immunity_seed-1) - 1)*3 + 1] <- start[(age_seed + (immunity_seed-1) - 1)*3 + 1] - 1
  #  start[(age_seed + (immunity_seed-1) - 1)*3 + 2] <- start[(age_seed + (immunity_seed-1) - 1)*3 + 2] + 1
  
  #beta_scales <- rep(beta_scales, each=Nimmunity)
  if(ver == "fast" & class(C) != "list"){
    y <- ode(y=start,t=ts,func=general_sir, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
  } else if(class(C) == 'list') {
    if(is.null(kappa)){
      y <- ode(y=start,t=ts,func=general_sir_timevarying_C, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
    } else {
      y <- ode(y=start,t=ts,func=general_sir_timevarying_C_symp, parms=c(beta,kappa,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
    }
  } else {
    y <- ode(y=start,t=ts,func=general_sir_explicit, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
  }
  ## Pull out the solved model.
  #rt_ts <- compute_Rt_series(y, pars = c(beta,alphas,Tg), C = C_list, Nage = Nage, Nimmunity = Nimmunity)

  y <- as.data.frame(y)
  y <- y[,2:ncol(y)]
  ## Label which disease state each column is
  if(!is.null(kappa)){
    colnames(y) <- rep(c("S","I","R","D","inc"), Nage*Nimmunity)
    if(return_compartments){
      colnames(y) <- expand_grid(age=1:Nage,immunity=1:Nimmunity,compartment=c("S","I","R","D","inc")) %>% 
        mutate(name = paste0(compartment,"_",age,"_",immunity)) %>% pull(name)
      y$time <- ts
      #y <- left_join(y, rt_ts)
      return(y)
    }
  } else {
    colnames(y) <- rep(c("S","I","R","inc"), Nage*Nimmunity)
    if(return_compartments){
      colnames(y) <- expand_grid(age=1:Nage,immunity=1:Nimmunity,compartment=c("S","I","R","inc")) %>% 
        mutate(name = paste0(compartment,"_",age,"_",immunity)) %>% pull(name)
      y$time <- ts
      #y <- left_join(y, rt_ts)
      return(y)
    }
  }
  

  
  ## Pull out the recovered population to get final size
  recovered <- y[,which(colnames(y) == "R")]
  prevalence <- y[,which(colnames(y) == "I")]
  total_prev <- rowSums(prevalence)
  peak_time <- ts[which.max(total_prev)]
  ## Final size is number in recovered at end minus number recovered at start
  final_incidence <-  recovered[nrow(recovered),] - recovered[1,]
  final_incidence <- as.numeric(final_incidence/long_Ns)
  final_incidence[is.nan(final_incidence)] <- 0
  A <- matrix(final_incidence,nrow=Nage,ncol=Nimmunity,byrow=T)
  if(!return_peak){
    return(A)
  } else {
    return(list(A, peak_time))
  }
}

# compute Rt at a single time/state
compute_Rt_one <- function(t, y, pars, C, Nage, Nimmunity, Tg_index = length(pars)) {
  # pars: numeric vector as in your ODE: pars[1]=beta, pars[2:(end-1)]=alphas, last = Tg
  beta <- pars[1]
  Tg   <- pars[Tg_index]
  alphas <- pars[2:(Tg_index-1)]
  # reconstruct S and I from y (your ODE uses 4 rows x (Nage*Nimmunity) columns stored column-major)
  m <- Nage * Nimmunity
  # y is vector length 4*m in the same ordering as in ODE: rbind(dS,dI,dR,dInc) => rows of 4 blocks
  # but original state y passed to ODE was vectorised from rbind(S,I,R,Inc) too.
  S <- matrix(y, ncol = m, nrow = 4)[1, ]   # first row = S
  I <- matrix(y, ncol = m, nrow = 4)[2, ]   # second row = I   (I not used directly here except if want group Re)
  # pick contact matrix for this time: your ODE uses C[[floor(t)]]
  Cmat <- C[[floor(t)]]
  # build A = diag(alphas * S) %*% t(Cmat)
  diag_vec <- as.numeric(alphas * S)
  A <- diag_vec * t(Cmat)   # uses broadcast: each row j of t(C) multiplied by diag_vec[j]
  # K = beta * Tg * A
  K <- beta * Tg * A
  # spectral radius: largest absolute eigenvalue
  ev <- eigen(K, symmetric = FALSE, only.values = TRUE)$values
  rho <- max(Mod(ev))
  return(as.numeric(rho))
}

# vectorised: compute Rt series from ODE output (deSolve style)
# 'out' is a matrix/data.frame returned by deSolve::ode: first column 'time', other columns are state variables in same order used by ODE
compute_Rt_series <- function(out, pars, C, Nage, Nimmunity, time_col = 1, Tg_index = length(pars)) {
  times <- out[, time_col]
  Rt_vec <- numeric(length(times))
  # state columns: assume out columns after time are states matching length(y0) and ordered like your ODE expects
  state_cols <- seq_len(ncol(out))[-time_col]
  for (i in seq_along(times)) {
    yrow <- as.numeric(out[i, state_cols])
    Rt_vec[i] <- compute_Rt_one(times[i], yrow, pars, C, Nage, Nimmunity, Tg_index = Tg_index)
  }
  data.frame(time = times, Rt = Rt_vec)
}
