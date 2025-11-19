

my_prior <- function(pars){
  a <- dbeta(pars["prop_immune_youngest"],1,5,1)
  b <- dbeta(pars["prop_immune_younger"],2,1.5,1)
  c <- dbeta(pars["prop_immune_older"],10,6,1)
  d <- dbeta(pars["prop_immune_oldest"],10,6,1)
  e <- dbeta(pars["alpha1"],1,1,1)
  f <- dunif(pars["reporting_rate"], 0, 0.01) #dbeta(pars["reporting_rate"],2,1000,1)
  g <- dnorm(pars["seed_size"],log(50),2,1)
  h <- dnorm(pars["R0"],rlnorm(1.8),0.25,1)
  i <- dlnorm(pars["obs_sd"],rlnorm(3),0.5)

  tot_prior <- a + b + c + d + e + f + g + h
  return(tot_prior)
}


posterior_func <- function(parTab, data, PRIOR_FUNC, 
                           #symp_frac,  ## Vector for symptomatic fraction by age group
                           N_props,  ## Proportion in each age group from polymod object
                           C_term, ## Contact matrix for term time, used to calibrate R0
                           C_list, ## List of contact matrices over time
                           return_dat=FALSE,...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################
  
  ## Somewhat redundant example
  parameter_names <- parTab$names
  
  tmp_dat <- flu_dat_wide %>% select(-date) %>% as.matrix()
  
  ## Vector of symptomatic fraction by age group
  
  cost_function <- function(pars){
    names(pars) <- parameter_names
    ## Pull out model parameters
    prop_immune_youngest <- pars[1]
    prop_immune_younger <- pars[2]
    prop_immune_older <- pars[3]
    prop_immune_oldest <- pars[4]
    alphas <- c(1,pars[5])
    reporting_rate <- pars[6]
    seed_size <- exp(pars[7])
    obs_sd <- pars[8]
    R0 <- pars[9]
    symp_frac <- c(pars[10:13])
    #print(R0)
    ## Get proporiton immune vector
    prop_immune <- c(prop_immune_youngest,prop_immune_younger, prop_immune_older,prop_immune_oldest)
    
    ## Housekeeping for setting up contact matrix
    N_age_classes <- length(N_props)
    N_immunity_classes <- 2
    
    N_props_long <- c(N_props*(1-prop_immune), N_props*(prop_immune))  ## Number of people in each age group and age class
    #N_props_long <- N_props
    N <- matrix(N_props_long*N_tot,ncol=length(alphas),nrow=N_age_classes)
    
    ## Get beta parameter
    beta_par <- get_beta(C_term,N_props,gamma, R0)
    beta_scales <- rep(1,N_age_classes)
    y_base <- epi_ode_size(C_list, beta_par, gamma, N, ts=ts,
                           alphas=alphas, age_seed=3,immunity_seed=1,seed_size=seed_size,return_compartments=TRUE)
    use_cols <- which(colnames(y_base) %like% "inc")
    
    
    ret <- as.matrix(y_base[,use_cols])
    
    ## Get daily incidence
    ret <- apply(ret, 2, function(x) c(0, diff(x)))
    
    ## Get expected symptomatic cases scaled by reporting rate
    # sum every 2 columns for the two immune states
    ret2 <- ret[, seq(1, ncol(ret), 2)] + ret[, seq(2, ncol(ret), 2)]
    ret2 <- t(reporting_rate*symp_frac*t(ret2))
    
    # then sum every 7 rows to get weekly data
    row_groups <- gl(nrow(ret2) %/% 7 + (nrow(ret2) %% 7 > 0), 7, nrow(ret2))
    ret_sum <- aggregate(ret2, by = list(row_groups), FUN = sum)[, -1]
    
    ## Final size 
    #y_base <- y_base %>% pivot_longer(-time)
    #y_base <- y_base %>% mutate(compartment = str_split(name, "_", simplify=TRUE)[,1],
    #                            age = as.integer(str_split(name, "_", simplify=TRUE)[,2]),
    #                            immunity = as.integer(str_split(name, "_", simplify=TRUE)[,3]))
    #total_inf <- y_base %>% filter(time == max(time),compartment == "inc")  %>% pull(value) %>% sum()
    
    #final_size <- total_inf/sum(N)
    
    if(return_dat){
      #if(aggregate){
      #  ret_sum <- aggregate_to_groups(ret_sum)
      #}
      return(ret_sum)
    } else {
      #ret_test <- aggregate_to_groups(ret_sum)
      #lik <- sum(dpois(round(c(tmp_dat)),c(as.matrix(ret_sum)),log=TRUE)) #+ dnorm(final_size,0.15,0.01,log=TRUE)
      lik <- sum(dnbinom(round(c(tmp_dat)),mu=c(as.matrix(ret_sum)),size=obs_sd,log=TRUE))
      if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
      
      return(lik)
    }
  }
  return(cost_function)
}
