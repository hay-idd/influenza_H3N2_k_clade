sigmoid_func <- function (x, k=0.1, x0=0, A=0) {
  y <- A + (1-A)/(1 + exp(-k * (x - x0)))
}


## Expand contact matrix for immunity classes
## For a given age-specific contact matrix and population size matrix (rows for age group, columns for immunity classes)
## Generates a normalised contact matrix giving contact rates between each age and immunity class combination
setup_C <- function(C1, Ns, beta_scales=rep(1, nrow=Ns)){
  Nimmunity <- ncol(Ns)
  Nage <- nrow(Ns)
  M <-  kron(C1,ones(Nimmunity,Nimmunity)) #' Non-normalised contact matrix scaled for each immunity class
  
  propns <- Ns/rowSums(Ns) #' Age/immunity population size as proportion of age population size
  propns_1 <- repmat(matrix(rep(propns, each=Nimmunity),
                            ncol = Nimmunity,byrow=FALSE),n=1,m=Nage)
  C <- M*propns_1 #' Generate scaled contact rates for age and immunity groups
  
  beta_scales_C <- repmat(matrix(rep(beta_scales, each=Nimmunity),ncol=Nimmunity*Nage,byrow=TRUE),n=Nage*Nimmunity,m=1)
  C <- C*beta_scales_C
  N_long <- c(t(Ns))
  C <- t(apply(C, 1, function(x) x/N_long))
  C[!is.finite(C)] <- 0
  #C <- t(C)
  
  return(C)
}
setup_C_explicit <- function(C1, Ns, beta_scales=rep(1, nrow=Ns)){
  Nimmunity <- ncol(Ns)
  Nage <- nrow(Ns)
  M <-  kron(C1,ones(Nimmunity,Nimmunity)) #' Non-normalised contact matrix scaled for each immunity class
  
  propns <- Ns/rowSums(Ns) #' Age/immunity population size as proportion of age population size
  propns_1 <- repmat(matrix(rep(propns, each=Nimmunity),
                            ncol = Nimmunity,byrow=FALSE),n=1,m=Nage)
  C <- M*propns_1 #' Generate scaled contact rates for age and immunity groups
  
  beta_scales_C <- repmat(matrix(rep(beta_scales, each=Nimmunity),ncol=Nimmunity*Nage,byrow=TRUE),n=Nage*Nimmunity,m=1)
  C <- C*beta_scales_C
  return(C)
}
## Get's the transmissibility parameter required for a desired R0
get_beta <- function(C, propns, inf_period=5, r0){
  fi <- kron(propns, ones(1, length(propns)))
  fj <- kron(t(propns), ones(length(propns),1))
  M <- C * fi/fj
  R0 <- max(eigen(M)$values)
  beta <- r0/(R0*inf_period)
  beta
}

## Get's the transmissibility parameter required for a desired R0
get_beta_vector <- function(C, propns, inf_period=5, r0, beta_scales){
  fi <- kron(propns, ones(1, length(propns)))
  fj <- kron(t(propns), ones(length(propns),1))
  C1 <- t(beta_scales * t(C))
  M <- C1 * fi/fj
  R0 <- max(eigen(M)$values)
  beta <- r0/(R0*inf_period)
  beta
}

get_sigmoid_func_outputs <- function(ages, pars){
  final <- matrix(nrow=nrow(pars),ncol=length(ages))
  for(j in 1:nrow(pars)){
    A <- pars$minimum_transmissibility_scale[j]
    k <- pars$transmissibility_scale[j]
    x0 <- pars$transmissibility_midpoint[j]
    y <- sigmoid_func(ages, k, x0, A)
    final[j,] <- y
  }
  colnames(final) <- ages
  final <- cbind(pars, final)
  final <- final %>% pivot_longer(as.character(ages),names_to="age_group",values_to="beta")
  return(final)
}

## Thins the polymod contact matrix for the desired age range
thin_polymod <- function(polymod, age_lower=0, age_upper=20, thin_propn=0.2,
                         same_age=TRUE){
  
  polymod1 <- polymod
  tmp_part <- polymod$participants
  tmp_contacts <- polymod$contacts
  tmp_contacts <- tmp_part %>% 
    select(part_id, part_age) %>% 
    right_join(tmp_contacts) %>% 
    mutate(row_id=1:n())
  ## Get contacts with age within range or overlapping range
  
  if(same_age){
    tmp_contacts_tothin <- tmp_contacts %>%
      filter((part_age <= age_upper & part_age >= age_lower) & 
               same_age & ((cnt_age_exact <= age_upper & cnt_age_exact >= age_lower) |
                           (cnt_age_est_min < age_upper & cnt_age_est_min > age_lower) |
                           (cnt_age_est_max > age_lower & cnt_age_est_max < age_upper)))
  } else {
    tmp_contacts_tothin <- tmp_contacts %>%
      filter((part_age <= age_upper & part_age >= age_lower))
  }
  tmp_contacts_subset <- tmp_contacts_tothin %>% 
    group_by(part_id) %>%
    sample_frac(1-thin_propn)
  
  tmp_contacts_rest <- tmp_contacts %>% filter(!(row_id %in% tmp_contacts_tothin$row_id))
  tmp_contacts <- bind_rows(tmp_contacts_subset, tmp_contacts_rest) %>% select(-part_age, -row_id)
  polymod1$contacts <- tmp_contacts
  polymod1
}

# function to compute weight vector given one holiday interval
weight_for_hol <- function(dates, start, end, smooth_time) {
  st <- as.Date(start); en <- as.Date(end)
  st_buf <- st - smooth_time
  en_buf <- en + smooth_time
  # normalized t in [0,1] for transition windows
  t_before <- pmin(pmax(as.numeric(dates - st_buf) / smooth_time, 0), 1)  # 0 at st_buf -> 1 at st
  t_after  <- pmin(pmax(as.numeric(dates - en) / smooth_time, 0), 1)      # 0 at en -> 1 at en+smooth_time
  # easing: eased(t) goes 0 -> 1 smoothly
  eased_before <- (1 - cos(pi * t_before)) / 2
  eased_after  <- (1 - cos(pi * t_after)) / 2
  # weight pieces:
  # - before window: weight = 1 -> during window (st_buf..st) weight transitions 1 -> 0
  w_before <- 1 - eased_before    # 1 -> 0 as t_before 0 -> 1
  # - inside holiday: weight = 0
  w_inside  <- rep(0, length(dates))
  # - after window: weight transitions 0 -> 1
  w_after <- eased_after          # 0 -> 1 as t_after 0 -> 1
  # assemble: for each date pick the correct piece
  w <- ifelse(dates >= st & dates <= en, w_inside,
              ifelse(dates >= st_buf & dates < st, w_before,
                     ifelse(dates > en & dates <= en_buf, w_after, 1)))
  as.numeric(w)
}


setup_holiday_tibble <- function(start_date,end_date, half_term_start,half_term_end,winter_holiday_start,winter_holiday_end,  smooth_time = 7){
  
  all_days <- tibble(date = seq(start_date, end_date, by = "day"))
  
  # label
  school_days <- all_days %>%
    mutate(
      label = case_when(
        date >= half_term_start   & date <= half_term_end   ~ "half_term",
        date >= winter_holiday_start & date <= winter_holiday_end ~ "winter_holiday",
        TRUE ~ "term_time"
      )
    )
  
  ## Smooths out contact matrix transitions around holidays
  
  ## Generate interpolations of contact matrix weightings around holidays
  # compute weights per holiday and take the minimum (closest holiday dominates)
  hols <- tibble(
    start = c(half_term_start, winter_holiday_start),
    end   = c(half_term_end, winter_holiday_end)
  )
  
  weight_list <- lapply(seq_len(nrow(hols)), function(i) {
    weight_for_hol(school_days$date, hols$start[i], hols$end[i], smooth_time)
  })
  
  # combine by taking the minimum weight across holidays (so any holiday proximity reduces weight)
  weights_combined <- Reduce(pmin, weight_list)
  
  # append to tibble
  school_days_weighted <- school_days %>%
    mutate(weight = weights_combined)
  school_days_weighted
}

setup_period_tibble <- function(start_date, end_date, ..., smooth_time = 7, period_names = NULL) {
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  pts   <- as.Date(unlist(list(...)))
  if (length(pts) == 0) return(tibble::tibble(date = dates, label = "term_time", weight = 1))
  starts <- pts[seq(1, length(pts), 2)]
  ends   <- pts[seq(2, length(pts), 2)]
  n      <- length(starts)
  if (is.null(period_names)) period_names <- paste0("period", seq_len(n))
  wlist  <- lapply(seq_len(n), function(i) weight_for_hol(dates, starts[i], ends[i], smooth_time))
  weight <- Reduce(pmin, wlist)
  label  <- rep("term_time", length(dates))
  for (i in seq_len(n)) label[dates >= starts[i] & dates <= ends[i] & label == "term_time"] <- period_names[i]
  tibble::tibble(date = dates, label = label, weight = weight)
}

aggregate_to_groups <- function(ret, age0 = NULL){
  cols <- c("[0,5)","[5,15)","[15,25)","[25,35)","[35,45)","[45,55)","[55,65)","[65,75)","75+")
  colnames(ret) <- cols
  #stopifnot(all(cols %in% colnames(ret)))
  a05 <- ret[ , "[0,5)"]
  # compute 1-4: subtract explicit age0 if supplied, otherwise assume uniform -> keep 4/5
  a1_4 <- if(!is.null(age0)) {
    stopifnot(length(age0) == nrow(ret))
    a05 - age0
  } else a05 * 4/5
  out <- cbind(
    `1-4`   = a1_4,
    `5-14`  = ret[ , "[5,15)"],
    `15-44` = ret[ , "[15,25)"] + ret[ , "[25,35)"] + ret[ , "[35,45)"],
    `45-64` = ret[ , "[45,55)"] + ret[ , "[55,65)"],
    `65+`   = ret[ , "[65,75)"] + ret[ , "75+"]
  )
  as.data.frame(out)
}
