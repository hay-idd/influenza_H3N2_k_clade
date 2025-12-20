####function to generate daily incidence from weekly incidence data 
#enter the dataframe with date and cases as columns:
scale_weekly_to_daily <- function(weekly_data_frame,roll_d) {
  #assign the average form the weekly case data to each day, and take the rolling mean 
  
  daily_df <- weekly_data_frame %>%
    rowwise() %>%
    mutate(
      daily_dates = list(seq(date, by = "day", length.out = 7)),
      daily_cases = list(rep(cases/7, 7))
    ) %>%
    unnest(cols = c(daily_dates, daily_cases)) %>%
    select(daily_dates, daily_cases) %>%  # Keep only daily columns
    rename(date = daily_dates, cases = daily_cases)
  
  
  # Compute 7-day rolling mean of cases
  daily_df$rolling_mean <- rollmean(daily_df$cases, k = roll_d, fill = NA, align = "right")
  
  # Compute growth rate: ratio of current rolling mean to previous rolling mean
  daily_df$growth_rate <- daily_df$rolling_mean / lag(daily_df$rolling_mean)
  
  #Remove NA rows 
  daily_df <- na.omit(daily_df)
  
  return(daily_df)
  
  

}



rolling_mean_to_weekly<-function(daily_df,anchor_date){

#anchor_date <- as.Date("2011-09-12")  # custom start
weekly_df <- daily_df %>%
  mutate(
    # Calculate offset in days from anchor_date
    offset_days = as.numeric(date - anchor_date),
    week_index = floor(offset_days / 7),
    date = anchor_date + (week_index * 7)
  ) %>%
  group_by(date) %>%
  summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")
}



#remove the zeros and ones ##this
remove_tails<-function(df){
  new_df<-NULL
  df$season_label=as.factor(df$season_label)
  xx=levels(df$season_label)
  nn<-length(xx)
  
  for (i in 1:nn) {
    #get the data by season 
    sub_df<-subset(df,season_label==xx[i])
    #get the cumulative sum:
   #yy<-cumsum(sub_df$I)
  yy<-cumsum(sub_df$cases)
    #get the start of the epidemic:
    inx_start<-which.max(yy>1)
    
    
    dup_mask <- duplicated(yy) | duplicated(yy, fromLast = TRUE)
    
    if (any(dup_mask)) {
      zz <- max(yy[dup_mask])
      inx_end <- match(zz, yy)
    } else {
      # fallback: largest value in yy
      zz <- max(yy)
      inx_end <- which.max(yy)  # first index of global max
    }
    
    #remove the tails from the epidemic;
    sub_df_updated<-sub_df[inx_start:inx_end,]
    new_df<-rbind(new_df,sub_df_updated)
  }
  return(new_df)
}



#function to estimate Rt by season; 

estimate_Rt_by_season<-function(dat,window_size,dtout,mean_si=3.6,std_si=1.6){
  Rt_df_hhh<-NULL
  dat$season_label=as.factor(dat$season_label)
  dat$season_label <- droplevels(dat$season_label)
  xx=levels(dat$season_label)
  nn<-length(xx)
  
  for(i in 1:nn){
    print(xx[i])
    dat_hx<-subset(dat,season_label==xx[i])
    dat_h<- data.frame("date"=dat_hx$date,I = round(dat_hx$I))
    
    
    t_start <- seq(2, nrow(dat_h) - window_size)       # start of each window
    t_end <- t_start + window_size                   # end of each window
    t_midpoint <- (t_start + t_end)/2
    
    
    #using the orginal function
    res_total_hx <- estimate_R(
      incid = dat_h,
      dt = 1L,
      dt_out = dtout,
      method = "parametric_si",
      config = make_config(list(
        mean_si = mean_si,
        std_si = std_si,
        t_start = t_start,
        t_end = t_end
      ))
    )

    Rt_df_hx <- data.frame(
      date = dat_h$date[res_total_hx$R$t_end],
      season   = xx[i],
      Rt_mean = res_total_hx$R$`Mean(R)`,
      Rt_lower = res_total_hx$R$`Quantile.0.025(R)`,
      Rt_upper = res_total_hx$R$`Quantile.0.975(R)`
    )
    
    Rt_df_hhh<-rbind(Rt_df_hhh,Rt_df_hx)
    
  }
  
  return(Rt_df_hhh)
}


estimate_Rt_by_season_aggregated<-function(dat,window_size,dtout,mean_si=3.6,std_si=1.6){
  
  Rt_df_hhh<-NULL
  dat$season_label <- droplevels(dat$season_label)
  dat$season_label=as.factor(dat$season_label)
  xx=levels(dat$season_label)
  nn<-length(xx)
  
  for(i in 1:nn){
    dat_hx<-subset(dat,season_label==xx[i])
    dat_h<- data.frame("date"=dat_hx$date,I = round(dat_hx$I))
    
    
    t_start <- seq(2, nrow(dat_h) - window_size)       # start of each window
    t_end <- t_start + window_size                   # end of each window
    
    
    #using the orginal function in the CRAN version 
    # res_total_hx <- estimate_R(
    #   incid = dat_h,
    #   method = "parametric_si",
    #   config = make_config(list(
    #     mean_si = mean_si,
    #     std_si = std_si,
    #     t_start = t_start,
    #     t_end = t_end
    #   ))
    # )
    
    #using aggregated data funtion
    res_total_hx <- estimate_R(
      incid = dat_h$I,
      dt = 7L,
      dt_out = dtout,
      method = "parametric_si",
      config = make_config(list(
        mean_si = mean_si,
        std_si = std_si,
        t_start = t_start,
        t_end = t_end
      ))
    )
    
    
    Rt_df_hx <- data.frame(
      date = dat_h$date[res_total_hx$R$t_end],
      season   = xx[i],
      Rt_mean = res_total_hx$R$`Mean(R)`,
      Rt_lower = res_total_hx$R$`Quantile.0.025(R)`,
      Rt_upper = res_total_hx$R$`Quantile.0.975(R)`
    )
    
    Rt_df_hhh<-rbind(Rt_df_hhh,Rt_df_hx)
    
  }
  
  return(Rt_df_hhh)
}
