
library(EpiEstim)
library(epitrix)
library(distcrete)
library(ggplot2)

setwd("~/Documents/GitHub/recent_influenza/Untitled")

source("source_code_Rt_aggregated.R")
source("source_code_epiestim_r_function.R")
source("source_code_epiestim_processes.R")


#read the data:fluenet data 


flu_net_data<-read.csv("data/final/flunet_h3_cases_historic.csv")

#get the total sums:
incidence_weekly<-flu_net_data$flu_pos_sum

#define seriel interval, Mean SI = 3.6 days
#SD SI = 1.6 days
# Cowling et al., Epidemiology, 2009

mean_si <- 3.6 
std_si  <- 1.6 

method <- "parametric_si"
config <- make_config(list(mean_si = mean_si,
                           std_si = std_si))


dat<-data.frame("dataes"=flu_net_data$date,"I"=flu_net_data$flu_pos_sum)



res_total <- estimate_R(
  incid = dat$I,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si
  )),
  dt=7L
)

XX=res_total$R
plot(res_total)
# Plot the default graph
plot(res_total)



#this is taking long to run.
flu_net_data$H3_sum[is.na(flu_net_data$H3_sum)] <- 0

res_h3 <- estimate_R(
  incid = flu_net_data$H3_sum,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si
  )),
  dt=7L
)


#england cases: 
england_cases<-read.csv("data/final/england_h3_cases_historic.csv")

englan_cases<-estimate_R(
  incid = england_cases$Influenza.A.H3N2,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si
  )),
  dt=7L,
  iter=20L
)

plot(englan_cases)

total_cases<-rowSums(england_cases[,2:5])

total_eng_cases<-estimate_R(
  incid = total_cases,
  method = "parametric_si",
  config = make_config(list(
    mean_si = mean_si,
    std_si = std_si
  )),
  dt=7L,
  iter=20L
)

plot(total_eng_cases)

ggplot(englan_cases$R, aes(x = t_start, y = `Mean(R)`)) +
  geom_line(color = "blue") +
  geom_line(data = total_eng_cases$R,aes(x = t_start, y = `Mean(R)`), color="red")
  labs(x = "Time", y = "R", title = "Estimated R") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # <-- horizontal line
  theme_minimal() +
  theme(legend.position = "bottom")

# Create a ggplot
ggplot(englan_cases$R, aes(x = t_start, y = `Mean(R)`)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = `Quantile.0.025(R)`, ymax = `Quantile.0.975(R)`), alpha = 0.2) +
  labs(x = "Time", y = "R", title = "Estimated R") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # <-- horizontal line
  theme_minimal() +
  theme(legend.position = "bottom")



# Extract t_end from res_total
t_end <- res_total$R$t_end

years <- 2012:2025
week_sequence <- unlist(lapply(years, function(y) paste(y, seq(1, 52), sep = "-")))

