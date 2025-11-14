#cleaning the digitised data:
library(dplyr)
library(ggplot2)
#set the path 
setwd("/Users/palahakoon/Documents/GitHub/influenza_H3N2_k_clade/data/rcgp_digitised/raw_data")



clean_digitised_data<-function(data,n_start_week){
  # Define start and end dates
  start_date <- as.Date("2019-10-07")
  end_date <- as.Date("2025-11-03")
  
  
  #get the decimal scales in the time axis:
  time_scale<-data$V1
  
  #get the y axis:
  y_axis<-data$V2
  
  #the week that start the data set is 
  new_start_date <- start_date + weeks(n_start_week)
  total_weeks <- as.numeric(difftime(end_date, new_start_date, units = "weeks"))
  
  # Create sequence of weekly dates
  week_dates <- seq(from = new_start_date, by = "week", length.out = total_weeks)
  
  # Scale decimals to week indices (proportional mapping)
  scaled_indices <- scales::rescale(time_scale, to = c(1, total_weeks))
  # Round to nearest integer for indexing
  week_indices <- round(scaled_indices)
  
  # Map decimals to proportional weeks
  mapped_dates <- week_dates[week_indices]
  
  all_data<-data.frame(mapped_dates,y_axis)
  
  # Remove duplicates, keep first occurrence, digitizer seems to be re-reading teh same point multiple times
  all_data <- all_data %>% distinct()
  
  #if values are negative, make them zero. I think the digitiser must have read the ticks in the figure 
  
  all_data$y_axis[all_data$y_axis<0]=0
  colnames(all_data)[colnames(all_data) == "mapped_dates"] <- "Date"
  
  return(all_data)
  
}

all_clean_dat<-NULL

#read the dataset
#less than 5 age group:%of positive samples by by viral strain
data_1<-read.csv("less_than_5_1.csv",header = F)
#the week that start the data set is 3: 
k<-3
clean_dat<-clean_digitised_data(data_1,k)
clean_dat$age_group<-rep("less than 5")
all_clean_dat<-rbind(all_clean_dat,clean_dat)

data_2<-read.csv("5_18_1.csv",header=F)
k=0
clean_dat_2<-clean_digitised_data(data_2,k)
clean_dat_2$age_group<-rep("5to18")
all_clean_dat<-rbind(all_clean_dat,clean_dat_2)

data_3<-read.csv("19_64_1.csv",header=F)
k<-3
clean_dat_3<-clean_digitised_data(data_3,k)
clean_dat_3$age_group<-rep("19to64")
all_clean_dat<-rbind(all_clean_dat,clean_dat_3)

data_4<-read.csv("65_over_1.csv",header=F)
k<-4
clean_dat_4<-clean_digitised_data(data_4,k)
clean_dat_4$age_group<-rep("65_over")
all_clean_dat<-rbind(all_clean_dat,clean_dat_4)

colnames(all_clean_dat)[colnames(all_clean_dat) == "y_axis"] <- "percent_positive_samples_H3"

ggplot(all_clean_dat,aes(x=Date,y=percent_positive_samples_H3))+
  geom_point()+
  facet_grid(~age_group)

#save the data
write_csv(all_clean_dat,"/Users/palahakoon/Documents/GitHub/influenza_H3N2_k_clade/data/rcgp_digitised/clean_data/percent_positive_samples_H3_by_age.csv")



#do the same for the number of samples:

all_clean_dat<-NULL

#read the dataset
#less than 5 age group:%of positive samples by by viral strain
data_1<-read.csv("less_than_5_2.csv",header = F)
#the week that start the data set is 3: 
k<-4
clean_dat<-clean_digitised_data(data_1,k)
clean_dat$age_group<-rep("less than 5")
all_clean_dat<-rbind(all_clean_dat,clean_dat)

data_2<-read.csv("5_18_2.csv",header=F)
k=1
clean_dat_2<-clean_digitised_data(data_2,k)
clean_dat_2$age_group<-rep("5to18")
all_clean_dat<-rbind(all_clean_dat,clean_dat_2)

data_3<-read.csv("19_64_2.csv",header=F)
k<-0
clean_dat_3<-clean_digitised_data(data_3,k)
clean_dat_3$age_group<-rep("19to64")
all_clean_dat<-rbind(all_clean_dat,clean_dat_3)

data_4<-read.csv("65_over_1.csv",header=F)
k<-5
clean_dat_4<-clean_digitised_data(data_4,k)
clean_dat_4$age_group<-rep("65_over")
all_clean_dat<-rbind(all_clean_dat,clean_dat_4)


colnames(all_clean_dat)[colnames(all_clean_dat) == "y_axis"] <- "number_of_influenza_samples"

ggplot(all_clean_dat,aes(x=Date,y=number_of_influenza_samples))+
  geom_point()+
  facet_grid(~age_group)


#save the data
write_csv(all_clean_dat,"/Users/palahakoon/Documents/GitHub/influenza_H3N2_k_clade/data/rcgp_digitised/clean_data/number_of_influenza_samples_by_age.csv")




