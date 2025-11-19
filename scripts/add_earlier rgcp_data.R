library(tidyverse)
library(data.table)
setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")

setwd("data/rcgp_digitised/raw_data_2/")

files <- list.files()

## Total samples
files_total_samples <- files[which(files %like% "total_samples_tested") ]
read_csv(files_total_samples[1])
total_samples <- lapply(files_total_samples,function(x) read_csv(x, col_names=FALSE))

all_dates <- tibble(x_floor=1:nrow(total_samples[[1]]) + 20,date=dates)
age_groups <- c("19-64","5-18","65+","<5")
dates <- seq(as.Date("2022-01-03"),as.Date("2023-12-25"),by="weeks")
for(i in 1:length(total_samples)){
  total_samples[[i]]$t <- 1:nrow(total_samples[[i]])
  total_samples[[i]]$date <- dates[total_samples[[i]]$t]
  total_samples[[i]]$age_group <- age_groups[i]
}

## Need to add entry for 2023-10-30 of about 60
total_samples[[3]] <- total_samples[[3]] %>% mutate(date = if_else(date >= "2023-10-30",date+7,date),
                                                    t = if_else(date >= "2023-10-30",t+1,t))
total_samples[[3]] <- bind_rows(total_samples[[3]], data.frame(X1=23.25,X2=60,t=96,date=as.Date("2023-10-30"),age_group=age_groups[3])) %>% arrange(t)
total_samples <- do.call("bind_rows",total_samples) %>% select(date,t,X2,age_group) %>% rename(total_samples=X2)
total_samples <- total_samples %>% mutate(total_samples=pmax(0,total_samples))
ggplot(total_samples) + geom_line(aes(x=date,y=total_samples,col=age_group))

## Total flu samples
files_flu_samples <- files[which(files %like% "samples_influenza") ]
flu_samples <- lapply(files_flu_samples,function(x) read_csv(x, col_names=FALSE))

## Expand x values by multiplying by 5, and then flooring and averaging duplicates
flu_samples <- lapply(flu_samples,function(x){
  x <- x %>%mutate(x_floor = floor(X1*5)) %>% group_by(x_floor) %>% arrange(X1) %>% group_by(x_floor) %>% summarize(X2=mean(X2))
  x <- x %>% mutate(X2 = pmax(X2,0))
  x
})


for(i in 1:length(flu_samples)){
  flu_samples[[i]] <- full_join(all_dates,flu_samples[[i]]) %>% mutate(X2 = if_else(is.na(X2),0,X2))
  flu_samples[[i]]$age_group <- age_groups[i]
}
flu_samples <- do.call("bind_rows",flu_samples) %>% select(date,X2,age_group) %>% rename(flu_samples=X2)
ggplot(flu_samples) + geom_line(aes(x=date,y=flu_samples,col=age_group))
comb <- left_join(total_samples,flu_samples)
comb <- comb %>% mutate(total_samples = if_else(total_samples == 0 & flu_samples > 0, 1, total_samples))
comb <- comb %>% 
  mutate(smooth_total_samples = zoo::rollmean(total_samples, k=3, fill=NA, align="center"),
         smooth_flu_samples = zoo::rollmean(flu_samples, k=3, fill=NA, align="center")) %>%
  mutate(prop_flu = flu_samples/total_samples,
         smooth_prop=smooth_flu_samples/smooth_total_samples) 
ggplot(comb) + 
  geom_line(aes(x=date,y=prop_flu,col="Raw"))+
  geom_line(aes(x=date,y=smooth_prop,col="Smooth")) +
  facet_wrap(~age_group)

comb <- comb %>% mutate(date = date + 7)
comb <- comb %>% mutate(age_group=if_else(age_group == "<5","0-4",age_group))
## Merge in ILI data
desired_age_groups <- c("0-4","5-18","19-64","65+")
ili_cases_comb_expanded_grouped <- combine_age_groups_ILI(ili_cases_comb_expanded,desired_age_groups)
final_2022_2023_dat <- ili_cases_comb_expanded_grouped %>% select(Year,Week,date,group,ILI,ILI_per_100k) %>% rename(age_group=group) %>% filter(date >= "2022-06-01", date <= "2023-06-01")  %>%
  left_join(comb %>% select(date, age_group,flu_samples,prop_flu,smooth_total_samples,smooth_flu_samples,smooth_prop) %>% mutate(date = date - 1)) %>% drop_na() %>% mutate(ILI_flu = ILI * smooth_prop)

setwd("~/Documents/GitHub/influenza_H3N2_k_clade/")
write_csv(final_2022_2023_dat,file="data/final/flu_2022_2023.csv")


