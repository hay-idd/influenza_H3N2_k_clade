library(tidyverse)
library(ggplot2)

files <- list.files("data/RCGP/")

files <- files[files != "pos_rcgp.csv"]
all_dat <-NULL
index <- 1
for(file in files){
  all_dat[[index]] <- read_csv(paste0("data/RCGP/",file),col_names=FALSE) %>% mutate(age = file)
  index <- index + 1
}


all_dat <- do.call('bind_rows',all_dat)
all_dat$X1 <- round(all_dat$X1)
all_dat$X2 <- round(all_dat$X2,2)

colnames(all_dat) <- c("week","inc","age")
all_dat$inc <- pmax(all_dat$inc,0.01)
age_key <- c("ili_1-4yr.csv"="1-4 years", "ili_15-64yr.csv"="15-64 years", "ili_1yr.csv"="<1 year", "ili_5-14yr.csv"="5-14 years", 
             "ili_65+.csv" = "65+ years", "ili_rcgp.csv"="Combined")

all_dat$age <- age_key[all_dat$age]
all_dat$age <- factor(all_dat$age,levels=c("<1 year", "1-4 years",  "5-14 years", "15-64 years","65+ years", 
                                           "Combined"))
all_dat <- all_dat %>% group_by(age) %>% mutate(gr = log(inc/lag(inc)))
all_dat <- all_dat %>% group_by(age) %>% mutate(smooth_gr = zoo::rollmean(gr,3,fill=NA,align="right"))

## Convert dates
all_dat$date <- NA
all_dat <- all_dat %>% mutate(date = if_else(week == 44, "2025-10-25",NA))
all_dat$date <- ymd(all_dat$date)
all_dat$date <- ymd("2025-10-25") - (44-all_dat$week)*7 + 7


ggplot(all_dat%>% filter(date >= "2025-08-01")) + 
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+
  geom_label(x=ymd("2025-10-01"),y=15,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=15,label="Half term",size=5)+
  scale_y_continuous(limits=c(0,17)) +
  geom_point(aes(x=date,y=inc,col=age))+ 
  geom_line(aes(x=date,y=inc,col=age),linewidth=0.75)+ 
  scale_color_viridis_d() +
  scale_x_date(breaks="1 month") +
  theme_psi() +
  ylab("Weekly ILI per 100,000") +
  xlab("Date (end of reporting week)")

ggplot(all_dat%>% filter(date >= "2025-08-01")) + 
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+
  geom_label(x=ymd("2025-10-01"),y=log(15),label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=log(15),label="Half term",size=5)+
  geom_point(aes(x=date,y=log(inc),col=age))+ 
  geom_line(aes(x=date,y=log(inc),col=age),linewidth=0.75)+ 
  scale_color_viridis_d() +
  scale_x_date(breaks="1 month") +
  theme_psi() +
  ylab("Log weekly ILI per 100,000") +
  xlab("Date (end of reporting week)")

ggplot(all_dat %>% filter(date >= "2025-08-01")) + 
  geom_hline(yintercept=0,linetype="dashed") +
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+
  geom_label(x=ymd("2025-10-01"),y=2,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=2,label="Half term",size=5)+
  scale_color_viridis_d() +
  geom_point(aes(x=date,y=smooth_gr,col=age)) + 
  geom_line(aes(x=date,y=smooth_gr,col=age),linewidth=0.75) + 
  scale_x_date(breaks="1 month")+
  theme_psi() +
  ylab("Log growth rate of weekly ILI per 100,000") +
  xlab("Date (end of reporting week)")

## % flu pos from resp datamart
flu_pos <- read_csv("data/respmart_flu_pos_by_age.csv") %>% filter(Season == "2025 to 2026") %>% rename(week = `Week number`, date = Date) %>% select(-Season) %>%
  mutate(date = dmy(date)) %>%
  pivot_longer(-c(date,week)) %>% rename(age = name, pos = value) %>% mutate(pos = pos/100)

age_key2 <- c("Up to 5 years"="0-4 years", "5 to 14 years"="5-14 years", "15 to 44 years"="15-44 years", "45 to 64 years"="45-65 years", 
              "65 to 79 years"="65-79 years", "80 years and above"="80+ years")
flu_pos$age <- age_key2[flu_pos$age]


## Need to combine the different age groups
age_key_new <- c("0-4 years"="0-4 years",
             "1-4 years"="0-4 years",
             "<1 year"="0-4 years", 
             "5-14 years"="5-14 years", 
             "15-44 years"="15-64 years", 
             "45-65 years"="15-64 years", 
             "15-64 years"="15-64 years",
             "65-79 years"="65+ years", 
             "65+ years"="65+ years",
             "80+ years"="65+ years",
             "Combined" = "Combined")
flu_pos$age_new <- age_key_new[as.character(flu_pos$age)]
flu_pos <- flu_pos %>% group_by(age_new,week) %>% summarize(pos = mean(pos), date = first(date))

all_dat$age_new <- age_key_new[as.character(all_dat$age)]
all_dat <- all_dat %>% group_by(age_new,week) %>% summarize(inc = mean(inc), date = first(date))

all_dat_comb <- left_join(all_dat %>% select(-date), flu_pos %>% select(-date))

all_dat_comb <- all_dat_comb %>% mutate(ILI_plus = inc*pos)


all_dat_comb$date <- ymd("2025-10-25") - (44-all_dat$week)*7 + 7

all_dat_comb <- all_dat_comb %>% group_by(age_new) %>% mutate(gr = log(ILI_plus/lag(ILI_plus)))
all_dat_comb <- all_dat_comb %>% group_by(age_new) %>% mutate(smooth_gr = zoo::rollmean(gr,3,fill=NA,align="right"))
all_dat_comb <- all_dat_comb %>% filter(age_new != "Combined")

all_dat_comb$age_new <- factor(all_dat_comb$age_new, levels=c("0-4 years","5-14 years","15-64 years","65+ years"))
## ILI+

ggplot(all_dat_comb %>% filter(date >= "2025-08-01")) +
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+  
  geom_label(x=ymd("2025-10-01"),y=2.2,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=2.2,label="Half term",size=5)+
  scale_y_continuous(limits=c(0,2.3)) +
  geom_point(aes(x=date,y=ILI_plus,col=age_new)) +
  geom_point(aes(x=date,y=ILI_plus,col=age_new)) +
  geom_line(aes(x=date,y=ILI_plus,col=age_new),linewidth=0.75) +
  ylab("Weekly ILI+ (% pos * ILI)")+
  xlab("Date (end of reporting week)") +
  scale_color_viridis_d("") +
  theme_psi()

ggplot(all_dat_comb %>% filter(date >= "2025-08-01")) +
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+  
  geom_label(x=ymd("2025-10-01"),y=1.1,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=1.1,label="Half term",size=5)+
  scale_y_continuous(limits=c(-8,1.2)) +
  geom_point(aes(x=date,y=log(ILI_plus),col=age_new)) +
  geom_line(aes(x=date,y=log(ILI_plus),col=age_new),linewidth=0.75) +
  ylab("Log weekly ILI+ (% pos * ILI)")+
  xlab("Date (end of reporting week)") +
  scale_color_viridis_d("") +
  theme_psi()


ggplot(all_dat_comb%>% filter(date >= "2025-08-01")) +
  geom_hline(yintercept=0,linetype="dashed") +
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+  
  geom_label(x=ymd("2025-10-01"),y=2,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=2,label="Half term",size=5)+
  #geom_point(aes(x=date,y=gr,col=age_new),alpha=0.5) +
  geom_line(aes(x=date,y=smooth_gr,col=age_new),linewidth=0.75) +
  #geom_line(aes(x=date,y=gr,col=age_new),linewidth=0.75,alpha=0.5) +
  #geom_smooth(aes(x=date,y=gr,col=age_new),linewidth=1,span=0.5,se=FALSE) +
 # coord_cartesian(ylim=c(-2,2))+
  ylab("Log weekly growth rate of ILI+ (% pos * ILI)") +
  scale_color_viridis_d("") +
  xlab("Date (end of reporting week)") +
  theme_psi()
 #geom_line(aes(x=date + 7,y=gr,col=age_new),alpha=0.5) +
  #geom_smooth(aes(x=date + 7,y=gr,col=age_new),method="loess",se=FALSE,span=0.5)


ggplot(all_dat_comb%>% filter(date >= "2025-08-01")) +
  geom_hline(yintercept=0,linetype="dashed") +
  annotate("rect",xmin=ymd("2025-09-02"),xmax=ymd("2025-10-24"),ymin=-Inf,ymax=Inf,fill="red",alpha=0.25)+
  annotate("rect",xmin=ymd("2025-10-24"),xmax=ymd("2025-11-03"),ymin=-Inf,ymax=Inf,fill="green",alpha=0.25)+  
  geom_label(x=ymd("2025-10-01"),y=2,label="Autumn term",size=5)+
  geom_label(x=ymd("2025-11-01"),y=2,label="Half term",size=5)+
  geom_point(aes(x=date,y=gr,col=age_new)) +
 # geom_line(aes(x=date,y=smooth_gr,col=age_new),linewidth=0.75) +
  geom_line(aes(x=date,y=gr,col=age_new),linewidth=0.75) +
  #geom_smooth(aes(x=date,y=gr,col=age_new),linewidth=1,span=0.5,se=FALSE) +
  # coord_cartesian(ylim=c(-2,2))+
  ylab("Log weekly growth rate of ILI+ (% pos * ILI)") +
  scale_color_viridis_d("") +
  xlab("Date (end of reporting week)") +
  theme_psi()
