flu_cumu_inc <- 0.05
non_flu_cumu_inc <- 0.1
non_flu_cumu_inc_plateau <- 0.05
vacc_prob_overall <- 0.25
ve <- 0.5
tmax <- 50
growth_rate <- 0.05
N <- 50000
## Declining non-flu ILI
non_flu <- exp(-growth_rate * 0:tmax)
non_flu <- (non_flu_cumu_inc-non_flu_cumu_inc_plateau)*non_flu/sum(non_flu) + non_flu_cumu_inc_plateau/(tmax+1)

## Growing flu
flu <- exp(growth_rate * 0:tmax)
flu <- flu_ratio * flu/sum(flu)
flu_unvacc <- flu_cumu_inc * flu/sum(flu)
flu_vacc <- flu_cumu_inc * flu * (1-ve)/sum(flu)
flu_probs <- cbind(flu_unvacc,flu_vacc)

plot(flu_unvacc,type='l',col="blue",ylim=c(0,0.01))
lines(flu_vacc,col="green")
lines(non_flu,col='red')

## Vaccine coverage
vacc_prob <- rep(vacc_prob_overall/(tmax+1),tmax+1)
plot(vacc_prob)

## Run sim
vaccination_times <- sample(c(-1,0:tmax),N,prob=c(1-vacc_prob_overall,vacc_prob),replace=TRUE)
vacc_times_dat <- data.frame(id=1:N,vacc_time = vaccination_times,vacc_status=1-(vaccination_times == -1))
vaccination_states <- matrix(1, nrow=N, ncol=tmax+1)
## Simulate flu infections
for(i in 1:N){
  if(vaccination_times[i] == -1){
    vaccination_states[i,] <- rep(1, tmax+1)
  } else {
    vaccination_states[i,] <- c(rep(1,vaccination_times[i]),rep(2,tmax - vaccination_times[i] +1))
  }
}

inf_states_flu <- matrix(0, nrow=N, ncol=tmax+1)
for(i in 1:N){
  vacc_states <- vaccination_states[i,]
  for(t in 0:tmax){
    inf_states_flu[i,t+1] <- rbinom(1,1,flu_probs[t+1,vacc_states[t+1]])
  }
}
inf_states_flu <- reshape2::melt(inf_states_flu)
colnames(inf_states_flu) <- c("id","flu_time","flu_inf")
inf_states_flu <- inf_states_flu[inf_states_flu$flu_inf==1,]

inf_states_non_flu <- matrix(0, nrow=N,ncol=tmax+1)
for(i in 1:N){
  for(t in 0:tmax){
    inf_states_non_flu[i,t+1] <- rbinom(1,1,non_flu[t+1])
  }
}
inf_states_non_flu <- reshape2::melt(inf_states_non_flu)
colnames(inf_states_non_flu) <- c("id","non_flu_time","non_flu_inf")
inf_states_non_flu <- inf_states_non_flu[inf_states_non_flu$non_flu_inf == 1,]

inf_states_flu <- left_join(inf_states_flu, vacc_times_dat) %>% mutate(vacc_status_at_inf = ifelse(flu_time >= vacc_time & vacc_time != -1, 1,0)) %>% as_tibble()
inf_states_non_flu <- left_join(inf_states_non_flu, vacc_times_dat) %>% mutate(vacc_status_at_inf = ifelse(non_flu_time >= vacc_time & vacc_time != -1, 1,0)) %>% as_tibble()

## Influenza cases non-vacc vs. vacc
ab <- inf_states_flu %>% pull(vacc_status_at_inf) %>% table()

## Non influenza cases non-vacc vs. vacc
cd <- inf_states_non_flu %>% pull(vacc_status_at_inf) %>% table()

a <- ab[2] ## Vaccinated, influenza positive
b <- ab[1] ## Unvaccinated, influenza positive
c <- cd[2] ## Vaccinated, influenza negative
d <- cd[1] ## Unvaccinated, influenza negative

print(1 - (a/b)/(c/d))
ggplot(inf_states_non_flu) + geom_histogram(aes(x=non_flu_time,fill=as.factor(vacc_status_at_inf)),binwidth=1)
ggplot(inf_states_flu) + geom_histogram(aes(x=flu_time,fill=as.factor(vacc_status_at_inf)),binwidth=1)
