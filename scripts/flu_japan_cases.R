library(tidyverse)

file_input <- "~/Desktop/flu_japan_cases.csv"

df <- read_csv(file_input)

ggplot(df) +
  geom_point(aes(Date, cases)) +
  theme_classic() +
  scale_y_log10()
#ggsave("~/Desktop/flu_japan_cases.pdf", height = 6, width = 6)

# r = exponential growth rate
coefs <- df %>%
  mutate(days = difftime(Date, min(Date), units = "days") %>% as.integer()) %>%
  filter(days > 8) %>%
  {lm(data = ., log(cases) ~ days)} %>%
  coef
r <- coefs[["days"]]

# Generation time distribution from https://pmc.ncbi.nlm.nih.gov/articles/PMC11370535/
# They report only mean and sd; assume it's gamma...
mean <- 3.2
sd <- 2.1
alpha <- mean^2 / sd^2
beta <- alpha / mean

# ...check a gamma looks like their plot
ggplot(tibble(t = seq(0, 10, by = 0.01),
              y = dgamma(t, shape = alpha, rate = beta))) +
  geom_line(aes(t, y))

R <- 1 / integrate(function(t) {dgamma(t, shape = alpha, rate = beta) * exp(-r * t)},
          0, Inf)$value 
R
