# Fit replicator frequencies with one reference fitness fixed to 0
# Data must have columns: time (numeric), clade (factor/character), freq (observed frequency in [0,1])
# Example: data <- read.csv("subclade_freqs.csv")

library(dplyr)
library(ggplot2)

data <- h3_dat
colnames(data) <- c("clade","time","n","n_tot","freq")
# ---------- USER SETTINGS ----------
ref_clade <- "J.2.3" #"K"  # set to a clade name (e.g. "A") or NULL to auto-pick the first observed clade
# -----------------------------------

# ----- sanity check & setup -----
data <- data %>% arrange(time, clade) %>% ungroup()
#data <- tidyr::complete(data, time, clade, fill = list(n = 0, n_tot = 1, freq = 0))%>% arrange(time, clade) %>% ungroup()
clades <- unique(as.character(data$clade))
if (is.null(ref_clade)) ref_clade <- clades[1]
if (!(ref_clade %in% clades)) stop("Reference clade not found in data")

# initial frequencies at earliest time (assume these exist for every clade)
t0 <- min(data$time)
p_init_df <- data %>% filter(time == t0) %>% select(clade, freq)
if (nrow(p_init_df) != length(clades)) {
  warning("Initial frequencies at t0 do not include every clade. Missing clades will get a tiny pseudocount.")
  # build p_init with small pseudocounts for missing
  p_init <- setNames(rep(1e-6, length(clades)), clades)
  for (r in seq_len(nrow(p_init_df))) p_init[p_init_df$clade[r]] <- p_init_df$freq[r]
  p_init <- p_init / sum(p_init)
} else {
  p_init <- setNames(p_init_df$freq, p_init_df$clade)
}

# reorder clades so ref is first (convenience)
clades <- c(ref_clade, setdiff(clades, ref_clade))
K <- length(clades)

# ----- helper that maps reduced params -> full f vector -----
# par_reduced length = K-1, f_ref = 0
make_full_f <- function(par_red, clades) {
  f_full <- numeric(length(clades))
  names(f_full) <- clades
  f_full[1] <- 0
  if (length(par_red) > 0) f_full[-1] <- par_red
  f_full
}

# ----- replicator prediction given full f -----
predict_freqs <- function(f_full, times, p_init) {
  # returns matrix (rows = times, cols = clades) of predicted frequencies
  times_unique <- sort(unique(times))
  preds <- sapply(times_unique, function(t) {
    num <- p_init * exp(f_full * t)
    num / sum(num)
  })
  # transpose to data frame in long form
  preds_df <- data.frame(time = rep(times_unique, each = length(f_full)),
                         clade = rep(names(f_full), times = length(times_unique)),
                         pred = as.vector(preds))
  preds_df
}

# ----- objective: squared error in logit space -----
objective <- function(par_red, data, p_init, clades) {
  f_full <- make_full_f(par_red, clades)
  preds_df <- predict_freqs(f_full, data$time, p_init)
  # align order
  preds_vec <- preds_df$pred[match(paste(data$time, data$clade),
                                   paste(preds_df$time, preds_df$clade))]
  obs <- pmax(pmin(data$freq, 1-1e-8), 1e-8) # avoid 0/1
  sum((qlogis(obs) - qlogis(preds_vec))^2)
}

# ----- initial guess: zeros for non-ref clades -----
init_par <- rep(0, K-1)

# ----- run optimization (BFGS to try get Hessian) -----
opt <- optim(par = init_par,
             fn  = objective,
             data = data,
             p_init = p_init[clades],
             clades = clades,
             method = "BFGS",
             control = list(maxit = 5000),
             hessian = TRUE)

# ----- collect results -----
f_hat_full <- make_full_f(opt$par, clades)
names(f_hat_full) <- clades

# approximate SEs from Hessian (if invertible)
se <- rep(NA, length(f_hat_full)); names(se) <- clades
if (!is.null(opt$hessian) && nrow(opt$hessian) == (K-1)) {
  h <- opt$hessian
  inv_h <- tryCatch(solve(h), error = function(e) NULL)
  if (!is.null(inv_h)) {
    # variance approx for reduced params; set NA for ref
    var_red <- diag(inv_h)
    se[-1] <- sqrt(pmax(0, var_red))
  }
}

# ----- predictions for plotting -----
data_complete <- tidyr::complete(data, time, clade, fill = list(n = 0, n_tot = 1, freq = 0))%>% arrange(time, clade) %>% ungroup()
pred_df <- predict_freqs(f_hat_full, data_complete$time, p_init[clades])
data_complete$pred <- pred_df$pred[match(paste(data_complete$time, data_complete$clade), paste(pred_df$time, pred_df$clade))]

# ----- output -----
cat("Reference clade (fixed):", ref_clade, " (f = 0)\n\n")
generation_interval <- 3/7
est_table <- data.frame(clade = clades, f_hat = f_hat_full, se = se,growth_advantage=exp(f_hat_full*generation_interval))
print(est_table)

# ----- plot observed vs predicted -----
ggplot(data_complete, aes(x = time, color = clade)) +
  geom_point(aes(y = freq), alpha = 0.6) +
  geom_line(aes(y = pred), size = 1) +
  labs(y = "frequency", title = "Replicator-fit: observed (points) vs predicted (lines)") +
  theme_minimal()


