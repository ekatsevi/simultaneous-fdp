###############################################################
#
# Simulate FDP bounds in the ordered setting 
# (to reproduce Figure 5).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 5...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n = 2500                           # number of hypotheses
alpha = 0.1                        # confidence level 
num_repeats = 100                  # number of Monte Carlo repetitions
a = 1                              # additive FDP regularization
num_non_nulls = 100                # number of non-null hypotheses
lambda = 0.1                       # parameter for accumulation function
theta_vals = c(15, 35, 55)         # control locations of non-nulls
mu_vals = c(2,3,4)                 # signal strength
bounds = c("proposed", "GS_Simes", # FDP bounds considered
           "GS_Fisher", "truth")

# set up structures to store the results
num_theta_vals = length(theta_vals)
num_mu_vals = length(mu_vals)
num_bounds = length(bounds)

FDP_bar_reps = array(NA, 
                     c(num_theta_vals, num_mu_vals, num_bounds, n, num_repeats), 
                     list(theta_vals, mu_vals, bounds, 1:n, 1:num_repeats))

non_nulls_vals = matrix(NA, num_non_nulls, num_theta_vals)
colnames(non_nulls_vals) = theta_vals

# compute multiplier c based on Theorem 2
B = 1/(1-lambda)
c = -log(1/alpha)/log(1 - (1-alpha^B)/B)

# run simulation for all values of theta and mu
for(theta_idx in 1:num_theta_vals){
  theta = theta_vals[theta_idx]
  cat(sprintf("Working on theta = %0.2f...\n", theta))
  # define non-nulls based on theta
  p = dexp(1:n, rate = theta/n)/sum(dexp(1:n, rate = theta/n))
  non_nulls = sort(sample(x = 1:n, size = num_non_nulls, replace = FALSE, prob = p))
  non_nulls_vals[,as.character(theta)] = non_nulls
  non_nulls_bool = rep(FALSE, n)
  non_nulls_bool[non_nulls] = TRUE
  for(mu_idx in 1:num_mu_vals){
    mu = mu_vals[mu_idx]
    cat(sprintf("Working on mu = %0.2f...\n", mu))
    for(rep in 1:num_repeats){
      # generate data based on normal means model
      Z = rnorm(n)
      Z[non_nulls] = Z[non_nulls] + mu
      P = pnorm(Z, lower.tail = FALSE)
      hommel = hommelFast(P)
      
      # compute all FDP bounds
      FDP_bar_reps[theta_idx, mu_idx, "proposed", , rep] = pmin(1, c*(a+1/(1-lambda)*cumsum(P > lambda))/(1:n))
      FDP_bar_reps[theta_idx, mu_idx, "GS_Simes", , rep] = 1 - curveSimes(hommel, order = 1:n, alpha = alpha, plot = FALSE)/(1:n)
      FDP_bar_reps[theta_idx, mu_idx, "GS_Fisher", , rep] = 1 - curveFisher(P, order = 1:n, alpha = alpha, plot = FALSE)/(1:n)
      FDP_bar_reps[theta_idx, mu_idx, "truth", , rep] = cumsum(!non_nulls_bool)/(1:n)
    }
  }
}

# compute mean FDP bounds and upper quantile of FDP
FDP_bar = apply(FDP_bar_reps[,,c("proposed", "GS_Simes", "GS_Fisher"),,], c(1,2,3,4), mean)
FDP = apply(FDP_bar_reps[,,"truth",,], c(1,2,3), function(vec)(quantile(vec,1-alpha)))

# massage results for plotting
mu_vals_df = tibble(mu_vals)
names(mu_vals_df) = "mu"
non_nulls_df = as_tibble(non_nulls_vals) %>% gather("theta", "index") %>% merge(mu_vals_df) %>% as_tibble()

FDP_df = as_tibble(melt(FDP))
names(FDP_df) = c("theta", "mu", "index", "FDP")
FDP_df$bound = "truth"

FDP_bar_df = as_tibble(melt(FDP_bar))
names(FDP_bar_df) = c("theta", "mu", "bound", "index", "FDP")

FDP_bar_df = FDP_bar_df %>% select(theta, mu, index, FDP, bound) %>% rbind(FDP_df)

FDP_bar_df$theta = factor(FDP_bar_df$theta, levels = theta_vals, 
                              labels = c("Weak ordering", "Medium ordering", "Strong ordering")) 
FDP_bar_df$mu = factor(FDP_bar_df$mu, levels = mu_vals, labels = c("Weak signal", "Medium signal", "Strong signal"))
FDP_bar_df$bound = factor(FDP_bar_df$bound, 
                          levels = c("proposed", "GS_Simes",  "GS_Fisher", "truth"), 
                          labels = c("Proposed", "GS (Simes)", "GS (Fisher)", "(True FDP)"))

non_nulls_df$theta = factor(non_nulls_df$theta, levels = theta_vals, labels = c("Weak ordering", "Medium ordering", "Strong ordering"))
non_nulls_df$mu = factor(non_nulls_df$mu, levels = mu_vals, labels = c("Weak signal", "Medium signal", "Strong signal"))

# plot results
p = FDP_bar_df %>%
  filter(index < 250) %>%
  ggplot(aes(x = index, y = FDP, group = bound)) + 
  geom_line(aes(colour = bound, linetype = bound)) + 
  geom_rug(data = non_nulls_df %>% filter(index < 250), aes(x = index), sides = "b", inherit.aes = FALSE) +
  facet_grid(theta ~ mu) + 
  scale_colour_manual(values = c("magenta", "blue",  "cyan", "black")) +
  scale_linetype_manual(values = c("solid", "longdash", "longdash", "dotted")) +
  theme_custom() + theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab("Hypothesis Index") + ylab("FDP bound")
plot(p) 
ggsave(filename = sprintf("%s/Figure5.pdf", figures_dir), plot = p, device = cairo_pdf, width = textwidth, height = textwidth)
