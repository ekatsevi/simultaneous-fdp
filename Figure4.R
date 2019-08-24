###############################################################
#
# Simulate FDP bounds in the sorted setting 
# (to reproduce Figure 4).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 4...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n = 2500                                # number of hypotheses
alpha = 0.1                             # confidence level
num_repeats = 100                       # number of Monte Carlo repetitions
a = 1                                   # additive FDP regularization
non_nulls_vals = c(100, 200, 300)       # number of non-nulls
mu_vals = c(2,3,4)                      # signal strength
bounds = c("proposed", "GS_Simes",      # FDP bounds considered
           "GS_Fisher", "DKW", "truth") 

# set up array to store results
num_non_nulls_vals = length(non_nulls_vals)
num_mu_vals = length(mu_vals)
num_bounds = length(bounds)
  
FDP_bar_reps = array(NA, 
                     c(num_non_nulls_vals,num_mu_vals, num_bounds, n, num_repeats), 
                     list(non_nulls_vals, mu_vals, bounds, 1:n, 1:num_repeats))

# compute multipler c based on Theorem 1
c = log(1/alpha)/log(1 + log(1/alpha))

# run simulation for all values of num_non_nulls and mu
for(non_null_idx in 1:num_non_nulls_vals){
  non_nulls = non_nulls_vals[non_null_idx]
  cat(sprintf("Working on num_non_nulls = %d...\n", non_nulls))
  non_nulls_bool = rep(FALSE, n)
  non_nulls_bool[1:non_nulls] = TRUE
  for(mu_idx in 1:num_mu_vals){
    mu = mu_vals[mu_idx]
    cat(sprintf("Working on mu = %0.2f...\n", mu))
    for(rep in 1:num_repeats){
      # generate data based on normal means model
      Z = rnorm(n)
      Z[non_nulls_bool] = Z[non_nulls_bool] + mu
      
      P = pnorm(Z, lower.tail = FALSE)
      P_sorted = sort(P)
      P_order = order(P)
      hommel = hommelFast(P)
      
      FDP_bar_reps[non_null_idx, mu_idx, "proposed", , rep] = pmin(1, c*(a + n*P_sorted)/(1:n))
      FDP_bar_reps[non_null_idx, mu_idx, "GS_Simes", , rep] = 1 - curveSimes(hommel, order = P_order, alpha = alpha, plot = FALSE)/(1:n)
      FDP_bar_reps[non_null_idx, mu_idx, "GS_Fisher", , rep] = 1 - curveFisher(P, order = P_order, alpha = alpha, plot = FALSE)/(1:n)
      FDP_bar_reps[non_null_idx, mu_idx, "DKW", , rep] = pmin(1, (sqrt(n/2*log(1/alpha)) + n*P_sorted)/(1:n))
      FDP_bar_reps[non_null_idx, mu_idx, "truth", , rep] = FDP_path
    }
  }
}
# compute mean FDP bounds and upper quantile of FDP
FDP_bar = apply(FDP_bar_reps[,,c("proposed", "GS_Simes", "GS_Fisher", "DKW"),,], c(1,2,3,4), mean)
FDP = apply(FDP_bar_reps[,,"truth",,], c(1,2,3), function(vec)(quantile(vec,1-alpha)))

# massage results for plotting
FDP_df = as_tibble(melt(FDP))
names(FDP_df) = c("non_nulls", "mu", "index", "FDP")
FDP_df$bound = "truth"

FDP_bar_df = as_tibble(melt(FDP_bar))
names(FDP_bar_df) = c("non_nulls", "mu", "bound"  ,"index", "FDP")

FDP_bar_df = FDP_bar_df %>% select(non_nulls, mu, index, FDP, bound) %>% rbind(FDP_df)
  
FDP_bar_df$non_nulls = factor(FDP_bar_df$non_nulls, levels = non_nulls_vals, 
                              labels = sapply(non_nulls_vals, function(non_nulls)(sprintf("%d non-nulls", non_nulls)))) 
FDP_bar_df$mu = factor(FDP_bar_df$mu, levels = mu_vals, labels = c("Weak signal", "Medium signal", "Strong signal"))
FDP_bar_df$bound = factor(FDP_bar_df$bound, 
                          levels = c("proposed", "GS_Simes",  "GS_Fisher", "DKW", "truth"), 
                          labels = c("Proposed", "GS (Simes)", "GS (Fisher)", "DKW",  "(FDP)"))

# plot results
p = FDP_bar_df %>%
  filter(index < 300) %>%
  ggplot(aes(x = index, y = FDP, group = bound)) + 
  geom_line(aes(colour = bound, linetype = bound)) + 
  facet_grid(non_nulls ~ mu) + 
  scale_colour_manual(values = c("magenta", "blue",  "cyan", "orangered", "black")) +
  scale_linetype_manual(values = c("solid", "longdash", "longdash", "dotdash", "dotted")) +
  theme_custom() + theme(legend.position = "bottom", legend.title = element_blank()) + 
  xlab("Hypothesis Index") + ylab("FDP bound")
plot(p)
ggsave(filename = sprintf("%s/Figure4.pdf", figures_dir), plot = p, device = cairo_pdf, width = textwidth, height = textwidth)