###############################################################
#
# Examine the ability of FDP estimates to bound the true FDP
# (to reproduce Figure 5).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 5...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n_vals = c(1000, 2500, 5000)         # number of hypotheses
alpha = 0.1                          # confidence level
num_repeats = 1000                   # number of Monte Carlo repetitions
delta_vals = 10^(seq(-5, -1,         # minimum p-value considered
                 length.out = 1000))

# set up array to store results
num_delta_vals = length(delta_vals)
num_n_vals = length(n_vals)
results = array(0, c(num_delta_vals, num_n_vals, num_repeats))

# run simulation
for(rep in 1:num_repeats){
  # generate independent uniform p-values
  P = runif(max(n_vals))
  for(n_idx in 1:num_n_vals){
    # compute degree to which FDP exceeds FDP-hat over the interval [delta, 1]
    n = n_vals[n_idx]
    P_sorted = sort(P[1:n])
    overshoot = 1:n/(n*P_sorted)
    results[,n_idx,rep] = sapply(delta_vals, function(delta)(max(overshoot[P_sorted > delta])))
  }
}

# find 1-alpha quantile of overshoot
overshoots = apply(results, c(1,2), function(vec)(quantile(vec, 1-alpha)))

# massage results for plotting
colnames(overshoots) = n_vals
overshoots = as_tibble(overshoots)
overshoots$delta = delta_vals

df_to_plot = overshoots %>% gather("n", "overshoot", -delta)

closest_delta_vals = sapply(alpha/n_vals, function(val)(delta_vals[which.min(abs(delta_vals - val))]))
overshoot_vals = sapply(1:num_n_vals, function(n_idx)(df_to_plot %>% filter(near(delta, closest_delta_vals[n_idx]), n == n_vals[n_idx]) %>% pull(overshoot)))
points_to_plot = tibble(closest_delta_vals, overshoot_vals, n_vals)
names(points_to_plot) = c("delta", "overshoot", "n")

df_to_plot$n = as.factor(df_to_plot$n)
points_to_plot$n = as.factor(points_to_plot$n)

# plot results
p = df_to_plot %>% ggplot(aes(x = delta, y = overshoot, group = n)) + geom_line(aes(colour = n, linetype = n)) + 
  geom_point(data = points_to_plot, aes(x = delta, y = overshoot, group = n, colour = n), inherit.aes = FALSE) + 
  scale_x_continuous(trans = "log", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) + 
  scale_y_continuous(trans = "log", breaks = c(1,2,4,8)) +  
  scale_colour_manual(values = c("red", "purple", "blue")) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash")) + 
  xlab(TeX("Minimum allowable p-value $\\delta$")) + ylab(TeX('$\\sup_{t \\geq \\delta} FDP(t)/\\widehat{FDP}(t)$')) + 
  theme_custom() + theme(legend.position = "right")
plot(p)
ggsave(filename = sprintf("%s/Figure5.pdf", figures_dir), plot = p, device = cairo_pdf, width = 0.9*textwidth, height = 0.6*textwidth)