###############################################################
#
# Simulate FDP bounds in under p-value correlations 
# (to reproduce Figure 7).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 7...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n = 2500                              # number of hypotheses
alpha = 0.1                           # confidence level
num_repeats = 1000                    # number of Monte Carlo repetitions
a = 1                                 # additive FDP regularization
p_star_online = 0.05                  # p-value threshold for online setting
lambda_preorder_acc = 0.1             # accumulation function parameter for preordered setting
p_star_preorder_sel = 0.1             # p-value threshold for preordered setting
rho_vals = seq(-0.9, 0.9, by = 0.1)   # correlation parameter
methods = c("sort", "preorder_acc",   # FDP bounds considered
            "preorder_sel", "online")

# set up array to store the results
num_rho_vals = length(rho_vals)
num_methods = length(methods)
FDP_bar_reps = array(NA, 
                     c(num_rho_vals, num_methods, num_repeats), 
                     list(rho_vals, methods, 1:num_repeats))

# compute constants c for each bound
c_vals = numeric(length(methods))
names(c_vals) = methods

# constant for sorted case, based on Theorem 1
c_vals["sort"] = log(1/alpha)/log(1 + log(1/alpha))

# constant for pre-ordered case (acc), based on Theorem 2
B = 1/(1-lambda_preorder_acc)
c_vals["preorder_acc"] = -log(1/alpha)/log(1 - (1-alpha^B)/B)

# constant for pre-ordered case (sel), based on Theorem 2
B = p_star_preorder_sel/(1-p_star_preorder_sel)
c_vals["preorder_sel"] = log(1/alpha)/log(1 + (1-alpha^B)/B)

# constant for online case, based on Theorem 4
c_vals["online"] = log(1/alpha)/log(1 + log(1/alpha))

# run simulation for all values of rho
for(rho_idx in 1:num_rho_vals){
  rho = rho_vals[rho_idx]
  cat(sprintf("Working on rho = %0.2f...\n", rho))
  # compute correlation matrix and its Cholesky factor
  Sigma = rho^(abs(outer(1:n, 1:n, "-")))
  Sigma_half = t(chol(Sigma))
  for(rep in 1:num_repeats){
    # generate correlated data from global null
    Z = Sigma_half%*%rnorm(n)
    P = pnorm(Z, lower.tail = FALSE)
    P_sorted = sort(P)
    P_order = order(P)
    
    # compute all FDP bounds
    FDP_bar_reps[rho_idx, "sort", rep] = min(pmin(1, floor(c_vals["sort"]*(a + n*P_sorted))/(1:n)))
    FDP_bar_reps[rho_idx, "preorder_acc", rep] = min(pmin(1, floor(c_vals["preorder_acc"]*(a+1/(1-lambda_preorder_acc)*cumsum(P > lambda_preorder_acc)))/(1:n)))
    FDP_bar_reps[rho_idx, "preorder_sel", rep] = min(pmin(1, floor(c_vals["preorder_sel"]*(a+p_star_preorder_sel/(1-p_star_preorder_sel)*cumsum(P > p_star_preorder_sel)))/pmax(1,cumsum(P <= p_star_preorder_sel))))
    FDP_bar_reps[rho_idx, "online", rep] = min(pmin(1,floor(c_vals["online"]*(a + (1:n)*p_star_online))/pmax(1,cumsum(P < p_star_online))))
  }
}
# compute lower alpha quantile of FDP_bar
FDP_bar = as_tibble(apply(FDP_bar_reps, c(1,2), function(vec)(quantile(vec,alpha))))

# massage results for plotting
FDP_bar$rho = rho_vals
df_to_plot = FDP_bar %>% 
  gather(method, FDP_bar, -rho) %>% 
  mutate(method = factor(method, 
                         levels = methods, 
                         labels = c("sort", "preorder-acc", "preorder-sel", "online")))

# plot results
p = df_to_plot %>% 
  ggplot(aes(x = rho, y = 1/FDP_bar, group = method)) + 
  geom_line(aes(colour = method, linetype = method)) + 
  xlab(TeX("Correlation $\\rho$")) + ylab(TeX('max $FDP/\\bar{FDP}$')) + 
  scale_y_continuous(limits = c(1,2), breaks = c(1,1.25,1.5, 1.75, 2)) +
  scale_linetype_manual(values = c("solid", "longdash", "longdash", "dotdash")) + 
  theme_custom() + theme(legend.position = "right", legend.title = element_blank())
plot(p)
ggsave(filename = sprintf("%s/Figure7.pdf", figures_dir), plot = p, device = cairo_pdf, width = textwidth, height = 0.6*textwidth)