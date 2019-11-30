###############################################################
#
# Examine the ability of FDP estimates to bound the true FDP
# (to reproduce Figure 5).
# 
# Author: Eugene Katsevich
# Date:   11/30/2019
###############################################################

cat(sprintf("Working on Figure 5...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n = 2500                    # number of hypotheses
n_1 = 100                   # number of non-nulls
mu = 4                      # signal strength
reps = 1000                 # number of Monte Carlo repetitions
q_subset = c(0.01,.025,.05, # special set of values of q
             .075,.1,.125,
             .15,.175,.2)

# minimum allowable FDR levels
min_q_vals = 10^(seq(-3, -1,        
                     length.out = 1000))

# non-null configuration
nonnulls = logical(n)
nonnulls[1:n_1] = TRUE

# array to store results
max_discrepancies = matrix(0, reps, length(min_q_vals), dimnames = list(NULL, min_q_vals))
max_discrepancies_special = numeric(reps)

# repeatedly generate p-values from the independent normal means model, 
# and calculate the maximum ratio between the realized and nominal FDP values.
for(rep in 1:reps){
  Z = rnorm(n)
  Z[nonnulls] = Z[nonnulls] + mu
  P = pnorm(Z, lower.tail = FALSE)
  P.adjusted = p.adjust(P, "fdr")
  FDP_ratio = function(q)(sum(!nonnulls[P.adjusted <= q])/max(1, sum(P.adjusted <= q))/q)
  discrepancies = sapply(min_q_vals, FDP_ratio)
  max_discrepancies[rep,] = rev(cummax(rev(discrepancies)))
  max_discrepancies_special[rep] = max(sapply(q_subset, FDP_ratio))
}

# massage results for plotting
df = melt(max_discrepancies)
names(df) = c("rep", "min_q_val", "discrepancy")
df = as_tibble(df)
special_df = tibble(min_q_val = min(q_subset),
                    "Mean" = mean(max_discrepancies_special), 
                    "90% quantile" = quantile(max_discrepancies_special, 0.9)) %>% 
  gather(key = Measure, value = discrepancy, -min_q_val)

# plot results
p = df %>% 
  group_by(min_q_val) %>% 
  summarise("Mean" = mean(discrepancy), 
            "90% quantile" = quantile(discrepancy, 0.9)) %>% 
  gather(key = Measure, value = discrepancy, -min_q_val) %>% 
  ggplot(aes(x = min_q_val, y = discrepancy, group = Measure, colour = Measure, linetype = Measure)) + 
  geom_line() + theme_bw() + scale_x_log10(breaks = c(0.001, 0.01, 0.1), labels = c("0.001", "0.01", "0.1")) + 
  scale_y_log10(breaks = c(1,2,4,8)) + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  geom_point(data = special_df, show.legend = FALSE) + xlab("Minimum nominal FDR level") + 
  ylab("Factor FDP exceeds nominal FDR") +  theme_custom()
plot(p)
ggsave(plot = p, filename = sprintf("%s/Figure5.pdf", figures_dir),
       device = cairo_pdf, width = 0.9*textwidth, height = 0.6*textwidth)
