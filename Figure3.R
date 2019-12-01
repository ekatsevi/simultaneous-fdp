###############################################################
#
# Compare proposed bounds to DKW, Robbins bounds
# (to reproduce Figure 3).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 3...\n"))

# set up workspace
source("setup.R")

# simulation parameters
n = 2500           # number of hypotheses
alpha = 0.1        # confidence level
t = seq(1e-5,1,    # p-value cutoffs
        by = 1e-5)

# compute constant c from Theorem 1
c = log(1/alpha)/log(1 + log(1/alpha))

# compute each bound
Robbins = floor(n*t/alpha)
Proposed = floor(c*(1+n*t))
DKW = floor(sqrt(n/2*log(1/alpha)) + n*t)
Quantile = qbinom(1-alpha, n, t)
df = data.frame(t, Robbins, Proposed, DKW, Quantile)
names(df)[5] = "(Pointwise FDP Quantile)"

# plot bounds
df_melt = melt(df, 1)
p1 = ggplot(data = df_melt, aes(x = t, y = value, group = variable)) + scale_x_log10() + scale_y_log10() + 
  geom_line(aes(color = variable, linetype = variable)) + theme_custom() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  geom_vline(xintercept = alpha/n, linetype = "dotted") + geom_vline(xintercept = alpha, linetype = "dotted") + 
  scale_colour_manual(name = "Bound", values = c("blue", "magenta", "orangered", "black")) + 
  scale_linetype_manual(name = "Bound", values = c("longdash", "solid", "dotdash", "dotted"))  +
  ylab("Bound on V(t)")
plot(p1)

# compute where bounds are tightest
reps = 10000
Robbins = numeric(reps)
Proposed = numeric(reps)
DKW = numeric(reps)
for(rep in 1:reps){
  p = sort(runif(n))
  Robbins[rep] = p[which.max((1:n)/(n*p/alpha))]
  Proposed[rep] = p[which.max((1:n)/(c*(1+n*p)))]
  DKW[rep] = p[which.max((1:n)/(sqrt(n/2*log(1/alpha)) + n*p))]
}

frac_in_range = mean(alpha/n <= Proposed & Proposed <= alpha)
cat(sprintf("Tightest bound achieved in interesting region %0.2f of the time.", frac_in_range))

# plot results
df = data.frame(Robbins, Proposed, DKW)
df_melt = melt(df)

p2 = ggplot(data = df_melt, aes(x = value, group = variable)) + 
  geom_histogram(aes(fill = variable, y = ..density..), position = "identity", alpha = 0.5, bins = 50) + 
  scale_x_log10() + 
  geom_vline(xintercept = alpha/n, linetype = "dotted") + geom_vline(xintercept = alpha, linetype = "dotted") + 
  scale_fill_manual(name = "Bound", values = c("blue", "magenta", "orangered")) + ylab("Frequency") +
  coord_cartesian(ylim=c(0, 1), xlim = c(1e-6,1)) + xlab("value of t where bound is tightest") + 
  theme_custom() + theme(legend.position = "bottom") 
plot(p2)

# combine plots
g <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
p1 = p1 + theme(legend.position="none")
p2 = p2 + theme(legend.position="none")
p1b <- ggplot_build(p1)
p2b <- ggplot_build(p2)
g1 = ggplot_gtable(p1b)
g2 = ggplot_gtable(p2b)
g <- cbind(g1, g2, size = "first")
grid.newpage()
grid.draw(g)
lheight <- sum(legend$height)
p = grid.arrange(g, legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
plot(p)
ggsave(filename = sprintf("%s/Figure3.pdf", figures_dir), plot = p, device = cairo_pdf, width = textwidth, height = 0.6*textwidth)