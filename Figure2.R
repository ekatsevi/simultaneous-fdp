###############################################################
#
# Plot FDP-bar and FDP-hat for UK Biobank platelet count data,
# as well as enrichment of "Blood coagulation" and "Platelet
# activation" GO terms using GREAT. 
# (to reproduce Figure 2).
# 
# Author: Eugene Katsevich
# Date:   8/16/2019
###############################################################

cat(sprintf("Working on Figure 2...\n"))

# set up workspace
source("setup.R")

# analysis parameters
alpha = 0.05              # confidence level for simultaneous FDP bound
q = 0.1                   # FDR target level
phenotype = "platelet"    # phenotype of interest

# download KnockoffZoom data, if necessary
download_KZ_data(data_dir)

# read in processed data
df = read_KZ_data(data_dir, phenotype)

# compute knockoffs FDP multipler
C = -log(alpha)/log(2-alpha)

# compute FDP-hat and FDP-bar
FDP_hat_vec = get_FDP_hat(df$W)
FDP_bar_vec = get_FDP_bar(df$W, C)
df$FDP_hat = FDP_hat_vec
df$FDP_bar = FDP_bar_vec
df$num_rejections = 1:nrow(df)

# find set of interesting points along the path
values = df %>% filter(num_rejections <= 1000) %>% pull(FDP_bar)
k_vals = get_informative_points(values, min_val_improvement = 0.01, min_k_improvement = 50)
informative_points = data.frame(k_vals, values[k_vals])
names(informative_points) = c("num_rejections", "FDP_bar")

# write regions corresponding to each interesting rejection set to file
if(!dir.exists(sprintf("%s/bed_files", data_dir))){
  dir.create(sprintf("%s/bed_files", data_dir))
}
for(k in k_vals){
  filename = sprintf("%s/bed_files/regions_platelet_%d.bed", data_dir, k)
  df %>% 
    filter(Phenotype == "platelet", num_rejections <= k, !is.na(BP.min)) %>% 
    select(CHR, BP.min, BP.max) %>%
    mutate(CHR = sprintf("chr%d", CHR)) %>%
    write_tsv(path = filename, col_names = FALSE)
}

# [compute GREAT enrichment for each interesting rejection set via 
# http://great.stanford.edu/public/html/, using the input files in /data_dir/bed_files
# and save results to /data_dir/GREAT_output]

# read GREAT results
great_results = lapply(k_vals, function(k)(read_great_results(data_dir, k)))
great_results = do.call("rbind", great_results)

# pull out fold enrichment for "Blood coagulation" and "Platelet activation" terms
informative_points$"coagulation" = unname(unlist(great_results %>% filter(Desc == "blood coagulation") %>% select(RegionFoldEnrich)))
informative_points$"platelet" = unname(unlist(great_results %>% filter(Desc == "platelet activation") %>% select(RegionFoldEnrich)))

# transform fold enrichment to FDP scale for plotting
M_enrichment = max(c(informative_points$coagulation, informative_points$platelet))
m_enrichment = min(c(informative_points$coagulation, informative_points$platelet))
M_FDP = 0.25
m_FDP = 0
mult_factor = (M_enrichment - m_enrichment)/(M_FDP - m_FDP)
add_factor = M_enrichment - M_FDP*mult_factor

df_to_plot_1 = informative_points %>% 
  select(num_rejections, FDP_bar, coagulation, platelet) %>%
  mutate_at(c("coagulation", "platelet"), function(enrichment)((enrichment-add_factor)/mult_factor)) %>%
  gather(key, value, coagulation, platelet, FDP_bar) %>%
  filter(key %in% c("coagulation", "platelet"))

# pull out first 1500 rejection sets on the path
df_to_plot_2 = df %>% 
  filter(num_rejections <= 1500) %>%
  select(num_rejections, FDP_hat, FDP_bar) %>%
  gather(key, value, FDP_hat, FDP_bar)

# combine FDP information and enrichment information
df_to_plot = rbind(df_to_plot_1, df_to_plot_2) %>%
  mutate(key = factor(key, 
                      levels = c("FDP_bar", "FDP_hat", "coagulation", "platelet"),
                      labels = c("FDP bound", "FDP estimate", 
                                 "Blood coagulation", "Platelet activation")))

df_to_plot_points = df_to_plot %>% 
  filter(num_rejections %in% informative_points$num_rejections, 
         key %in% c("FDP bound", "Blood coagulation", "Platelet activation"))

# plot results
p = df_to_plot %>% 
  ggplot(aes(x = num_rejections, y = value, group = key, colour = key)) + 
  geom_line(aes(colour = key, linetype = key)) + 
  geom_point(data = df_to_plot_points, aes(x = num_rejections, y = value, group = key, fill = key), shape = 21, colour = "transparent", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_colour_manual(values = c("magenta", "darkblue", "goldenrod3", "darkorange2"), 
                      labels = c("FDP bound", "FDP estimate", "Blood coagulation", "Platelet activation")) + 
  scale_linetype_manual(values = c("solid", "dotdash", "longdash", "longdash")) + 
  scale_fill_manual(values = c("magenta", "goldenrod3", "darkorange2")) + 
  scale_y_continuous(limits = c(0, 0.25), sec.axis = sec_axis(~.*mult_factor + add_factor, name = "Fold enrichment")) + 
  guides(fill=FALSE, colour = guide_legend(nrow=1, byrow = TRUE)) +
  theme_custom() + 
  theme(legend.title = element_blank(), legend.position = "bottom", legend.spacing.x = unit(0.05, 'cm')) + 
  xlab("Number of rejections") + ylab("FDP")
plot(p)
ggsave(filename = sprintf("%s/Figure2.pdf", figures_dir), plot = p, device = cairo_pdf, width = textwidth, height = 0.75*textwidth)