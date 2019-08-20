# clear working directory
rm(list = ls())

# set paths to figures and stored data
data_dir = "data"
figures_dir = "figures"
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}

# load R packages
library(tidyverse)
library(cherry)
library(reshape2)
library(downloader)
require(gtable)
require(grid)
library(gridExtra)
library(cowplot)
library(latex2exp)
library(kableExtra)

# source helper functions
source("UKBB_utils.R")

# set seed for reproducibility
set.seed(123)

# set some plotting parameters
textwidth = 4.9823
font.size = 10
theme_custom <- function(){ 
  theme_bw() + theme(text = element_text(size = font.size), 
                     legend.text = element_text(size = font.size),
                     strip.text = element_text(size = font.size))
}