# download KnockoffZoom summary statistics
download_KZ_data = function(data_dir){
  discoveries_dir=sprintf("%s/ukb_discoveries", data_dir)
  statistics_dir=sprintf("%s/ukb_stats", data_dir)
  discoveries_file="ukb_discoveries.zip"
  statistics_file="ukb_stats.zip"
  remote="https://msesia.github.io/knockoffzoom/data"

  if(!dir.exists(discoveries_dir)){
    cat(sprintf("Downloading KnockoffZoom discoveries...\n"))
    url = sprintf("%s/%s", remote, discoveries_file)
    dest = sprintf("%s/%s", data_dir, discoveries_file)
    download(url, dest, mode="wb") 
    unzip(dest, exdir = data_dir)    
    file.remove(dest)
  }

  if(!dir.exists(statistics_dir)){
    cat(sprintf("Downloading KnockoffZoom statistics...\n"))
    url = sprintf("%s/%s", remote, statistics_file)
    dest = sprintf("%s/%s", data_dir, statistics_file)
    download(url, dest, mode="wb") 
    unzip(dest, exdir = data_dir)    
    file.remove(dest)
  }
}

# helper function to read KZ data for one phenotype
read_phenotype = function(phenotype, data_dir){
  statistics = read_delim(file = sprintf("%s/ukb_stats/ukb_%s_kz_res2.txt", data_dir, phenotype), delim = " ", col_types = "iiciiddc")
  discoveries = read_delim(file = sprintf("%s/ukb_discoveries/ukb_%s_kz_res2.txt", data_dir, phenotype), delim = " ", col_types = "icidiiiii")
  R2.max = 0.99     # maximum allowable correlation between knockoff and original variable
  statistics = statistics %>% filter(is.na(R2)|R2<=R2.max) %>% select(CHR, Group, W)
  discoveries = discoveries %>% select(CHR, Group, BP.min, BP.max)
  df = left_join(statistics, discoveries, by = c("CHR", "Group"))
  # sort by decreasing order of magnitude
  mag_order = order(abs(df$W), decreasing = TRUE)
  df = df[mag_order,]
  df$Phenotype = phenotype
  return(df)
}

# read and process low-resolution phenotype data
read_KZ_data = function(data_dir, phenotypes){
  # read all phenotypes
  dfs_list = lapply(phenotypes, function(phenotype)(read_phenotype(phenotype, data_dir)))
  # combine into one data frame
  df = do.call("rbind", dfs_list)
  return(df)
}

# definition of informative points
get_informative_points = function(values, min_val_improvement = 0.01, min_k_improvement = 5){
  k_vals = c()
  m = length(values)
  while(TRUE){
    still_eligible = !logical(m)
    if(length(k_vals) > 0){
      ineligible_k = max(k + min_k_improvement - 1, max(which(values <= min_val + min_val_improvement)))
      still_eligible[1:ineligible_k] = FALSE
    }
    if(!any(still_eligible)){
      break
    }
    min_val = min(values[still_eligible])
    k = which(still_eligible)[which.min(values[still_eligible])]
    k_vals = c(k_vals, k)
  }
  return(k_vals)
}

# compute FDP-hat based on knockoff statistics W
get_FDP_hat = function(W){
  return((1+cumsum(W<0))/cumsum(W>0))
}

# compute FDP-bar based on knockoff statistics W
get_FDP_bar = function(W, C){
  return(pmin(1,floor(C*(1+cumsum(W<0)))/cumsum(W>0)))
}

# find number of rejections made for given phenotype 
# when controlling FDR at level q
get_FDR_thresh = function(df, phenotype, q){
  df_phenotype = df %>% filter(Phenotype == phenotype)
  eligible = df_phenotype$num_rejections[df_phenotype$FDP_hat <= q]
  if(length(eligible) > 0){
    return(max(eligible))
  } else{
    return(0)
  }
}

# find number of rejections made for given phenotype 
# when controlling FDP at level q
get_FDP_thresh = function(df, phenotype, q){
  df_phenotype = df %>% filter(Phenotype == phenotype)
  eligible = df_phenotype$num_rejections[df_phenotype$FDP_bar <= q]
  if(length(eligible) > 0){
    return(max(eligible))
  } else{
    return(0)
  }
}

# Read GREAT results
read_great_results = function(data_dir, k){
  great_results = read_tsv(sprintf("%s/GREAT_output/greatExportAll_%d.tsv", data_dir, k), 
                           skip = 3, n_max = 500, col_types = "-cc-d--d----------------")
  # great_results = great_results %>% filter(Desc %in% GO_terms) %>% select("ID", "Desc", "BinomP", "RegionFoldEnrich")
  great_results = great_results %>% select("ID", "Desc", "BinomP", "RegionFoldEnrich")
  great_results$set = k
  return(great_results)
}
