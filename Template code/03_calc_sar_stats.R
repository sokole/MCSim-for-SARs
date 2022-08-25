#' author: Eric R. Sokol <eric.r.sokol@colorado.edu>
#' created: 2022-08-24

# clean out workspace
rm(list=ls())
gc()

# load libraries and functions
library(tidyverse)
source("00_FUNCTIONS.R")

# where is the simulation data?
sim_results_dir <- "SIM_RESULTS_Diatoms - 400_L_300 Counts"

# read in parameter mapping table
sim_parameter_map <- read_csv(paste0(sim_results_dir,"/sim_parameter_map.csv"))

# initialize sar stats table
sar_stats_tab <- data.frame()

for(i in 1:nrow(sim_parameter_map)){
  simoutput <- NULL
  comm_sim_wide <- NULL
  comm_sim_last_t <- NULL
  sar_stats_tab_i <- NULL
  try({
    # read in sim output
    simoutput <- readRDS(
      file = paste0(sim_results_dir,"/sim_",sim_parameter_map$sim_id[i],".RDS"))
    
    # extract community data and make wide
    comm_sim_wide <- simoutput$J.long %>%
      as.data.frame() %>%
      mutate(count = as.numeric(count)) %>%
      pivot_wider(id_cols = c("timestep","site"),
                  names_from = spp,
                  values_from = count)
    
    # filter to last timestep and only include the 
    # abundance data
    comm_sim_last_t <- comm_sim_wide %>%
      filter(timestep == max(comm_sim_wide$timestep)) %>%
      select(-c(timestep,site))
    
    # calc sar stats
    sar_stats_tab_i <- data.frame(
      sim_id = sim_parameter_map$sim_id[i],
      sar_stats(comm_sim_last_t))
    
    # add to sar_stats table
    sar_stats_tab <- bind_rows(
      sar_stats_tab, sar_stats_tab_i)
  })
  
  message(i, " / ",nrow(sim_parameter_map), " sar calcs complete")
}

# join with param table
sar_stats_tab_out <- left_join(
  sim_parameter_map,
  sar_stats_tab,
  by = "sim_id")

# write output
write_csv(
  sar_stats_tab_out,
  file = paste0(sim_results_dir,"/sar_stats.csv"))


# plot frequencies of best SAR mod types
sar_stats_tab_out %>% 
  ggplot(
    aes(BestSARModAICc, 
        fill = as.factor(niche_mult))) +
  geom_bar() +
  theme_bw()
  