#' author: Eric R. Sokol <eric.r.sokol@colorado.edu>
#' created: 2022-08-24

# clean out workspace
rm(list=ls())
gc()

# load libraries and functions
library(tidyverse)
library(ade4)
library(MCSim)
source("00_FUNCTIONS.R")

#######################
# USER INPUTS

# how many simulations to run
n_reps <- 10

# fixed sim parameters (same across all sims)
sim_JM <- 1e6 #metacommunity size
sim_ntimesteps <- 10 #generations - you want at least 50 generations
sim_nu <- 0.00001 #metacommunity invasion rate

# param ranges
sim_m <- 10^runif(n = n_reps, min = -3, max = -.1) #immigration rate
sim_w <- 10^runif(n = n_reps, min = 0, max = 9) #dispersal kernel slope
sim_niche_mult <- sample(c(.15, 1, 20), size = n_reps, replace = TRUE) #niche breadth multiplier

#######################

# read source data to constrain sim
sim_input_data <- readRDS("sim_input_data.RDS")
sim_scenario_id <- sim_input_data$orig_source_file_names$comm_infile %>% 
  gsub("(?i)\\.csv$","",.)

# where to write sims to
outdirname <- paste0('SIM_RESULTS_',sim_scenario_id)

# table to record sim parameters
sim_parameter_map <- data.frame()

# loop for each sim
for(i in 1:n_reps){
  try({
    simulation.landscape <- NULL
    simoutput <- NULL
    sim_parameter_map_i <- data.frame()
    
    # -- make the landscape
    # -- I arbitrarily chose m = 0.05 and JM = 1e6
    simulation.landscape <- MCSim::make.landscape(
      site.coords = sim_input_data$geo_data,
      Ef = sim_input_data$niche_data$site.niche$site.Ef,
      m = sim_m[i],
      JM = sim_JM)
    
    # -- IMPORTANT NOTE on the simulation
    # -- W.r = 5e6 produces a steep dispersal kernel such that dispersal 
    # --    effectively only occurs between adjacent sites.
    simoutput <- MCSim::metasim(
      landscape = simulation.landscape,
      output.dir.path = outdirname,
      scenario.ID = sim_scenario_id,
      trait.Ef = sim_input_data$niche_data$species.niche$trait.Ef,
      trait.Ef.sd = sim_niche_mult[i] * 
        sim_input_data$niche_data$species.niche$trait.Ef.sd.rescaled,
      J.t0 = sim_input_data$comm_data,
      n.timestep = sim_ntimesteps,
      W.r = sim_w[i],
      nu = sim_nu,
      save.sim = FALSE)
    
    # record parameter settings
    sim_parameter_map_i <- data.frame(
      rep = i,
      sim_scenario_id = sim_scenario_id,
      sim_id = simoutput$sim.result.name,
      m = sim_m[i],
      w = sim_w[i],
      niche_mult = sim_niche_mult[i],
      JM = sim_JM, 
      ntimesteps = sim_ntimesteps,
      nu = sim_nu)
    
    # add to mapping table
    sim_parameter_map <- bind_rows(
      sim_parameter_map, sim_parameter_map_i)
  })
  
  # save simulation as .RDS
  if(!is.null(simoutput)){
    saveRDS(
      simoutput,
      file = paste0(outdirname,"/sim_",simoutput$sim.result.name,".RDS"))}
  
  message(i, " / ", n_reps, " sims complete")
}

write_csv(
  sim_parameter_map,
  file = paste0(outdirname,"/sim_parameter_map.csv"))

