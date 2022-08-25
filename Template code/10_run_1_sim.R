#' author: Eric R. Sokol <eric.r.sokol@colorado.edu>
#' created: 2022-08-24

# clean out workspace
rm(list=ls())
gc()

# load libraries and functions
library(tidyverse)
library(ade4)
source("00_FUNCTIONS.R")

#######################
# metacomm params

sim_m <- 0.10
sim_w <- 1000
sim_JM <- 1e6
sim_ntimesteps <- 50
sim_niche_mult <- 0.15
sim_nu <- 0.00001

sim_scenario_id <- "my_test_sim"
#######################

# read in source data
comm_data_in <- read_csv("Diatoms - 400_L_300 Counts.csv") %>% as.data.frame()
env_data_in <- read_csv("Diatoms - 400_L_300 Enviros.csv") %>% as.data.frame()

# maite sure the samples are ordered the same way
# get count data, remove taxa with all 0's
comm_data_counts <- comm_data_in %>% 
  arrange(SampleID) %>%
  select_if(function(col)is.numeric(col)&&(sum(col)>0))

# get the centered and standardized env data
env_data_CST <- env_data_in %>%
  arrange(SampleID) %>%
  select(ends_with("_CST"))

# get spatial data
geo_data <- env_data_in %>% 
  arrange(SampleID) %>%
  select(Longitude,Latitude)

# extract niche and trait info using ade4
data_niche <- get_niches(
  d.comm = comm_data_counts,
  d.env = env_data_CST)


# set random seed
set.seed(1)

# -- make the landscape
# -- I arbitrarily chose m = 0.05 and JM = 1e6
simulation.landscape <- MCSim::make.landscape(
  site.coords = geo_data,
  Ef = data_niche$site.niche$site.Ef,
  m = sim_m,
  JM = sim_JM)

# -- IMPORTANT NOTE on the simulation
# -- W.r = 5e6 produces a steep dispersal kernel such that dispersal 
# --    effectively only occurs between adjacent sites.
simoutput <- MCSim::metasim(
  landscape = simulation.landscape,
  output.dir.path = 'SIM_RESULTS',
  scenario.ID = sim_scenario_id,  
  trait.Ef = data_niche$species.niche$trait.Ef,
  trait.Ef.sd = sim_niche_mult * data_niche$species.niche$trait.Ef.sd.rescaled,
  J.t0 = comm_data_counts,
  n.timestep = sim_ntimesteps,
  W.r = sim_w,
  nu = sim_nu,
  save.sim = FALSE)



MCSim::plot.dot.plots(simoutput)
MCSim::plot.coenoclines(simoutput)
MCSim::plot.standardized.disp.kernel(simoutput)


comm_sim_wide <- simoutput$J.long %>%
  as.data.frame() %>%
  mutate(count = as.numeric(count)) %>%
  pivot_wider(id_cols = c("timestep","site"),
              names_from = spp,
              values_from = count)

comm_sim_last_t <- comm_sim_wide %>%
  filter(timestep == max(comm_sim_wide$timestep)) %>%
  select(-c(timestep,site))

sar_stats(comm_sim_last_t)

