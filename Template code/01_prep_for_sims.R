#' author: Eric R. Sokol <eric.r.sokol@colorado.edu>
#' created: 2022-08-24

# clean out workspace
rm(list=ls())
gc()

# load packages, fxns
library(tidyverse)
library(ade4)
source("00_FUNCTIONS.R")

#################################################
# USER INPUTS

# Set data source here
comm_infile <- "Diatoms - 400_L_300 Counts.csv"
env_infile <- "Diatoms - 400_L_300 Enviros.csv"
  
#################################################

# read in source data
comm_data_in <- read_csv(comm_infile) %>% as.data.frame()
env_data_in <- read_csv(env_infile) %>% as.data.frame()

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

# save formatted input data for sims
# as a .RDS
sim_input_data <- list(
  orig_source_file_names = list(
    comm_infile = comm_infile,
    env_infile = env_infile),
  comm_data = comm_data_counts,
  env_data = env_data_CST,
  geo_data = geo_data,
  niche_data = data_niche)
                       
saveRDS(sim_input_data, "sim_input_data.RDS")




