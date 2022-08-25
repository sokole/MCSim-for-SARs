#' author: Eric R. Sokol <eric.r.sokol@colorado.edu>
#' created: 2022-08-24

# functions for MCSim simulations for Passy et al.

# make niche info to use in sim
# based off example at 
# https://rpubs.com/sokole/MCSim-intro
get_niches <- function(
    d.comm, 
    d.env, 
    site.labels = NA_character_,
    niche.factor = 1, #rescale niche breadths
    scale.env.data = FALSE){
  
  # if not already scaled, center and scale the environmental data
  if(scale.env.data)  d.env <- scale(as.data.frame(d.env))
  if(is.na(site.labels)) site.labels <- row.names(d.env)
  
  # -- calculate RAs from densities
  d.comm.ra <- d.comm / rowSums(d.comm)
  
  # -- calculate niches for species
  dudi.pca.env <- dudi.pca(d.env, scale = TRUE, scan = FALSE, nf=1)
  niche.species <- niche(dudi.pca.env, Y = d.comm.ra, scann = FALSE)
  d.niche <- data.frame(
    niche.pos=niche.species$li,
    as.data.frame(niche.param(niche.species)))
  
  # -- calculate niche positions for each of the sites
  d.site.niche  <-  data.frame(
    site.label = site.labels,
    site.Ef = dudi.pca.env$li$Axis1)
  
  # -- trait matrix
  species.niche <- d.niche %>%
    select(Axis1, Tol) %>%
    rename(
      trait.Ef = Axis1,
      trait.Ef.sd = Tol) %>%
    mutate(
      trait.Ef.sd.rescaled = trait.Ef.sd * niche.factor)
  
  return(list(
    site.niche = d.site.niche,
    species.niche = species.niche))
}


# fxn nblike
nblike <- function(theta, n) #calculating aggregation "k", which is parameter from negative binomial distribution, measures aggregation.
{
  occi <- colSums(n > 0)
  ni=colSums(n[,which(colSums(n)>0)])
  k     <- theta[1]
  m     <-nrow(n)
  propi <- 1-(1+ni/(m*k))^-k
  logl  <- sum(occi*log(propi)+(m-occi)*log(1-propi))
  return(-logl)
}


# Joe's code to get SAR stats
# returns a 1-row table of stats used in Passy et al. paper

sar_stats <- function(
    comm_data_count = data.frame()){
  
  library(vegan)
  library(sars)
  library(AICcmodavg)
  


    # for(i in 1:length(Sample.Windows)){
    #Create the MASTER community dataset.  This is the primary metacommunity table.
  #comm_data_Envir    <- read.csv(paste0(Source1, Sample.Windows[i], " Enviros.csv"))
  ## The enviro medians have already been calculated and won't change. they aren't needed for 
  ## SAR fits
  
  lwindows_output <- NULL
  var.out         <- NULL
  var.out.long    <- NULL
  df.specaccum.FINAL <- NULL
  
  # scalewind <- scale.list[i]
  # Win.Now   <- win.list[i]   
  
  # comm.df     <- comm_data_count[,-1]
  # comm.df     <- comm.df[order(comm.df$SampleID),]
  # rownames(comm.df)  <- comm.df$SampleID
  # 
  # comm.list3  <- comm.df[,-1]
  # comm.list3  <- comm.list3[,(which(colSums(comm.list3) > 0))]
  
  # assume data has been sanitized to only have species counts (no sample info)
  comm.list3 <- comm_data_count
  #Diversity Measures -----
  # The full code is 3 nested for loops, which is what there is some funky indenting
  
  comm.mat<-as.data.frame(matrix(c(
    mean(specnumber(comm.list3)), # Mean Richness
    median(specnumber(comm.list3)), # Median Richness
    mean(na.omit(vegan::diversity(comm.list3)/log(specnumber(comm.list3)))),  #Mean Local Pielou J 
    median(na.omit(vegan::diversity(comm.list3)/log(specnumber(comm.list3)))), #Median Local level Pielou J evenness
    #Regional Diversity Measures
    specnumber(colSums(comm.list3)), #Reg rich
    vegan::diversity(colSums(comm.list3))/log(specnumber(colSums(comm.list3))), #Regional Pielou J
    (1/optim(0.1, nblike, n=comm.list3, method="BFGS")$par)), # Aggregation
    nrow=1)) 
  
  colnames(comm.mat) <- c("MeanRichness","MedianRichness",
                          "MeanJ","MedianJ",
                          "RegRich","RegJ", 
                          "Aggregation")
  
  ##SAC Model fits (Old SARs) ------
  
  df.specaccum      <- vegan::specaccum(comm.list3, method="exact") #Species accumulator function
  df.specaccum.out  <- as.data.frame(cbind(df.specaccum$sites, df.specaccum$richness)) #Expected Richness Outlet
  df.specaccum.out2 <- df.specaccum.out
  colnames(df.specaccum.out2) <- c("A", "S")
  df.sar.power      <- sar_power(df.specaccum.out)  #Actual SAR Model Fits
  df.sar.micmen     <- nls(V2 ~ SSmicmen(V1,VM,K), data=df.specaccum.out)
  df.sar.loga.LM    <- lm((S) ~ log10(A), data = df.specaccum.out2) ## More correct Logarithmic SAR JLM 9/8/21
  df.sar.power.coef <- df.sar.power$par
  
  ## Power SAC -----
  df.sar.power.coef.boPval <- summary(df.sar.power)$Parameters[1,4]
  df.sar.power.coef.b1Pval <- summary(df.sar.power)$Parameters[2,4]
  df.sar.power.aicc        <- summary(df.sar.power)$AICc
  S.Power       <- df.sar.power.coef[1] * df.specaccum.out2$A^df.sar.power.coef[2]
  Power.lmR     <- summary(lm(S.Power ~ (df.specaccum.out2$S)))
  Power.R2      <- Power.lmR$r.squared
  Power.R2adj   <- Power.lmR$adj.r.squared
  
  ## These are old SARs (SACs) The SAR. prefix is old. The true SAR version of this has .SAR suffix
  SAR.Loga.LM.coef        <- df.sar.loga.LM$coefficients
  SAR.Loga.LM.coef.boPval <- summary(df.sar.loga.LM)$coefficients[1,4]
  SAR.Loga.LM.coef.b1Pval <- summary(df.sar.loga.LM)$coefficients[2,4]
  SAR.Loga.LM.aicc        <- AICc(df.sar.loga.LM)
  SAR.Loga.LM.R2          <- summary(df.sar.loga.LM)$r.squared
  SAR.Loga.LM.R2a         <- summary(df.sar.loga.LM)$adj.r.squared
  
  # MM SAC ----
  df.sar.micmen.coef        <- coef(df.sar.micmen)
  df.sar.micmen.coef.boPval <- summary(df.sar.micmen)$parameters[1,4]
  df.sar.micmen.coef.b1Pval <- summary(df.sar.micmen)$parameters[2,4]
  df.sar.micmen.aicc        <- AICcmodavg::AICc(df.sar.micmen)
  S.MM       <- (df.sar.micmen.coef[1]*df.specaccum.out2$A) / (df.sar.micmen.coef[2] + df.specaccum.out2$A)
  MM.lmR     <- summary(lm(S.MM ~ (df.specaccum.out2$S)))
  MM.R2      <- MM.lmR$r.squared
  MM.R2adj   <- MM.lmR$adj.r.squared
  
  aicc.mods        <- c(SAR.Loga.LM.aicc, df.sar.power.aicc, df.sar.micmen.aicc)
  names(aicc.mods) <- c("L","P","M")
  del.aicc         <- aicc.mods - min(aicc.mods)
  champ.aicc       <- names(del.aicc)[which(del.aicc <= 2)]
  champ.aicc.mods  <- paste(champ.aicc, collapse="")
  
  # SAR Results Matrix
  
  sar.mat          <- as.data.frame(matrix(c(df.sar.power.coef, df.sar.power.coef.boPval, df.sar.power.coef.b1Pval,
                                             Power.R2, Power.R2adj,
                                             
                                             SAR.Loga.LM.coef, 
                                             SAR.Loga.LM.coef.boPval, SAR.Loga.LM.coef.b1Pval,
                                             SAR.Loga.LM.R2, SAR.Loga.LM.R2a,
                                             
                                             df.sar.micmen.coef, 
                                             df.sar.micmen.coef.boPval, df.sar.micmen.coef.b1Pval,
                                             MM.R2, MM.R2adj, 
                                             
                                             champ.aicc.mods), 
                                           nrow=1))
  
  colnames(sar.mat) <- c("SarPowerb0", "SarPowerb1", "SarPowerb0_Pval", "SarPowerb1_Pval",
                         "SarPowerR2", "SarPowerR2adj",
                         
                         "SAR.Loga.LMb0","SAR.Loga.LMb1","SAR.Loga.LMb0_Pval","SAR.Loga.LMb1_Pval",
                         "SAR.Loga.LMR2","SAR.Loga.LMR2adj",
                         
                         "SarMicmenb0","SarMicmenb1","SarMicmenb0_Pval","SarMicmenb1_Pval",
                         "SarMicmenR2","SarMicmenR2adj",
                         
                         "BestSARModAICc")
  
  out.vars.long <- as.data.frame(cbind(comm.mat, sar.mat))
  
  # write.csv(out.vars.long,  paste0(Resultz, Sample.Windows[i], " SAR Output Eric Addition.csv"))
  return(out.vars.long)
}
