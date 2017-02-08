setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load workspace ####
load(file = "src/Figure_FF_Sensitivity_Analysis.RData")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
source("src/revaccination.R")
require(ggplot2)
require(data.table)
require(deSolve)

#### Define Static Conditions ####
years = 20
# years = 20 # When using perfect vaccine
times <- seq(0,356*years)

N_pop <- 10000

#### Define Variable Conditions ####
N_vax <- N_pop*3

max_V_months = 12*(4) 
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

R0_conditions <- c(1.25, 1.5, 2, 2.5)

mig_conditions <- c(0, 1/(365*20), 1/(365*4.3), 1/(365*2), 1/(365*1))

VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE,max_V_months = max_V_months)
# VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect",bound = TRUE,max_V_months = max_V_months) # When using perfect vaccine

# Define a dynamic strategy for each whereby the mass vaccination coverage for mass_maintain is chosen using the R0 and the routine vaccination count is chosen by the R0, mig_rate, and birth rate.
# Work here.....

vac_mass_frac_conditions <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.4, 0.2, 0.1)     # This also applies to the initial mass coverage for "Mass" and "Mass_Maintain"
# vac_mass_frac_conditions <- c(1, 0.4)     # When using perfect vaccine

vac_mass_freq_conditions <- c(365, 365*2, 365*3) # Applies only to "Mass"
vac_routine_count_conditions <- c(0.5*round(N_vax/(365*years)),0.75*round(N_vax/(365*years)), round(N_vax/(365*years)),1.5*round(N_vax/(365*years)), 2*round(N_vax/(365*years)), 4*round(N_vax/(365*years))) # Applies to "Routine" and "Mass_Maintain"

strategies <- c("Routine", "Mass","Mass_Maintain")

routine_sims <- length(vac_routine_count_conditions) * length(mig_conditions) * length(R0_conditions)
mass_sims <- length(vac_mass_frac_conditions) * length(vac_mass_freq_conditions) * length(mig_conditions) * length(R0_conditions)
mass_maintain_sims <- length(vac_mass_frac_conditions) * length(vac_routine_count_conditions) * length(mig_conditions) * length(R0_conditions)

sims <- sum(routine_sims, mass_sims, mass_maintain_sims)

#### Initialize data frame ####
fig_FF_df <- data.table(strategy = rep(NA, sims), vac_mass_frac_conditions = NA, vac_mass_freq_conditions = NA, vac_routine_count_conditions = NA, mig_conditions = 999, R0_conditions = 999, Re = NA, prob_outbreak_10 = NA, DHI = NA, Vax = NA)

fig_FF_df$strategy <- c(rep(strategies[1], routine_sims), rep(strategies[2], mass_sims), rep(strategies[3], mass_maintain_sims))

fig_FF_df$vac_mass_frac_conditions <- c(rep(0, routine_sims), rep(vac_mass_frac_conditions, each = mass_sims/length(vac_mass_frac_conditions)), rep(vac_mass_frac_conditions, each = mass_maintain_sims/length(vac_mass_frac_conditions)))
fig_FF_df$vac_mass_freq_conditions <- c(rep(0, routine_sims), rep(vac_mass_freq_conditions, times = mass_sims/length(vac_mass_freq_conditions)), rep(0, mass_maintain_sims))
fig_FF_df$vac_routine_count_conditions <- c(rep(vac_routine_count_conditions, routine_sims/length(vac_routine_count_conditions)), rep(0, mass_sims), rep(vac_routine_count_conditions, times = mass_maintain_sims/length(vac_routine_count_conditions)))

# Add migration and R0 settings
fig_FF_df[fig_FF_df$strategy == "Routine",]$mig_conditions <- rep(mig_conditions, each = length(vac_routine_count_conditions)*length(R0_conditions))
fig_FF_df[fig_FF_df$strategy == "Routine",]$R0_conditions <- rep(rep(R0_conditions, each = length(vac_routine_count_conditions)), times = length(mig_conditions))

fig_FF_df[fig_FF_df$strategy == "Mass",]$mig_conditions <- rep(rep(mig_conditions, each = length(vac_mass_freq_conditions)*length(R0_conditions)), times =length(vac_mass_frac_conditions))
fig_FF_df[fig_FF_df$strategy == "Mass",]$R0_conditions <- rep(rep(rep(R0_conditions, each = length(vac_mass_freq_conditions)), times = length(mig_conditions)), times = length(vac_mass_frac_conditions))

fig_FF_df[fig_FF_df$strategy == "Mass_Maintain",]$mig_conditions <- rep(rep(mig_conditions, each = length(vac_routine_count_conditions)*length(R0_conditions)), times =length(vac_mass_frac_conditions))
fig_FF_df[fig_FF_df$strategy == "Mass_Maintain",]$R0_conditions <- rep(rep(rep(R0_conditions, each = length(vac_routine_count_conditions)), times = length(mig_conditions)), times = length(vac_mass_frac_conditions))

View(fig_FF_df)
#### Loop through conditions ####
for (row in seq_len(sims)){
  
  vac_routine_count = fig_FF_df$vac_routine_count_conditions[row]
  vac_mass_frac = fig_FF_df$vac_mass_frac_conditions[row]
  vac_mass_freq = fig_FF_df$vac_mass_freq_conditions[row]
  vac_mass_freq = fig_FF_df$vac_mass_freq_conditions[row]
  R0 = fig_FF_df$R0_conditions[row]
  init_frac = vac_mass_frac
  mig_in = fig_FF_df$mig_condition[row]
  mig_out = fig_FF_df$mig_condition[row]
  
  if (fig_FF_df[row]$strategy == "Routine"){
    vac_recip = c("routine_all")
  } else if (fig_FF_df[row]$strategy == "Mass"){
    vac_recip = c("mass_S")
  } else if (fig_FF_df[row]$strategy == "Mass_Maintain"){
    vac_recip = c("routine_all")
  }
  
  # Create params
  params <- list(beta=R0/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.00,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=1/(365*40), # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
                 mig_rates_constant = TRUE,      # TRUE if migration rates are constant
                 mig_in= mig_in,             # Rate of immigration
                 mig_out= mig_out,             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE= VE,                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                 vac_routine_count = vac_routine_count,          # Fraction of the current S pop that is vaccinated on each day
                 vac_mass_freq = vac_mass_freq,         # Days between mass re-vaccination campaigns
                 vac_mass_frac = vac_mass_frac,           # Fraction of the population revaccinated during mass revaccination campaigns
                 vac_birth_frac = 0,            # Fraction of babies vaccinated
                 vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
                 vac_max = N_vax,                 # Maximum number of vaccines to be given
                 vac_recip = vac_recip,  # Recipients of vaccination ("all", "S", "migrant", "birth")
                 vac_stopper = 1e10          # Don't vaccinate after this day
  )
  inits = rep(0, 7+params$n.comps.V)
  inits[1] = N_pop*(1-init_frac) # initially susceptible
  inits[2] = N_pop*init_frac # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
  
  # Run model
  output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
  
  fig_FF_df$Re[row] <- list(output$Re)
  fig_FF_df$prob_outbreak_10[row] <- list(output$prob_outbreak_10)
  fig_FF_df$Vax[row] <- list(output$Vax)
  
  fig_FF_df$DHI[row] <- (sum(fig_FF_df$Re[row][[1]] <= 1))/365
  
  cat("\n",row)
}

plot(fig_FF_df$Re[[1]], type = "l")
cat(N_vax - tail(fig_FF_df$Vax[[1]])[1], "remaining Vaccines")

# Add a marker for whether all vaccines were used in that trial or not
fig_FF_df$vax_rem <- 99999999
for (i in 1:nrow(fig_FF_df)){
  fig_FF_df[i]$vax_rem <- N_vax - tail(fig_FF_df$Vax[[i]])[1]
}

#### Find optimal strategy for each setting ####
settings <- data.frame(R0 = rep(R0_conditions, times = length(mig_conditions)), 
                       mig = rep(mig_conditions, each = length(R0_conditions)),
                       optimal_strategy = NA,
                       optimal_DHI = NA)
for (i in 1:nrow(settings)){
  mig <- settings[i, "mig"]
  R0 <- settings[i, "R0"]
  settings[i, "optimal_DHI"] <- max(fig_FF_df[fig_FF_df$mig_conditions == mig & fig_FF_df$R0_conditions == R0,]$DHI)
  settings[i, "optimal_strategy"] <- fig_FF_df[fig_FF_df$mig_conditions == mig & fig_FF_df$R0_conditions == R0 & fig_FF_df$DHI == settings[i, "optimal_DHI"],]$strategy[1]
}
settings

#### Inspect Re(t) for each setting ####
for (i in 1:nrow(settings)){
  mig <- settings[i, "mig"]
  R0 <- settings[i, "R0"]
  counter <- 1
  for (j in 1:nrow(fig_FF_df)){
    if (fig_FF_df[j,]$mig_conditions == mig & fig_FF_df[j,]$R0_conditions == R0){
      if (counter == 1){
        plot(fig_FF_df[j,]$Re[[1]], main = paste(mig, " ", R0), type = "l", ylim = c(0, 3))
      } else {
        lines(fig_FF_df[j,]$Re[[1]], type = "l")
      }
      counter = counter + 1
    }
  }
  readline(prompt="Press [enter] to continue")
}

fig_FF_df[fig_FF_df$mig_conditions == mig & fig_FF_df$R0_conditions == R0 & fig_FF_df$DHI > 5,]

#### Reduce file size ####
fig_FF_df$keep <- 0
for (i in 1:nrow(settings)){
  mig <- settings[i, "mig"]
  R0 <- settings[i, "R0"]

  max <- sort(fig_FF_df[fig_FF_df$mig_conditions == mig & 
              fig_FF_df$R0_conditions == R0 &
              fig_FF_df$strategy == "Routine",]$DHI, decreasing = TRUE)[1]
  which_max <- which(fig_FF_df[fig_FF_df$mig_conditions == mig & 
                    fig_FF_df$R0_conditions == R0 &
                    fig_FF_df$strategy == "Routine",]$DHI > floor(max))
  fig_FF_df[fig_FF_df$mig_conditions == mig & 
              fig_FF_df$R0_conditions == R0 &
              fig_FF_df$strategy == "Routine",]$keep[which_max] <- 1
  
  max <- sort(fig_FF_df[fig_FF_df$mig_conditions == mig & 
                          fig_FF_df$R0_conditions == R0 &
                          fig_FF_df$strategy == "Mass_Maintain",]$DHI, decreasing = TRUE)[1]
  which_max <- which(fig_FF_df[fig_FF_df$mig_conditions == mig & 
                                 fig_FF_df$R0_conditions == R0 &
                                 fig_FF_df$strategy == "Mass_Maintain",]$DHI > floor(max))
  fig_FF_df[fig_FF_df$mig_conditions == mig & 
              fig_FF_df$R0_conditions == R0 &
              fig_FF_df$strategy == "Mass_Maintain",]$keep[which_max] <- 1
  
  max <- sort(fig_FF_df[fig_FF_df$mig_conditions == mig & 
                          fig_FF_df$R0_conditions == R0 &
                          fig_FF_df$strategy == "Mass",]$DHI, decreasing = TRUE)[1]
  which_max <- which(fig_FF_df[fig_FF_df$mig_conditions == mig & 
                                 fig_FF_df$R0_conditions == R0 &
                                 fig_FF_df$strategy == "Mass",]$DHI > floor(max))
  fig_FF_df[fig_FF_df$mig_conditions == mig & 
              fig_FF_df$R0_conditions == R0 &
              fig_FF_df$strategy == "Mass",]$keep[which_max] <- 1
  cat(".")
}

fig_FF_df <- fig_FF_df[fig_FF_df$keep == 1,]

rm("output")

#### Save workspace ####
save.image(file = "src/Figure_FF_Sensitivity_Analysis.RData")
