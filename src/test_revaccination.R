#### Header ####
# Test revaccination function
# X-axis is n_v/S, or the fraction of S that is to be vaccinated that day. Log transformed axis
# Y-axis is the fraction of S actually moved V_1*/S
# Solid line is with the log inflation
# Dashed line is without the log inflation

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
source("src/revaccination.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
require(ggplot2)

#### Test Example ####
N = 1000 
n_v <- 100
revaccination(t = 10,
              N = N,
              S = N,
              V.states.vector = rep(0, 48),
              Vax = n_v,
              vac_routine_count = 0,
              vac_mass_freq = 10,
              vac_mass_frac = n_v/N,
              vac_mig_frac = 0,
              vac_recip = "mass_all",
              vac_max = n_v*100,
              birth_rate = 0,
              mig_in = 0,
              vac_birth_frac = 0)

revaccination(t = 10,
              N = N,
              S = N,
              V.states.vector = rep(0, 48),
              Vax = n_v,
              vac_routine_count = 0,
              vac_mass_freq = 10,
              vac_mass_frac = n_v/N,
              vac_mig_frac = 0,
              vac_recip = "mass_all",
              vac_max = n_v*100,
              birth_rate = 0,
              mig_in = 0,
              vac_birth_frac = 0,
              revaccination_log_trans = TRUE)

revaccination(t = 10,
              N = N,
              S = N,
              V.states.vector = rep(0, 48),
              Vax = n_v,
              vac_routine_count = 0,
              vac_mass_freq = 10,
              vac_mass_frac = n_v/N,
              vac_mig_frac = 0,
              vac_recip = "mass_all",
              vac_max = n_v*100,
              birth_rate = 0,
              mig_in = 0,
              vac_birth_frac = 0,
              revaccination_log_trans = FALSE)

#### Set parameters ####
n_v_options <- seq(0.01, 0.99, by = 0.01) * N
vac_mass_freq <- 10

df <- data.frame(n_v = rep(n_v_options, 2), log_trans = c(rep(TRUE, length(n_v_options)), rep(FALSE, length(n_v_options))), S_rate = NA, S_moved = NA)

#### Calculate S_rate for each n_v ####
for (row in 1:nrow(df)){
  out <- revaccination(t = 10,
                       N = N,
                       S = N,
                       V.states.vector = rep(0, 48),
                       Vax = df[row, "n_v"],
                       vac_routine_count = 0,
                       vac_mass_freq = vac_mass_freq,
                       vac_mass_frac = df[row, "n_v"]/N,
                       vac_mig_frac = 0,
                       vac_recip = "mass_all",
                       vac_max = df[row, "n_v"]*100,
                       birth_rate = 0,
                       mig_in = 0,
                       vac_birth_frac = 0,
                       revaccination_log_trans = df[row, "log_trans"])
    df[row, "S_rate"] <- out$mass_S_rate
}

ggplot(df, aes(x = n_v/N, y = S_rate, lty = log_trans)) + geom_line()

#### Test: Run model to calculate S_moved ####
max_V_months = 12*(4) 
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect",bound = TRUE,max_V_months = max_V_months)

N_pop = N
init_frac = 0
years = 15/365
times <- seq(0,356*years)

params <- list(beta=1.5/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.00,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=0, # Average birth and death rate
               nat_wane=0,         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= 0,             # Rate of immigration
               mig_out= 0,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE= VE,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/3000.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,          # Fraction of the current S pop that is vaccinated on each day
               vac_mass_freq = vac_mass_freq,         # Days between mass re-vaccination campaigns
               vac_mass_frac = .75,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
               vac_max = 1e20,                 # Maximum number of vaccines to be given
               vac_recip = "mass_all",  # Recipients of vaccination ("all", "S", "migrant", "birth")
               vac_stopper = 1e10,          # Don't vaccinate after this day
               revaccination_log_trans = TRUE 
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = N_pop*(1-init_frac) # initially susceptible
inits[2] = N_pop*init_frac # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

# Run model
output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
plot(output$S, type = "l")
output[10, "S"] - output[12, "S"]

#### Run model as a loop through all df conditions ####
# Note that errors arise for unknown reasons when the fraction to vaccination is between 0.75 and 0.99. Depends on what the N size is.... very odd....
for (row in 1:nrow(df)){
  params <- list(beta=1.5/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.00,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=0, # Average birth and death rate
                 nat_wane=0,         # Rate of natural immunity waning
                 mig_rates_constant = TRUE,      # TRUE if migration rates are constant
                 mig_in= 0,             # Rate of immigration
                 mig_out= 0,             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE= VE,                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/3000.5, # Average time in each vaccine compartment is one month
                 vac_routine_count = 0,          # Fraction of the current S pop that is vaccinated on each day
                 vac_mass_freq = 10,         # Days between mass re-vaccination campaigns
                 vac_mass_frac = df[row,"n_v"]/N,           # Fraction of the population revaccinated during mass revaccination campaigns
                 vac_birth_frac = 0,            # Fraction of babies vaccinated
                 vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
                 vac_max = 1e20,                 # Maximum number of vaccines to be given
                 vac_recip = "mass_all",  # Recipients of vaccination ("all", "S", "migrant", "birth")
                 vac_stopper = 1e10,          # Don't vaccinate after this day
                 revaccination_log_trans = df[row, "log_trans"] 
  )
  inits = rep(0, 7+params$n.comps.V)
  inits[1] = N_pop*(1-init_frac) # initially susceptible
  inits[2] = N_pop*init_frac # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
  
  # Run model
  output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
  df[row, "S_moved"] <- output[10, "S"] - output[12, "S"]
  
  cat(".")
}

df$log_trans <- factor(df$log_trans, levels = c(TRUE, FALSE), ordered = TRUE)
ggplot(df, aes(x = n_v/N, y = S_moved/N, lty = log_trans)) + 
  theme_classic() +
  stat_smooth(color = "black") +
  xlab("Desired Population Fraction to Vaccinate") +
  ylab("Actual Population Fraction Vaccinated") +
  scale_linetype(name = "Logarithmic\nAdjustment") +
  coord_fixed(1) +
  theme(legend.justification=c(1,0), legend.position=c(1,.1))

ggsave(file = "figures/Figure_II.pdf", width = 4, height = 4, units = "in")
