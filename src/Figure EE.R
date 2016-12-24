#### Figure EE ####
# y-axis is DHI
# x-axis is the offset between the seasonal peak and the timing of mass vaccination
# Different lines for amplitude = 1, 0.5, and 0

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load workspace ####
# load(file = "src/Figure_EE.RData")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
require(ggplot2)
library(data.table)

#### Set Static Conditions ####
years = 10
times <- seq(0,356*years)

max_V_months = 12*years
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

#### Changing Conditions ####
VE_conditions <- list(
  "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months),
  "Dukoral" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE, max_V_months = max_V_months)
)

beta_amp_conditions <- c(0, 0.5, 1)
beta_phase_shift_conditions <- seq(0, 365, length.out = 8)

sims <- length(VE_conditions) * length(beta_amp_conditions) * length(beta_phase_shift_conditions)

#### Initialize data frame ####
fig_EE_df <- data.table(VE_condition = rep(NA, sims), beta_amp_condition = rep(NA, sims), beta_phase_shift_condition = rep(NA, sims), Re = rep(NA, sims), prob_outbreak_10 = rep(NA, sims))

fig_EE_df$VE_condition_name <- rep(names(VE_conditions), each = sims/length(VE_conditions))
fig_EE_df$beta_amp_condition_name <- rep((beta_amp_conditions), times = sims/length(beta_amp_conditions))
fig_EE_df$beta_phase_shift_condition_name <- rep(rep((beta_phase_shift_conditions), each = length(beta_amp_conditions)), times = sims / (length(beta_amp_conditions)*length(beta_phase_shift_conditions))) 

fig_EE_df$VE_condition <- rep(VE_conditions, each = sims/length(VE_conditions))
fig_EE_df$beta_amp_condition <- rep(beta_amp_conditions, times = sims/length(beta_amp_conditions))
fig_EE_df$beta_phase_shift_condition <- rep(rep(beta_phase_shift_conditions, each = length(beta_amp_conditions)), times = sims / (length(beta_amp_conditions)*length(beta_phase_shift_conditions))) 

#### Loop through conditions ####
for (row in seq_len(sims)){
  
  # Create params
  params <- list(beta=0.75,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = fig_EE_df$beta_amp_condition[row],               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = fig_EE_df$beta_phase_shift_condition[row],          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=0*1/(365*40), # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
                 mig_in= 0,             # Rate of immigration
                 mig_out=0,             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE=fig_EE_df$VE_condition[row][[1]],                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                 vac_routine_freq = 0,          # Days between routine re-vaccination campaigns
                 vac_routine_frac = 0,          # Fraction of the population revaccinated during routine revaccination campaigns
                 vac_birth_frac = 0,            # Fraction of babies vaccinated
                 vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
                 vac_max = 5e50,                # Maximum number of vaccines to be given
                 vac_recip = c("all"),  # Recipients of vaccination ("all", "S", "migrant", "birth")
                 vac_stopper = 1e10          # Don't vaccinate after this day
  )
  inits = rep(0, 7+params$n.comps.V)
  inits[1] = 00000 # initially susceptible
  inits[2] = 100000 # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
  
  # Run model
  output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
  
  fig_EE_df$Re[row] <- list(output$Re)
  fig_EE_df$prob_outbreak_10[row] <- list(output$prob_outbreak_10)
  
  cat(".")
}

plot(fig_EE_df$Re[[1]], type = "l")
plot(fig_EE_df$prob_outbreak_10[[1]], type = "l")

#### Calculate DHI ####
fig_EE_df$DHI <- NA

for (row in 1:nrow(fig_EE_df)){
  DHI <- sum(fig_EE_df$Re[row][[1]] <= 1) - 1
  if (is.na(DHI)){
    fig_EE_df$DHI[row] <- 0
  } else{
    fig_EE_df$DHI[row] <- DHI
  }
}

#### Calculate days until Herd Immunity is first lost ####
fig_EE_df$HI_lost <- NA

for (row in 1:nrow(fig_EE_df)){
  # Find first day with herd immunity
  HI_gained <- which(fig_EE_df$Re[row][[1]] < 1)[1]
  HI_lost <- which(fig_EE_df$Re[row][[1]] >= 1)[1] - HI_gained + 1
  if (is.na(HI_lost)){
    fig_EE_df$HI_lost[row] <- 0
  } else{
    fig_EE_df$HI_lost[row] <- HI_lost
  }
}

#### Plot of DHI ####
ggplot(fig_EE_df, aes(x=beta_phase_shift_condition, y = HI_lost, color = factor(beta_amp_condition))) + geom_line()  + theme_bw() +
  facet_grid(.~VE_condition_name)

ggsave(file = "figures/Figure_EE.pdf", width = 4, height = 4, units = "in")

#### Save workspace ####
save.image(file = "src/Figure_EE.RData")