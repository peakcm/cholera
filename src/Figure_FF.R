#### Figure FF ####
# Show three different vaccination strategies
# 1. mass vaccination only
# 2. routine vaccination only
# 3. Mass then Maintain
# Fix the vaccine courses to a particular number (eg 2*N) and Measure how long herd immunity is maintained
# Or, set the goal of maintaining herd immunity for 10 years and see how many doses you need to acheive it
  # For routine vaccination, start counting the clock when you acheive herd immunity
# Does the % difference (in time or number of vaccines needed) depend on R_0 or migration rates?
  # Does the choice in strategy depend on the R_0 or migration rates?
  # How do we estimate the number of vaccines needed in the initial mass campaign?
    # Assume mass campaigns always target 100% of the suscpetibles?

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

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

#### Set Static Conditions ####
years = 10
times <- seq(0,356*years)

N_pop <- 10000
N_vax <- N_pop*3

max_V_months = 12*(4) 
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months)

#### Changing Conditions ####
vac_mass_frac_conditions <- c(1, 0.8, 0.7)     # This also applies to the initial mass coverage for "Mass" and "Mass_Maintain"
vac_mass_freq_conditions <- c(365) # Applies only to "Mass"
vac_routine_count_conditions <- c(round(N_vax/(365*years)), 1.5*round(N_vax/(365*years)), 2*round(N_vax/(365*years))) # Applies to "Routine" and "Mass_Maintain"

# mig_conditions <- c("low" = 1/(365*20), "moderate" = 1/(365*4.3), "high" = 1/(365*2)) # Applies to all
mig_conditions <- c("high" = 1/(365*2))

strategies <- c("Routine", "Mass","Mass_Maintain")

routine_sims <- length(vac_routine_count_conditions)*length(mig_conditions)
mass_sims <- length(vac_mass_frac_conditions) * length(vac_mass_freq_conditions) * length(mig_conditions)
mass_maintain_sims <- length(vac_mass_frac_conditions) * length(vac_routine_count_conditions) * length(mig_conditions)

sims <- sum(routine_sims, mass_sims, mass_maintain_sims)

#### Initialize data frame ####
fig_FF_df <- data.table(strategy = rep(NA, sims), vac_mass_frac_conditions = NA, vac_mass_freq_conditions = NA, vac_routine_count_conditions = NA, mig_conditions = NA, Re = NA, prob_outbreak_10 = NA, DHI = NA, Vax = NA)

fig_FF_df$strategy <- c(rep(strategies[1], routine_sims), rep(strategies[2], mass_sims), rep(strategies[3], mass_maintain_sims))

fig_FF_df$vac_mass_frac_conditions <- c(rep(0, routine_sims), rep(vac_mass_frac_conditions, each = mass_sims/length(vac_mass_frac_conditions)), rep(vac_mass_frac_conditions, each = mass_maintain_sims/length(vac_mass_frac_conditions)))
fig_FF_df$vac_mass_freq_conditions <- c(rep(0, routine_sims), rep(vac_mass_freq_conditions, times = mass_sims/length(vac_mass_freq_conditions)), rep(0, mass_maintain_sims))
fig_FF_df$vac_routine_count_conditions <- c(rep(vac_routine_count_conditions, routine_sims/length(vac_routine_count_conditions)), rep(0, mass_sims), rep(vac_routine_count_conditions, times = mass_maintain_sims/length(vac_routine_count_conditions)))
  
fig_FF_df$mig_conditions <- rep(mig_conditions, sims)


#### Loop through conditions ####
for (row in seq_len(sims)){
  
    vac_routine_count = fig_FF_df$vac_routine_count_conditions[row]
    vac_mass_frac = fig_FF_df$vac_mass_frac_conditions[row]
    vac_mass_freq = fig_FF_df$vac_mass_freq_conditions[row]
    init_frac = vac_mass_frac
    mig_in = fig_FF_df$mig_condition[row]
    mig_out = fig_FF_df$mig_condition[row]
    
  if (fig_FF_df[row]$strategy == "Routine"){
    vac_recip = c("routine_S")
  } else if (fig_FF_df[row]$strategy == "Mass"){
    vac_recip = c("mass_S")
  } else if (fig_FF_df[row]$strategy == "Mass_Maintain"){
    vac_recip = c("routine_S")
  }
  
  # Create params
  params <- list(beta=1.5/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.00,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=1/(365*40), # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
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
                 vac_recip = vac_recip  # Recipients of vaccination ("all", "S", "migrant", "birth")
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
  
  fig_FF_df$DHI[row] <- (sum(fig_FF_df$Re[row][[1]] <= 1) - 1)/365
  
  cat(".")
}

plot(fig_FF_df$Re[[1]], type = "l")
plot(fig_FF_df$prob_outbreak_10[[1]], type = "l")

#### Manually melt to a data frame ####
fig_FF_df_melt <- data.frame(times = rep(times, sims), Re = NA, prob_outbreak_10 = NA,  strategy = NA, vac_mass_frac_condition = NA, vac_mass_freq_condition = NA, vac_routine_count_condition = NA, mig_condition = NA, DHI = FALSE)

for (row in 1:sims){
  strategy = fig_FF_df$strategy[[row]]
  vac_mass_frac_condition = fig_FF_df$vac_mass_frac_conditions[[row]]
  vac_mass_freq_condition = fig_FF_df$vac_mass_freq_conditions[[row]]
  vac_routine_count_condition = fig_FF_df$vac_routine_count_conditions[[row]]
  mig_condition = fig_FF_df$mig_conditions[[row]]
  
  
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("strategy")] <- strategy
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("vac_mass_frac_condition")] <- vac_mass_frac_condition
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("vac_mass_freq_condition")] <- vac_mass_freq_condition
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("vac_routine_count_condition")] <- vac_routine_count_condition
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("mig_condition")] <- mig_condition
  
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("Re")] <- fig_FF_df$Re[[row]]
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("prob_outbreak_10")] <- fig_FF_df$prob_outbreak_10[[row]]
  
  fig_FF_df_melt[fig_FF_df_melt$strategy %in% fig_FF_df$strategy[row] &
                   fig_FF_df_melt$vac_mass_frac_condition %in% fig_FF_df$vac_mass_frac_conditions[row] &
                   fig_FF_df_melt$vac_mass_freq_condition %in% fig_FF_df$vac_mass_freq_conditions[row] &
                   fig_FF_df_melt$vac_routine_count_condition %in% fig_FF_df$vac_routine_count_conditions[row] &
                   fig_FF_df_melt$mig_condition %in% fig_FF_df$mig_conditions[row] &
                   fig_FF_df_melt$times == fig_FF_df$DHI[row], "DHI"] <- TRUE
  #     fig_FF_df_melt[fig_FF_df_melt$VE_condition_name %in% fig_FF_df$VE_condition_name[row] &
  #                     fig_FF_df_melt$R0_condition_name %in% fig_FF_df$R0_condition_name[row] &
  #                     fig_FF_df_melt$mig_condition_name %in% fig_FF_df$mig_condition_name[row] &
  #                     fig_FF_df_melt$times == fig_FF_df$DHI[row], "prob_outbreak_10"]
  cat(".")
}

# # Rename Factors
# fig_FF_df_melt$VE_condition_name <- factor(fig_FF_df_melt$VE_condition_name, levels = c("Shanchol", "Dukoral", "Perfect"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)
# fig_FF_df_melt$mig_condition_name <- factor(fig_FF_df_melt$mig_condition_name, levels = c("high", "low", "none"), labels = c("High", "Low", "None"), ordered = TRUE)
# fig_FF_df_melt$R0_condition_name <- factor(fig_FF_df_melt$R0_condition_name, levels = c("High", "Moderate", "Low"), labels = c("High (2)", "Moderate (1.5)", "Low (1)"), ordered = TRUE)

#### Plot Re over time ####
# The blue dotted line essentially shows the assume vaccine waning profile.
ggplot(fig_FF_df_melt[fig_FF_df_melt$vac_mass_freq_condition != 730 & fig_FF_df_melt$vac_mass_frac_condition %in% c(0, 0.8, 1),], aes(x = times/365, y = Re, color = factor(vac_routine_count_condition), linetype = factor(vac_mass_frac_condition ))) + geom_hline(yintercept = 1, col = "darkgrey") + geom_line() + facet_grid(strategy ~.) + xlab("Years since Vaccination") + ylab("Effective Reproductive Number") + theme_bw() + ylim(0,2) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6)) + scale_linetype_manual(name = "Mass Vaccination Coverage", values = c("dotted", "longdash", "solid")) + scale_color_discrete(name = "Routine Doses per Day") 

ggsave(file = "figures/Figure_FF.pdf", width = 5, height = 3, units = "in")
 