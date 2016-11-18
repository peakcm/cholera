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
N_vax <- N_pop*10

max_V_months = 12*(4)  # Go beyond "years" so the "perfect" vaccine remains perfect 
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

#### Changing Conditions ####
VE_conditions <- list(
  "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol", bound = TRUE, max_V_months = max_V_months)
)

R0_conditions <- c("Moderate" = 1.5, "High" = 2)
mig_conditions <- c("low" = 1/(365*20), "high" = 1/(365*2))

strategies <- c("Routine_1%", "Mass_100%_2yr","Mass_100%_routine_1%")

sims <- length(VE_conditions) * length(R0_conditions) * length(mig_conditions) * length(strategies)

#### Initialize data frame ####
fig_FF_df <- data.table(VE_condition = rep(NA, sims), R0_condition = rep(NA, sims), mig_condition = rep(NA, sims), strategy = rep(NA, sims), Re = rep(NA, sims), prob_outbreak_10 = rep(NA, sims))

fig_FF_df$VE_condition_name <- rep(names(VE_conditions), each = sims/length(VE_conditions))
fig_FF_df$R0_condition_name <- rep(names(R0_conditions), times = sims/length(R0_conditions))
fig_FF_df$mig_condition_name <- rep(rep(names(mig_conditions), each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(mig_conditions))) 
fig_FF_df$strategy_name <- rep(rep(strategies, each = length(R0_conditions)*length(VE_conditions)), times = sims / (length(VE_conditions)*length(R0_conditions)*length(strategies))) 

fig_FF_df$VE_condition <- rep(VE_conditions, each = sims/length(VE_conditions))
fig_FF_df$R0_condition <- rep(R0_conditions, times = sims/length(R0_conditions)) 
fig_FF_df$mig_condition <- rep(rep(mig_conditions, each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(mig_conditions))) 


#### Loop through conditions ####
for (row in seq_len(sims)){
  
  if (fig_FF_df[row]$strategy_name == "Routine_1%"){
    vac_routine_frac = 0.01
    vac_mass_frac = 0
    vac_mass_freq = 0
    vac_birth_frac = 0
    vac_mig_frac = 0
    vac_recip = c("routine_S")
    init_frac = 0
  } else if (fig_FF_df[row]$strategy_name == "Mass_100%_2yr"){
    vac_routine_frac = 0
    vac_mass_frac = 1
    vac_mass_freq = 365*2
    vac_birth_frac = 0
    vac_mig_frac = 0
    vac_recip = c("mass_S")
    init_frac = 1
  } else if (fig_FF_df[row]$strategy_name == "Mass_100%_routine_1%"){
    vac_routine_frac = 0.01
    vac_mass_frac = 0
    vac_mass_freq = 0
    vac_birth_frac = 0
    vac_mig_frac = 0
    vac_recip = c("routine_S")
    init_frac = 1
  }
  
  # Create params
  params <- list(beta=fig_FF_df$R0_condition[row]/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.00,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=1/(365*40), # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
                 mig_in= fig_FF_df$mig_condition[row],             # Rate of immigration
                 mig_out=fig_FF_df$mig_condition[row],             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE=fig_FF_df$VE_condition[row][[1]],                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                 vac_routine_frac = vac_routine_frac,          # Fraction of the current S pop that is vaccinated on each day
                 vac_mass_freq = vac_mass_freq,         # Days between mass re-vaccination campaigns
                 vac_mass_frac = vac_mass_frac,           # Fraction of the population revaccinated during mass revaccination campaigns
                 vac_birth_frac = vac_birth_frac,            # Fraction of babies vaccinated
                 vac_mig_frac = vac_mig_frac,            # Fraction of immigrants vaccinated upon arrival
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
  
  cat(".")
}

plot(fig_FF_df$Re[[1]], type = "l")
plot(fig_FF_df$prob_outbreak_10[[1]], type = "l")

#### Calculate DHI ####
fig_FF_df$DHI <- NA

for (row in 1:nrow(fig_FF_df)){
  fig_FF_df$DHI[row] <- sum(fig_FF_df$Re[row][[1]] <= 1) - 1
}

#### Manually melt data frame ####
fig_FF_df_melt <- data.frame(times = rep(times, sims), Re = NA, prob_outbreak_10 = NA,  VE_condition_name = NA, R0_condition_name = NA, mig_condition_name = NA, DHI = FALSE)

for (row in 1:sims){
  VE_condition_name = fig_FF_df$VE_condition_name[[row]]
  R0_condition_name = fig_FF_df$R0_condition_name[[row]]
  mig_condition_name = fig_FF_df$mig_condition_name[[row]]
  
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("VE_condition_name")] <- VE_condition_name
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("R0_condition_name")] <- R0_condition_name
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("mig_condition_name")] <- mig_condition_name
  
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("Re")] <- fig_FF_df$Re[[row]]
  fig_FF_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("prob_outbreak_10")] <- fig_FF_df$prob_outbreak_10[[row]]
  
  fig_FF_df_melt[fig_FF_df_melt$VE_condition_name %in% fig_FF_df$VE_condition_name[row] &
                   fig_FF_df_melt$R0_condition_name %in% fig_FF_df$R0_condition_name[row] &
                   fig_FF_df_melt$mig_condition_name %in% fig_FF_df$mig_condition_name[row] &
                   fig_FF_df_melt$times == fig_FF_df$DHI[row], "DHI"] <- TRUE
  #     fig_FF_df_melt[fig_FF_df_melt$VE_condition_name %in% fig_FF_df$VE_condition_name[row] &
  #                     fig_FF_df_melt$R0_condition_name %in% fig_FF_df$R0_condition_name[row] &
  #                     fig_FF_df_melt$mig_condition_name %in% fig_FF_df$mig_condition_name[row] &
  #                     fig_FF_df_melt$times == fig_FF_df$DHI[row], "prob_outbreak_10"]
  cat(".")
}

# Rename Factors
fig_FF_df_melt$VE_condition_name <- factor(fig_FF_df_melt$VE_condition_name, levels = c("Shanchol", "Dukoral", "Perfect"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)
fig_FF_df_melt$mig_condition_name <- factor(fig_FF_df_melt$mig_condition_name, levels = c("high", "low", "none"), labels = c("High", "Low", "None"), ordered = TRUE)
fig_FF_df_melt$R0_condition_name <- factor(fig_FF_df_melt$R0_condition_name, levels = c("High", "Moderate", "Low"), labels = c("High (2)", "Moderate (1.5)", "Low (1)"), ordered = TRUE)

#### Plot Re over time ####
# The blue dotted line essentially shows the assume vaccine waning profile.
ggplot(fig_FF_df_melt, aes(x = times/365, y = Re, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept =1, col = "darkgrey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Vaccination") + ylab("Effective Reproductive Number") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_discrete(name = "In/Out Migration Rate") + ylim(0,2) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))