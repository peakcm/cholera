#### Figure AA ####
# Panels for VE waning only, VE and high migration, and VE and high migration and birth/death
# Show curves for Shanchol, Dukoral, Exponential waning over 5 years
# Y axis is "X(t)"

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load workspace ####
load(file = "src/Figure_AA.RData")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
source("src/set_panel_size.R")
require(ggplot2)
require(data.table)

#### Set Static Conditions ####
years = 10
times <- seq(0,356*years)

max_V_months = 12*(years+5) # Go beyond "years" so the "perfect" vaccine remains perfect
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

#### Changing Conditions ####
VE_conditions <- list(
  "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months),
  "Dukoral" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE, max_V_months = max_V_months),
  "Perfect" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect", bound = TRUE, max_V_months = max_V_months)
)

birth_death_rate_conditions <- c("none" = 0, "high" = 1/(365*40))
mig_conditions <- c("none" = 0, "low" = 1/(365*20), "high" = 1/(365*2))

sims <- length(VE_conditions) * length(birth_death_rate_conditions) * length(mig_conditions)

#### Initialize data frame ####
fig_AA_df <- data.table(VE_condition = rep(NA, sims), birth_death_rate_condition = rep(NA, sims), mig_condition = rep(NA, sims), Re = rep(NA, sims), prob_outbreak_10 = rep(NA, sims))

fig_AA_df$VE_condition_name <- rep(names(VE_conditions), each = sims/length(VE_conditions))
fig_AA_df$birth_death_rate_condition_name <- rep(names(birth_death_rate_conditions), times = sims/length(birth_death_rate_conditions))
fig_AA_df$mig_condition_name <- rep(rep(names(mig_conditions), each = length(birth_death_rate_conditions)), times = sims / (length(birth_death_rate_conditions)*length(mig_conditions))) 

fig_AA_df$VE_condition <- rep(VE_conditions, each = sims/length(VE_conditions))
fig_AA_df$birth_death_rate_condition <- rep(birth_death_rate_conditions, times = sims/length(birth_death_rate_conditions))
fig_AA_df$mig_condition <- rep(rep(mig_conditions, each = length(birth_death_rate_conditions)), times = sims / (length(birth_death_rate_conditions)*length(mig_conditions))) 

#### Loop through conditions ####
for (row in seq_len(sims)){
  
  # Create params
  params <- list(beta=.5,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=fig_AA_df$birth_death_rate_condition[row], # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
                 mig_in= fig_AA_df$mig_condition[row],             # Rate of immigration
                 mig_out=fig_AA_df$mig_condition[row],             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE=fig_AA_df$VE_condition[row][[1]],                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                 vac_routine_freq = 0,          # Days between routine re-vaccination campaigns
                 vac_routine_frac = 0,          # Fraction of the population revaccinated during routine revaccination campaigns
                 vac_birth_frac = 0,            # Fraction of babies vaccinated
                 vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
                 vac_max = 5e50,                # Maximum number of vaccines to be given
                 vac_recip = c("all")  # Recipients of vaccination ("all", "S", "migrant", "birth")
  )
  inits = rep(0, 7+params$n.comps.V)
  inits[1] = 0000 # initially susceptible
  inits[2] = 100000 # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
  
  # Run model
  output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
  
  fig_AA_df$Re[row] <- list(output$Re)
  fig_AA_df$prob_outbreak_10[row] <- list(output$prob_outbreak_10)
  
  cat(".")
}

plot(fig_AA_df$Re[[1]], type = "l")
plot(fig_AA_df$prob_outbreak_10[[1]], type = "l")

#### Manually melt data frame ####
fig_AA_df_melt <- data.frame(times = rep(times, sims), Re = NA, prob_outbreak_10 = NA,  VE_condition_name = NA, birth_death_rate_condition_name = NA, mig_condition_name = NA)

for (row in 1:sims){
  VE_condition_name = fig_AA_df$VE_condition_name[[row]]
  birth_death_rate_condition_name = fig_AA_df$birth_death_rate_condition_name[[row]]
  mig_condition_name = fig_AA_df$mig_condition_name[[row]]
  
  fig_AA_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("VE_condition_name")] <- VE_condition_name
  fig_AA_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("birth_death_rate_condition_name")] <- birth_death_rate_condition_name
  fig_AA_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("mig_condition_name")] <- mig_condition_name
  
  fig_AA_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("Re")] <- fig_AA_df$Re[[row]]
  fig_AA_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("prob_outbreak_10")] <- fig_AA_df$prob_outbreak_10[[row]]
  
  cat(".")
}

# Rename Factors
fig_AA_df_melt$VE_condition_name <- factor(fig_AA_df_melt$VE_condition_name, levels = c("Shanchol", "Dukoral", "Perfect"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)
fig_AA_df_melt$mig_condition_name <- factor(fig_AA_df_melt$mig_condition_name, levels = c("high", "low", "none"), labels = c("High", "Low", "None"), ordered = TRUE)
fig_AA_df_melt$birth_death_rate_condition_name <- factor(fig_AA_df_melt$birth_death_rate_condition_name, levels = c("high", "none"), labels = c("High", "None"), ordered = TRUE)

#### Plot X(t) ####
# The no-birth, no-mig line shows the vaccine waning profile.
fig_AA_df_melt$VE_condition_name <- factor(fig_AA_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), ordered = TRUE)

ggplot(fig_AA_df_melt[fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$VE_condition_name %in% c("kOCV", "Perfect Vaccine"),], aes(x = times/365, y = Re, linetype = mig_condition_name)) + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Mass Vaccination") + ylab("X(t)") + theme_bw()  + scale_linetype_manual(name = "In/Out Migration Rate", values = c("solid", "longdash", "dotted")) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "figures/Figure_AA.pdf", width = 5, height = 3, units = "in")

fig_AA_df_melt$VE_condition_name <- factor(fig_AA_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)


#### Plot X(t) for multi-paneled Figure JJ ####
# The no-birth, no-mig line shows the vaccine waning profile.
fig_AA_df_melt$VE_condition_name <- factor(fig_AA_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), ordered = TRUE)

jj <- ggplot(fig_AA_df_melt[fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$VE_condition_name %in% c("kOCV", "Perfect Vaccine"),], aes(x = times/365, y = Re, linetype = mig_condition_name)) + geom_line() + facet_grid(.~VE_condition_name) +  ylab("X(t)") + theme_bw()  + scale_linetype_manual(name = "In/Out Migration Rate", values = c("solid", "longdash", "dotted")) + theme(text = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + scale_x_continuous(labels = NULL, breaks = c(0,2,4,6,8,10), name = NULL) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank())

plot(set_panel_size(p = jj, g = ggplotGrob(jj), margin = unit(.25,"in"), width=unit(2.5, "in"), height=unit(2, "in")))
set_panel_size(p = jj, g = ggplotGrob(jj), file = "figures/Figure_JJ_1.pdf", margin = unit(0.25,"in"), width=unit(2.5, "in"), height=unit(2, "in"))

fig_AA_df_melt$VE_condition_name <- factor(fig_AA_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)

#### Plot X(t) completely for supplement ####
# The no-birth, no-mig line shows the vaccine waning profile.
ggplot(fig_AA_df_melt, aes(x = times/365, y = Re, linetype = mig_condition_name, color = birth_death_rate_condition_name)) + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Mass Vaccination") + ylab("X(t)") + theme_bw() + scale_color_manual(values = c("black", "grey"), name = "Birth/Death Rate") + scale_linetype_manual(name = "In/Out Migration Rate", values = c("solid", "longdash", "dotted")) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "figures/Figure_AA_supplement.pdf", width = 5, height = 3, units = "in")

#### Losing herd immunity ####

# No mobility, Shanchol
fig_AA_df_melt[fig_AA_df_melt$VE_condition_name == "Whole Cell\n(eg Shanchol)" & fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$mig_condition_name == "None" & fig_AA_df_melt$Re > 0.99,"times"][1]/365

# High mobility, Shanchol
fig_AA_df_melt[fig_AA_df_melt$VE_condition_name == "Whole Cell\n(eg Shanchol)" & fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$mig_condition_name == "High" & fig_AA_df_melt$Re > 0.99,"times"][1]/365

# No mobility, Dukoral
fig_AA_df_melt[fig_AA_df_melt$VE_condition_name == "BS-Whole Cell\n(eg Dukoral)" & fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$mig_condition_name == "None" & fig_AA_df_melt$Re > 0.99,"times"][1]/365

# High mobility, Dukoral
fig_AA_df_melt[fig_AA_df_melt$VE_condition_name == "BS-Whole Cell\n(eg Dukoral)" & fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$mig_condition_name == "High" & fig_AA_df_melt$Re > 0.99,"times"][1]/365

# High mobility, Perfect vaccine
fig_AA_df_melt[fig_AA_df_melt$VE_condition_name == "Perfect Vaccine" & fig_AA_df_melt$birth_death_rate_condition_name == "None" & fig_AA_df_melt$mig_condition_name == "High" & fig_AA_df_melt$Re > 0.99,"times"][1]/365

#### Save workspace ####
save.image(file = "src/Figure_AA.RData")
