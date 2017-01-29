#### Figure BB ####
# Panels for VE waning only, VE and high migration, and VE and high migration and birth/death
# Show curves for R = 1, 1.5, 2
# Y axis is "Probability of an Outbreak

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load workspace ####
load(file = "src/Figure_BB.RData")

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

max_V_months = 12*(years+5)  # Go beyond "years" so the "perfect" vaccine remains perfect 
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

#### Changing Conditions ####
VE_conditions <- list(
  "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months),
  "Dukoral" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE, max_V_months = max_V_months),
  "Perfect" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect", bound = TRUE, max_V_months = max_V_months)
)

R0_conditions <- c("Low" = 1, "Moderate" = 1.5, "High" = 2)
mig_conditions <- c("low" = 1/(365*20), "high" = 1/(365*2))

sims <- length(VE_conditions) * length(R0_conditions) * length(mig_conditions)

#### Initialize data frame ####
fig_BB_df <- data.table(VE_condition = rep(NA, sims), R0_condition = rep(NA, sims), mig_condition = rep(NA, sims), Re = rep(NA, sims), prob_outbreak_10 = rep(NA, sims))

fig_BB_df$VE_condition_name <- rep(names(VE_conditions), each = sims/length(VE_conditions))
fig_BB_df$R0_condition_name <- rep(names(R0_conditions), times = sims/length(R0_conditions))
fig_BB_df$mig_condition_name <- rep(rep(names(mig_conditions), each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(mig_conditions))) 

fig_BB_df$VE_condition <- rep(VE_conditions, each = sims/length(VE_conditions))
fig_BB_df$R0_condition <- rep(R0_conditions, times = sims/length(R0_conditions))
fig_BB_df$mig_condition <- rep(rep(mig_conditions, each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(mig_conditions))) 

#### Loop through conditions ####
for (row in seq_len(sims)){
  
  # Create params
  params <- list(beta=fig_BB_df$R0_condition[row]/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                 beta_shape = "sinusoidal",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                 beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                 beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                 gamma=1/2,                     # Duration of disease
                 sigma=1/1.4,                   # Incubation period
                 birth_death_rate=1/(365*40), # Average birth and death rate
                 nat_wane=0*1/(365*10),         # Rate of natural immunity waning
                 mig_in= fig_BB_df$mig_condition[row],             # Rate of immigration
                 mig_out=fig_BB_df$mig_condition[row],             # Rate of emigration
                 foreign_infection=0.00,        # Proportion of immigrants who are infected
                 n.comps.V=n.comps.V,           # Number of V compartments
                 VE=fig_BB_df$VE_condition[row][[1]],                         # Vaccine efficacy over time
                 V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                 vac_routine_frac = 0.00,          # Fraction of the current S pop that is vaccinated on each day
                 vac_mass_freq = 0,         # Days between mass re-vaccination campaigns
                 vac_mass_frac = 0,           # Fraction of the population revaccinated during mass revaccination campaigns
                 vac_birth_frac = 0,            # Fraction of babies vaccinated
                 vac_mig_frac = 0,            # Fraction of immigrants vaccinated upon arrival
                 vac_max = 5e50,                 # Maximum number of vaccines to be given
                 vac_recip = c("all"),  # Recipients of vaccination ("all", "S", "migrant", "birth")
                 vac_stopper = 1e10          # Don't vaccinate after this day
  )
  inits = rep(0, 7+params$n.comps.V)
  inits[1] = 0000 # initially susceptible
  inits[2] = 100000 # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
  
  # Run model
  output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
  
  fig_BB_df$Re[row] <- list(output$Re)
  fig_BB_df$prob_outbreak_10[row] <- list(output$prob_outbreak_10)
  
  cat(".")
}

plot(fig_BB_df$Re[[1]], type = "l")
plot(fig_BB_df$prob_outbreak_10[[1]], type = "l")

#### Calculate DHI ####
fig_BB_df$DHI <- NA

for (row in 1:nrow(fig_BB_df)){
  fig_BB_df$DHI[row] <- sum(fig_BB_df$Re[row][[1]] <= 1) - 1
}

#### Manually melt data frame ####
fig_BB_df_melt <- data.frame(times = rep(times, sims), Re = NA, prob_outbreak_10 = NA,  VE_condition_name = NA, R0_condition_name = NA, mig_condition_name = NA, DHI = FALSE)

for (row in 1:sims){
  VE_condition_name = fig_BB_df$VE_condition_name[[row]]
  R0_condition_name = fig_BB_df$R0_condition_name[[row]]
  mig_condition_name = fig_BB_df$mig_condition_name[[row]]
  
  fig_BB_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("VE_condition_name")] <- VE_condition_name
  fig_BB_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("R0_condition_name")] <- R0_condition_name
  fig_BB_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("mig_condition_name")] <- mig_condition_name
  
  fig_BB_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("Re")] <- fig_BB_df$Re[[row]]
  fig_BB_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("prob_outbreak_10")] <- fig_BB_df$prob_outbreak_10[[row]]
  
  fig_BB_df_melt[fig_BB_df_melt$VE_condition_name %in% fig_BB_df$VE_condition_name[row] &
                  fig_BB_df_melt$R0_condition_name %in% fig_BB_df$R0_condition_name[row] &
                  fig_BB_df_melt$mig_condition_name %in% fig_BB_df$mig_condition_name[row] &
                  fig_BB_df_melt$times == fig_BB_df$DHI[row], "DHI"] <- TRUE
#     fig_BB_df_melt[fig_BB_df_melt$VE_condition_name %in% fig_BB_df$VE_condition_name[row] &
#                     fig_BB_df_melt$R0_condition_name %in% fig_BB_df$R0_condition_name[row] &
#                     fig_BB_df_melt$mig_condition_name %in% fig_BB_df$mig_condition_name[row] &
#                     fig_BB_df_melt$times == fig_BB_df$DHI[row], "prob_outbreak_10"]
  cat(".")
}

# Rename Factors
fig_BB_df_melt$VE_condition_name <- factor(fig_BB_df_melt$VE_condition_name, levels = c("Shanchol", "Dukoral", "Perfect"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)
fig_BB_df_melt$mig_condition_name <- factor(fig_BB_df_melt$mig_condition_name, levels = c("high", "low", "none"), labels = c("High", "Low", "None"), ordered = TRUE)
fig_BB_df_melt$R0_condition_name <- factor(fig_BB_df_melt$R0_condition_name, levels = c("High", "Moderate", "Low"), labels = c("High (2)", "Moderate (1.5)", "Low (1)"), ordered = TRUE)

#### Plot Re over time ####
# The blue dotted line essentially shows the assume vaccine waning profile.
ggplot(fig_BB_df_melt, aes(x = times/365, y = Re, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept =1, col = "darkgrey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Mass Vaccination") + ylab("Effective Reproductive Number") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_discrete(name = "In/Out Migration Rate") + ylim(0,2) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

# Plot figure JJ_2
fig_BB_df_melt_jj_2 <- fig_BB_df_melt
fig_BB_df_melt_jj_2$VE_condition_name <- fig_BB_df_melt_jj_2$VE_condition_name <- factor(fig_BB_df_melt_jj_2$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), ordered = TRUE)

jj_2 <- ggplot(fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$VE_condition_name %in% c("kOCV", "Perfect Vaccine"),], aes(x = times/365, y = Re, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept =1, col = "darkgrey") + geom_line() + facet_grid(.~VE_condition_name) + xlab("Years since Mass Vaccination") + ylab(expression(paste("Effective Reproductive Number ",(R[e])))) + theme_bw() + scale_color_discrete(name = expression(paste("Basic Reproductive Number ",(R[0])))) + scale_linetype_manual(name = "In/Out Migration Rate",values =  c("solid", "longdash")) + ylim(0,2) + theme(text = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + guides(lty = FALSE) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank()) + scale_x_continuous(name = NULL, labels = NULL, breaks = c(0,2,4,6,8,10)) + theme(strip.background = element_blank(),strip.text.x = element_blank())

plot(set_panel_size(p = jj_2, g = ggplotGrob(jj_2), margin = unit(.25,"in"), width=unit(2.5, "in"), height=unit(2, "in")))
set_panel_size(p = jj_2, g = ggplotGrob(jj_2), file = "figures/Figure_JJ_2.pdf", margin = unit(0.25,"in"), width=unit(2.5, "in"), height=unit(2, "in"))

#### Plot Prob of Outbreak ####
# Add an indicator on each curve marking the loss of Herd Immunity
# Alternatively, since the y=0.246 line is Re = 1, just draw a line
# Can add Re as a secondary axis with the probabilty of an outbreak

fig_BB_df_melt$VE_condition_name <- factor(fig_BB_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), ordered = TRUE)

ggplot(fig_BB_df_melt[fig_BB_df_melt$VE_condition_name %in% c("kOCV", "Perfect Vaccine"),], aes(x = times/365, y = prob_outbreak_10, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept = 0.246, col = "grey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Mass Vaccination") + ylab("Probability One Case sparks an Outbreak") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_manual(name = "In/Out Migration Rate", values = c("solid", "dashed")) + ylim(0,1) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "figures/Figure_BB.pdf", width = 5, height = 3, units = "in")

fig_BB_df_melt$VE_condition_name <- factor(fig_BB_df_melt$VE_condition_name, levels = c("Whole Cell\n(eg Shanchol)", "kOCV", "Perfect Vaccine"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)", "Perfect Vaccine"), ordered = TRUE)

# Plot figure JJ_3
fig_BB_df_melt_jj_3 <- fig_BB_df_melt_jj_2

jj_3 <- ggplot(fig_BB_df_melt_jj_3[fig_BB_df_melt_jj_3$VE_condition_name %in% c("kOCV", "Perfect Vaccine"),], aes(x = times/365, y = prob_outbreak_10, linetype = mig_condition_name, color = R0_condition_name)) + geom_line() + facet_grid(.~VE_condition_name) + xlab("Years since Mass Vaccination") + ylab("Probability one case\nsparks an outbreak (>10)") + theme_bw() + scale_color_discrete(name = expression(paste("Basic Reproductive Number ",(R[0])))) + scale_linetype_manual(name = "In/Out Migration Rate",values =  c("solid", "longdash")) + ylim(0,1) + theme(text = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10)) + guides(lty = FALSE, color = FALSE) + theme(panel.grid.major = element_blank()) + theme(panel.grid.minor = element_blank()) + scale_x_continuous(name = "Years since Mass Vaccination", breaks = c(0,2,4,6,8,10)) + theme(strip.background = element_blank(),strip.text.x = element_blank())

plot(set_panel_size(p = jj_3, g = ggplotGrob(jj_3), margin = unit(.25,"in"), width=unit(2.5, "in"), height=unit(2, "in")))
set_panel_size(p = jj_3, g = ggplotGrob(jj_3), file = "figures/Figure_JJ_3.pdf", margin = unit(0.25,"in"), width=unit(2.5, "in"), height=unit(2, "in"))

#### Calculate DHI for each condition ####
# kOCV, High mobility, High R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 1 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "High" & fig_BB_df_melt_jj_2$R0_condition_name == "High (2)",][1,"times"]/365
# kOCV, High mobility, Medium R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 1 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "High" & fig_BB_df_melt_jj_2$R0_condition_name == "Moderate (1.5)",][1,"times"]/365
# kOCV, High mobility, Low R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 0.99 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "High" & fig_BB_df_melt_jj_2$R0_condition_name == "Low (1)",][1,"times"]/365

# kOCV, Low mobility, High R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 1 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "Low" & fig_BB_df_melt_jj_2$R0_condition_name == "High (2)",][1,"times"]/365
# kOCV, Low mobility, Medium R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 1 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "Low" & fig_BB_df_melt_jj_2$R0_condition_name == "Moderate (1.5)",][1,"times"]/365
# kOCV, Low mobility, Low R0
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$Re > 0.99 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "Low" & fig_BB_df_melt_jj_2$R0_condition_name == "Low (1)",][1,"times"]/365

# When does the outbreak probability in a high mobility, high R0 setting exceed 0.5
fig_BB_df_melt_jj_2[fig_BB_df_melt_jj_2$prob_outbreak_10 > 0.50 & fig_BB_df_melt_jj_2$VE_condition_name == "kOCV" & fig_BB_df_melt_jj_2$mig_condition_name == "High" & fig_BB_df_melt_jj_2$R0_condition_name == "High (2)",][1,"times"]/365



#### Plot Prob of Outbreak for supplement ####
# Add an indicator on each curve marking the loss of Herd Immunity
# Alternatively, since the y=0.246 line is Re = 1, just draw a line
# Can add Re as a secondary axis with the probabilty of an outbreak

ggplot(fig_BB_df_melt, aes(x = times/365, y = prob_outbreak_10, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept = 0.246, col = "grey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years since Mass Vaccination") + ylab("Probability One Case sparks an Outbreak") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_manual(name = "In/Out Migration Rate", values = c("solid", "dashed")) + ylim(0,1) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "figures/Figure_BB_supplement.pdf", width = 5, height = 3, units = "in")
# ggsave(file = "figures/Figure_BB_seasonal.pdf", width = 5, height = 3, units = "in")

# #### Bar chart of DHI ####
# # A work in progress
# ggplot(fig_BB_df_melt[fig_BB_df_melt$DHI == 1,], aes(x = VE_condition_name, y = times/365, fill = R0_condition_name, color = mig_condition_name))  + geom_bar(position = "dodge", stat = "identity", width = .8) + ylab("Duration of Herd Immunity (Years)") + coord_flip() + scale_x_discrete(limits = rev(levels(fig_BB_df_melt$VE_condition_name)))


#### Save workspace ####
save.image(file = "src/Figure_BB.RData")
