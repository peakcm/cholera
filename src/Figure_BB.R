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
require(ggplot2)
require(data.table)

#### Set Static Conditions ####
years = 10
times <- seq(0,356*years)

max_V_months = 12*years
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
                 beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
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
                 V_step=V_comps_per_month/30.5  # Average time in each vaccine compartment is one month
  )
  inits = rep(0, 6+params$n.comps.V)
  inits[1] = 0000 # initially susceptible
  inits[2] = 100000 # initially vaccinated
  inits[params$n.comps.V+3] = 0 # initially infected
  
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
ggplot(fig_BB_df_melt, aes(x = times/365, y = Re, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept =1, col = "darkgrey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years") + ylab("Effective Reproductive Number") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_discrete(name = "In/Out Migration Rate") + ylim(0,2) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

#### Plot Prob of Outbreak ####
# Add an indicator on each curve marking the loss of Herd Immunity
# Alternatively, since the y=0.246 line is Re = 1, just draw a line
# Can add Re as a secondary axis with the probabilty of an outbreak
ggplot(fig_BB_df_melt, aes(x = times/365, y = prob_outbreak_10, linetype = mig_condition_name, color = R0_condition_name)) + geom_hline(yintercept = 0.246, col = "grey") + geom_line() + facet_grid(VE_condition_name ~.) + xlab("Years") + ylab("Probability of an Outbreak\n(at least 10 cases)") + theme_bw() + scale_color_discrete(name = "Basic Reproductive\nNumber") + scale_linetype_discrete(name = "In/Out Migration Rate") + ylim(0,1) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6))

ggsave(file = "figures/Figure_BB.pdf", width = 5, height = 3, units = "in")
# ggsave(file = "figures/Figure_BB_seasonal.pdf", width = 5, height = 3, units = "in")


#### Bar chart of DHI ####
# A work in progress
ggplot(fig_BB_df_melt[fig_BB_df_melt$DHI == 1,], aes(x = VE_condition_name, y = times/365, fill = R0_condition_name, color = mig_condition_name))  + geom_bar(position = "dodge", stat = "identity", width = .8) + ylab("Duration of Herd Immunity (Years)") + coord_flip() + scale_x_discrete(limits = rev(levels(fig_BB_df_melt$VE_condition_name)))


#### Save workspace ####
save.image(file = "src/Figure_BB.RData")
