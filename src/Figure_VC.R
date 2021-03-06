#### Figure VC on impact of vaccine coverage ####

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load workspace ####
load(file = "src/Figure_VC.RData")

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
require(RColorBrewer)

#### Set Static Conditions ####
years = 10
times <- seq(0,356*years)

max_V_months = 12*years
V_comps_per_month = 1
n.comps.V = max_V_months*V_comps_per_month

#### Changing Conditions ####
# VE_conditions <- list(
#   "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months),
#   "Dukoral" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE, max_V_months = max_V_months)
# )
VE_conditions <- list(
  "Shanchol" = Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months)
  )

R0_conditions <- seq(1, 3, 0.05)
VC_conditions <- seq(0.05, 1, 0.05)

sims <- length(VE_conditions) * length(R0_conditions) * length(VC_conditions)

#### Initialize data frame ####
fig_VC_df <- data.table(VE_condition = rep(NA, sims), R0_condition = rep(NA, sims), VC_condition = rep(NA, sims), Re = rep(NA, sims), prob_outbreak_10 = rep(NA, sims))

fig_VC_df$VE_condition_name <- rep(names(VE_conditions), each = sims/length(VE_conditions))
fig_VC_df$R0_condition_name <- rep((R0_conditions), times = sims/length(R0_conditions))
fig_VC_df$VC_condition_name <- rep(rep((VC_conditions), each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(VC_conditions))) 

fig_VC_df$VE_condition <- rep(VE_conditions, each = sims/length(VE_conditions))
fig_VC_df$R0_condition <- rep(R0_conditions, times = sims/length(R0_conditions))
fig_VC_df$VC_condition <- rep(rep(VC_conditions, each = length(R0_conditions)), times = sims / (length(R0_conditions)*length(VC_conditions))) 

#### Drop simulations I know won't have DHI ####
fig_VC_df$skip = 0
fig_VC_df[fig_VC_df$R0_condition > 1.6 & fig_VC_df$VC_condition < 0.5, "skip"] <- 1
fig_VC_df[fig_VC_df$R0_condition > 2 & fig_VC_df$VC_condition < 0.75, "skip"] <- 1

#### Loop through conditions ####
for (row in 1:nrow(fig_VC_df)){
  if (fig_VC_df[row, "skip"] != 1){
    # Create params
    params <- list(beta=fig_VC_df$R0_condition[row]/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
                   beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
                   beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
                   beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
                   gamma=1/2,                     # Duration of disease
                   sigma=1/1.4,                   # Incubation period
                   birth_death_rate=0, # Average birth and death rate
                   nat_wane=0,         # Rate of natural immunity waning
                   mig_rates_constant = TRUE,      # TRUE if migration rates are constant
                   mig_in= 0,             # Rate of immigration
                   mig_out=0,             # Rate of emigration
                   foreign_infection=0.00,        # Proportion of immigrants who are infected
                   n.comps.V=n.comps.V,           # Number of V compartments
                   VE=fig_VC_df$VE_condition[row][[1]],                         # Vaccine efficacy over time
                   V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
                   vac_routine_count = 0,          # Fraction of the current S pop that is vaccinated on each day
                   vac_mass_freq = 0,          # Days between routine re-vaccination campaigns
                   vac_mass_frac = 0,          # Fraction of the population revaccinated during routine revaccination campaigns
                   vac_birth_frac = 0,            # Fraction of babies vaccinated
                   vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
                   vac_max = 5e50,                # Maximum number of vaccines to be given
                   vac_recip = c("routine_S"),  # Recipients of vaccination ("all", "S", "migrant", "birth")
                   vac_stopper = 1e10          # Don't vaccinate after this day
    )
    inits = rep(0, 7+params$n.comps.V)
    inits[1] = 1000*(1-fig_VC_df$VC_condition[row]) # initially susceptible
    inits[2] = 1000*(fig_VC_df$VC_condition[row]) # initially vaccinated
    inits[params$n.comps.V+3] = 0 # initially infected
    inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment
    
    # Run model
    output <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
    
    fig_VC_df$Re[row] <- list(output$Re)
    fig_VC_df$prob_outbreak_10[row] <- list(output$prob_outbreak_10)

  }
  cat(".")
}

plot(fig_VC_df$Re[[1]], type = "l")
plot(fig_VC_df$prob_outbreak_10[[1]], type = "l")

#### Calculate DHI ####
fig_VC_df$DHI <- NA

for (row in 1:nrow(fig_VC_df)){
  DHI <- sum(fig_VC_df$Re[row][[1]] <= 1) - 1
  if (is.na(DHI)){
    fig_VC_df$DHI[row] <- 0
  } else{
    fig_VC_df$DHI[row] <- DHI
    }
}

#### Manually melt data frame ####
fig_VC_df_melt <- data.frame(times = rep(times, sims), Re = NA, prob_outbreak_10 = NA,  VE_condition_name = NA, R0_condition_name = NA, VC_condition_name = NA, DHI = FALSE)

for (row in 1:sims){
  VE_condition_name = fig_VC_df$VE_condition_name[[row]]
  R0_condition_name = fig_VC_df$R0_condition_name[[row]]
  VC_condition_name = fig_VC_df$VC_condition_name[[row]]
  
  fig_VC_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("VE_condition_name")] <- VE_condition_name
  fig_VC_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("R0_condition_name")] <- R0_condition_name
  fig_VC_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("VC_condition_name")] <- VC_condition_name
  
  fig_VC_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("Re")] <- fig_VC_df$Re[[row]]
  fig_VC_df_melt[(length(times)*(row-1)+1):(length(times)*row), c("prob_outbreak_10")] <- fig_VC_df$prob_outbreak_10[[row]]
  
  fig_VC_df_melt[fig_VC_df_melt$VE_condition_name %in% fig_VC_df$VE_condition_name[row] &
                  fig_VC_df_melt$R0_condition_name %in% fig_VC_df$R0_condition_name[row] &
                  fig_VC_df_melt$VC_condition_name %in% fig_VC_df$VC_condition_name[row] &
                  fig_VC_df_melt$times == fig_VC_df$DHI[row], "DHI"] <- TRUE
  #     fig_VC_df_melt[fig_VC_df_melt$VE_condition_name %in% fig_VC_df$VE_condition_name[row] &
  #                     fig_VC_df_melt$R0_condition_name %in% fig_VC_df$R0_condition_name[row] &
  #                     fig_VC_df_melt$mig_condition_name %in% fig_VC_df$mig_condition_name[row] &
  #                     fig_VC_df_melt$times == fig_VC_df$DHI[row], "prob_outbreak_10"]
  cat(".")
}

# Rename Factors
fig_VC_df_melt$VE_condition_name <- factor(fig_VC_df_melt$VE_condition_name, levels = c("Shanchol", "Dukoral"), labels = c("Whole Cell\n(eg Shanchol)", "BS-Whole Cell\n(eg Dukoral)"), ordered = TRUE)
# fig_VC_df_melt$mig_condition_name <- factor(fig_VC_df_melt$mig_condition_name, levels = c("high", "low", "none"), labels = c("High", "Low", "None"), ordered = TRUE)
# fig_VC_df_melt$R0_condition_name <- factor(fig_VC_df_melt$R0_condition_name, levels = c("High", "Moderate", "Low"), labels = c("High (2)", "Moderate (1.5)", "Low (1)"), ordered = TRUE)

#### Clean up workspace ####
rm(fig_VC_df)
rm(output)
fig_VC_df_melt <- fig_VC_df_melt[fig_VC_df_melt$DHI == 1,]

#### Heatmap of DHI ####
ggplot(fig_VC_df_melt[fig_VC_df_melt$R0_condition_name > 1 & fig_VC_df_melt$VE_condition_name == "Whole Cell\n(eg Shanchol)",], aes(x=VC_condition_name, y =R0_condition_name, fill = times/365)) + geom_tile()  + theme_bw() + ylab("Basic Reproductive Number") + xlab("Vaccine Coverage") + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6)) + scale_fill_gradientn(name="Duration of\nHerd Immunity\n(Years)", colours = rainbow(4))

# scale_fill_gradient2(name="Duration of\nHerd Immunity\n(Years)", low = "#80cdc1", mid = "white", midpoint = 2, high = "#a6611a") 

ggsave(file = "figures/Figure_VC.pdf", width = 3, height = 2, units = "in")

#### Heatmap of DHI for supplement ####
ggplot(fig_VC_df_melt[fig_VC_df_melt$R0_condition_name > 1,], aes(x=VC_condition_name, y =R0_condition_name, fill = times/365)) + geom_tile()  + theme_bw() + ylab("Basic Reproductive Number") + xlab("Vaccine Coverage") + facet_grid(VE_condition_name~.) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6)) + scale_fill_gradientn(name="Duration of\nHerd Immunity\n(Years)", colours = rainbow(4))

# scale_fill_gradient2(name="Duration of\nHerd Immunity\n(Years)", low = "#80cdc1", mid = "white", midpoint = 2, high = "#a6611a") 

ggsave(file = "figures/Figure_VC_supplement.pdf", width = 4, height = 4, units = "in")
# 
# #### Bar chart of DHI ####
# # A work in progress
# ggplot(fig_VC_df_melt[fig_VC_df_melt$DHI == 1,], aes(x = VE_condition_name, y = times/365, fill = R0_condition_name, color = mig_condition_name))  + geom_bar(position = "dodge", stat = "identity", width = .8) + ylab("Duration of Herd Immunity (Years)") + coord_flip() + scale_x_discrete(limits = rev(levels(fig_VC_df_melt$VE_condition_name)))

#### Save workspace ####
save.image(file = "src/Figure_VC.RData")
