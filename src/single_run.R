#### Test SIRV model with single conditions ####
setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
require(ggplot2)

#### Example input parms #####
times <- seq(0,356*10)

# Calculate elements of VE
max_V_months = 48
V_comps_per_month = 1

n.comps.V = max_V_months*V_comps_per_month
# VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE)
VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE,max_V_months = max_V_months)
# VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect", bound = TRUE, max_V_months = max_V_months)

params <- list(beta=0.6538415,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=0*1/(365*40), # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_in= 1/(365*2),             # Rate of immigration
               mig_out=1/(365*2),             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5  # Average time in each vaccine compartment is one month
)
x = rep(0, 6+params$n.comps.V)
x[1] = 0000 # initially susceptible
x[2] = 100000 # initially vaccinated
x[params$n.comps.V+3] = 0 # initially infected

#### Test function ####
test.run <- run_model(func = SIRV.generic, times = times, params = params)

# Check to make sure number in pop is stationary
summary(apply(test.run[,c("S","E","I","R","V_total")], 1, sum))

#### Plot results ####
ggplot(test.run, aes(x = time/365, y = Re)) + geom_line() + geom_hline(yintercept=1, col="red") + 
  theme_bw() + xlab("Years") + ylab("R Effective") + scale_x_continuous(breaks = seq(0, 10, 1))
# ggplot(test.run, aes(x = time/365)) +
#   geom_line(aes(y=V_total), col="black", lty="dashed") +
#   geom_line(aes(y=S), col="blue") +
#   # geom_bar(aes(y=10*E), stat="identity", col="pink", alpha=0.1) +
#   geom_bar(aes(y=10*I),stat="identity", col="darkred", alpha=0.1) +
#   geom_line(aes(y=R), col="forestgreen") +
#   theme_bw() + xlab("Years") + ylab("Number of People\nNote: 'I' are scaled by 10") + scale_x_continuous(breaks = seq(0, 10, 1))

ggplot(test.run, aes(x = time/365, y = prob_outbreak_10)) + geom_line() + theme_bw() + xlab("Years") + ylab("Probability of\nan Outbreak") + scale_x_continuous(breaks=seq(0,10,1)) + scale_y_continuous(limits = c(0,1))

#### Summarize results ####
(return_to_R0 <- which(round(test.run$Re,3) >= 0.95*round(params$beta/params$gamma, 3))[1] / 365)
(DHI <- sum(round(test.run$Re,3) < 1)/365)
