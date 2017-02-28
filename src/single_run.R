#### Test SIRV model with single conditions ####
setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load libraries and functions ####
source("src/calculate_Re.R")
source("src/calculate_VE.R")
source("src/Seasonality.R")
source("src/prob_outbreak_fcn.R")
source("src/migration_rate_calculator.R")
source("src/SIRV_model.R")
source("src/Run_SIRV_model.R")
source("src/revaccination.R")
require(ggplot2)

#### Example input parms #####
times <- seq(0,365*5)

# Calculate elements of VE
max_V_months = 48
V_comps_per_month = 0.5 # choose from 0.5, 1, 1.5, etc.

n.comps.V = max_V_months*V_comps_per_month
VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months)
# VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Dukoral",bound = TRUE,max_V_months = max_V_months)
# VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Perfect", bound = TRUE, max_V_months = max_V_months)

params <- list(beta=0.35,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4*1000,                   # Incubation period
               birth_death_rate=1/(365*40)*0,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= 1/(365*3.6)*0,             # Rate of immigration
               mig_out= 1/(365*4.6)*0,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = 365*0,         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 0,                 # Maximum number of vaccines to be given
               vac_recip = c("mass_S"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 1e10          # Don't vaccinate after this day
               
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 990 # initially susceptible
inits[2] = 0 # initially vaccinated
inits[params$n.comps.V+3] = 10 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

#### Test function ####
test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

#### Plot results ####
ggplot(test.run, aes(x = time/365, y = Re)) + geom_line() + geom_hline(yintercept=1, col="red") + 
  theme_bw() + xlab("Years") + ylab("R Effective") + scale_x_continuous(breaks = seq(0, 10, 1))

ggplot(test.run, aes(x = time/365)) +
  geom_line(aes(y=V_total), col="black", lty="dashed") +
  geom_line(aes(y=S), col="blue") +
  # geom_bar(aes(y=10*E), stat="identity", col="pink", alpha=0.1) +
  geom_bar(aes(y=10*I),stat="identity", col="darkred", alpha=0.1) +
  geom_line(aes(y=R), col="forestgreen") +
  theme_bw() + xlab("Years") + ylab("Number of People\nNote: 'I' are scaled by 10") + scale_x_continuous(breaks = seq(0, 10, 1))

ggplot(test.run, aes(x = time/365, y = prob_outbreak_10)) + geom_line() + theme_bw() + xlab("Years") + ylab("Probability of\nan Outbreak") + scale_x_continuous(breaks=seq(0,10,1)) + scale_y_continuous(limits = c(0,1))

ggplot(test.run, aes(x = time/365, y = Vax)) + geom_line() + theme_bw() + xlab("Years") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax)))

#### Summarize results ####
(return_to_R0 <- which(round(test.run$Re,3) >= 0.95*round(params$beta/params$gamma, 3))[1] / 365)
(DHI <- sum(round(test.run$Re,3) < 1)/365)


#### Plot simple epi curves ####
test.run_bigR <- test.run
test.run_smallR <- test.run

ggplot() +
  geom_line(data = test.run_bigR, aes(y=I, x = time/365), col="red") +
  geom_line(data = test.run_smallR, aes(y=I, x = time/365), col="blue") +
  geom_hline(yintercept = 10, col = "black", lty = "dashed") +
  theme_classic() + xlab("Time") + scale_x_continuous(breaks = seq(0, 1, 1), limits = c(0, 0.2)) +
  ylab("Cases") + scale_y_continuous(breaks = seq(0, 50, 10))
ggsave(file = "/Users/peakcm/Desktop/Epi_Curve.pdf", width = 3, height = 3, units = "in")
ggsave(file = "/Users/peakcm/Desktop/Epi_Curve_red.pdf", width = 3, height = 3, units = "in")
ggsave(file = "/Users/peakcm/Desktop/Epi_Curve_red_black.pdf", width = 3, height = 3, units = "in")


