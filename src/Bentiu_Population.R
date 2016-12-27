# Bentiu case study #
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
require(R0)

#### Load workspace ####
load(file = "src/Figure_GG.RData")

#### Load data ####
df <- read.csv("data/Bentiu.csv")
df_cases <- read.csv("data/Bentiu_cases.csv")

#### Clean data ####
df$Date <- as.Date(df$Date, format = "%m/%d/%y")
df_cases$Date <- as.Date(df_cases$Date, format = "%m/%d/%y")

#### Create weekly case data ####
df_cases$Cases_weekly <- 0
for (i in 1:nrow(df_cases)){
  if (i %in% seq(7, 700, by = 7)){
    df_cases[i, "Cases_weekly"] <- sum(df_cases[(i-6):i, "Cases"])
  }
  if (i == 31){
    df_cases[i, "Cases_weekly"] <- sum(df_cases[(i-2):i, "Cases"])
  }
}

#### Create rectangle ####
df_rect <- data.frame(xmin = as.Date("2016-10-01"), xmax = as.Date("2016-12-01"), ymin = 0, ymax = 50000)

#### Add fluxes ####
df$influx <- df$pop + df$Entries
df[df$Date == "2015-07-01", "influx"] <- NA  # Remove a weird observation
df$outflux <- df$pop - df$Exits

#### Plot Pop size, vaccines, and cases ####
ggplot(df) + 
  # geom_rect(data = df_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", color = "black") +
  # geom_bar(data = df_cases, aes(x = Date, y = Cases_weekly*500), stat = "identity", width = 2) + 
  geom_bar(aes(x = Date, y = Cholera.cases*1000), stat = "identity", fill = "lightgrey") + 
  geom_text(aes(x = Date, y = (Cholera.cases*1000), label = paste(Cholera.cases, "\nCases")), size = 2) + 
  geom_ribbon(data = df[is.na(df$influx)==0,], aes(x = Date, ymin = pop, ymax = influx), fill = "lightgreen") +
  geom_ribbon(data = df[is.na(df$outflux)==0,], aes(x = Date, ymin = outflux, ymax = pop), fill = "pink") +
  geom_line(aes(x = Date, y = pop), color = "darkblue") +
  geom_bar(data = df[is.na(df$Vaccine.Doses)==0,], aes(x = Date, y = Vaccine.Doses/2), stat = "identity", width = 6, color = "forestgreen", alpha = .6) +
  geom_text(data = df[is.na(df$Vaccine.Doses)==0,], aes(x = (Date + 90), y = (Vaccine.Doses/2 - 2000), label = paste(Vaccine.Doses, "\nDoses")), color = "forestgreen", size = 2) +
  scale_y_continuous(breaks = c(0, 5e4, 1e5, 1.5e5), labels = c("0", "50", "100", "150"), name = "Population (thousands)") +
  theme_bw() +
  theme(text = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))  +
  scale_x_date(date_labels = "%b '%y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Bentiu PoC Camp")

ggsave(file = "figures/Figure_GG.pdf", width = 4, height = 3, units = "in")


#### Run Observed simulation for Bentiu ####
# Time elapses from 
time_start <- as.numeric(as.Date("2014-06-01"))
time_end <- as.numeric(as.Date("2016-12-01"))
times <- seq(0,time_end-time_start)

# Calculate elements of VE
max_V_months = 48
V_comps_per_month = 0.5 # choose from 0.5, 1, 1.5, etc.

n.comps.V = max_V_months*V_comps_per_month
VE <- Create_VE(timesteps_per_month = V_comps_per_month, VE_shape = "Shanchol",bound = TRUE,max_V_months = max_V_months)

# Set migration rates
knot_1 <- as.numeric(as.Date("2015-12-01")) - time_start
knot_2 <- (as.numeric(as.Date("2016-05-01"))-time_start) - knot_1
base_rate <- 1/(365*4.3)
mig_in <- c(rep(1/(365*1.21)+base_rate, knot_1), rep(0+base_rate, knot_2), rep(0+base_rate, length(times)-knot_1-knot_2))
mig_out <-  c(rep(0, knot_1)+base_rate, rep(1/(365*1.21)+base_rate, knot_2), rep(0+base_rate, length(times)-knot_1-knot_2))

params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=1/(365*24.4),   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = FALSE,      # TRUE if migration rates are constant
               mig_in= mig_in,             # Rate of immigration
               mig_out= mig_out,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 174288,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 7310 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- seq(from = 4291, to = 40574, length.out = time_start - time_first)
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

ggplot() + geom_line(data = test.run, aes(x = time, y = N)) + theme_bw() + xlab("Years") + ylab("N")

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0.5,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggsave(file = "figures/Figure_GG_Xt.pdf", width = 4, height = 1.5, units = "in")

ggplot(test.run, aes(x = time, y = Re)) + geom_line() + geom_hline(yintercept=1, col="red") + 
  theme_bw() + xlab("Date") + ylab("R Effective") +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(test.run, aes(x = time, y = prob_outbreak_50)) + geom_line() + theme_bw() + xlab("Date") + ylab("Probability of\nan Outbreak (>50)") + scale_y_continuous(limits = c(0,1))+  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(test.run, aes(x = time)) +
  geom_line(aes(y=V_total), col="black", lty="dashed") +
  geom_line(aes(y=S), col="blue") +
  # geom_bar(aes(y=10*E), stat="identity", col="pink", alpha=0.1) +
  geom_bar(aes(y=10*I),stat="identity", col="darkred", alpha=0.1) +
  geom_line(aes(y=R), col="forestgreen") +
  theme_bw() + xlab("Date") + ylab("Number of People\nNote: 'I' are scaled by 10") +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Plot Pop size and simulated Pop size ####
ggplot(df) + 
  geom_vline(xintercept = as.numeric(as.Date("2014-06-01")), col = "lightblue", lty = "longdash") + 
  geom_vline(xintercept = as.numeric(as.Date("2015-12-01")), col = "lightblue", lty = "longdash") + 
  geom_vline(xintercept = as.numeric(as.Date("2016-05-01")), col = "lightblue", lty = "longdash") + 
  geom_line(aes(x = Date, y = pop), color = "black") +
  scale_y_continuous(breaks = c(0, 5e4, 1e5, 1.5e5), labels = c("0", "50", "100", "150"), name = "Population (thousands)") +
  theme_bw() +
  theme(text = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))  +
  scale_x_date(date_labels = "%b '%y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_line(data = test.run, aes(x = time, y = N), color = "blue", alpha = 0.5) +
  ggtitle("Bentiu PoC Camp")

ggsave(file = "figures/Figure_GG_Pop.pdf", width = 4, height = 3, units = "in")

#### Run Counterfactual simulation for Bentiu ####
VE_counterfactual <- rep(max(VE), length(VE))
base_rate_counterfactual <- 0
mig_in_counterfactual <- 0
mig_out_counterfactual <- 0
birth_death_rate_counterfactual <- 0

params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=birth_death_rate_counterfactual,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= mig_in_counterfactual,             # Rate of immigration
               mig_out= mig_out_counterfactual,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE_counterfactual,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 113426,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 66735 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- 100000
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

# Store with different R
test.run$R0 <- params$beta/params$gamma

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(test.run, aes(x = time, y = prob_outbreak_50)) + geom_line() + theme_bw() + xlab("Date") + ylab("Probability of\nan Outbreak (>50)") + scale_y_continuous(limits = c(0,1))+  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

prob_outbreak_fcn(R = 1, outbreak_size = 50)

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### Run Counterfactual simulation for Bentiu with empirical VE(t) ####
VE_counterfactual <- rep(max(VE), length(VE))
base_rate_counterfactual <- 0
mig_in_counterfactual <- 0
mig_out_counterfactual <- 0
birth_death_rate_counterfactual <- 0

params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=birth_death_rate_counterfactual,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= mig_in_counterfactual,             # Rate of immigration
               mig_out= mig_out_counterfactual,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 113426,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 66735 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Check to make sure enough vaccines were given (106,624)
test.run[nrow(test.run), "Vax"]

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- 100000
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

# Store with different R
test.run$R0 <- params$beta/params$gamma

test.run[test.run$time > "2016-10-16","Re"][1]

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### Run Counterfactual simulation for Bentiu with empirical N(t) ####
VE_counterfactual <- rep(max(VE), length(VE))
base_rate_counterfactual <- 0
mig_in <- c(rep(1/(365*1.21)+base_rate_counterfactual, knot_1), rep(0+base_rate_counterfactual, knot_2), rep(0+base_rate_counterfactual, length(times)-knot_1-knot_2))
mig_out <-  c(rep(0, knot_1)+base_rate_counterfactual, rep(1/(365*1.21)+base_rate_counterfactual, knot_2), rep(0+base_rate_counterfactual, length(times)-knot_1-knot_2))

birth_death_rate_counterfactual <- 0

params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=birth_death_rate_counterfactual,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = FALSE,      # TRUE if migration rates are constant
               mig_in= mig_in,             # Rate of immigration
               mig_out= mig_out,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE_counterfactual,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 113526,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 7310 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Check to make sure enough vaccines were given (106,624)
test.run[nrow(test.run), "Vax"]

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- 100000
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

# Store with different R
test.run$R0 <- params$beta/params$gamma

test.run[test.run$time > "2016-10-16","Re"][1]

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Run Counterfactual simulation for Bentiu with only birth/death ####
birth_death_rate <- 1/(365*24.4)
VE_counterfactual <- rep(max(VE), length(VE))
base_rate_counterfactual <- 0
mig_in <- 0
mig_out <- 0

birth_death_rate_counterfactual <- 0

params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=birth_death_rate,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= mig_in,             # Rate of immigration
               mig_out= mig_out,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE_counterfactual,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 113426,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 66735 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Check to make sure enough vaccines were given (106,624)
test.run[nrow(test.run), "Vax"]

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- 100000
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

# Store with different R
test.run$R0 <- params$beta/params$gamma

test.run[test.run$time > "2016-10-16","Re"][1]

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### Run Counterfactual simulation for Bentiu with only routine migration ####
params <- list(beta=1/2,                # Daily transmission parameter. From Guinea, beta=0.6538415
               beta_shape = "constant",       # Shape of the seasonal forcing function. "constant" or "sinusoidal"
               beta_amp = 0.05,               # Amplitude of sinusoidal seasonal forcing function (0 if no change, 1 if doubles)
               beta_phase_shift = 0,          # Phase shift in a sinusoidal seasonal forcing function
               gamma=1/2,                     # Duration of disease
               sigma=1/1.4,                   # Incubation period
               birth_death_rate=birth_death_rate_counterfactual,   # Average birth and death rate
               nat_wane=0*1/(365*10),         # Rate of natural immunity waning
               mig_rates_constant = TRUE,      # TRUE if migration rates are constant
               mig_in= base_rate,             # Rate of immigration
               mig_out= base_rate,             # Rate of emigration
               foreign_infection=0.00,        # Proportion of immigrants who are infected
               n.comps.V=n.comps.V,           # Number of V compartments
               VE=VE_counterfactual,                         # Vaccine efficacy over time
               V_step=V_comps_per_month/30.5, # Average time in each vaccine compartment is one month
               vac_routine_count = 0,        # Number of courses given each day
               vac_mass_freq = floor(365*11/12),         # Days between mass re-vaccination campaigns
               vac_mass_frac = 0.9,           # Fraction of the population revaccinated during mass revaccination campaigns
               vac_birth_frac = 0,            # Fraction of babies vaccinated
               vac_mig_frac = 0,              # Fraction of immigrants vaccinated upon arrival
               vac_max = 112526,                 # Maximum number of vaccine courses to be given
               vac_recip = c("mass_all"),     # Recipients of vaccination ("routine_S","routine_all", "mass_all", "mass_S", "migrant", "birth")
               vac_stopper = 370      # Don't vaccinate after this day
)
inits = rep(0, 7+params$n.comps.V)
inits[1] = 66735 # initially susceptible
inits[2] = 33265 # initially vaccinated
inits[params$n.comps.V+3] = 0 # initially infected
inits[7+params$n.comps.V] = inits[2] #Count those initially vaccinated in the Vax compartment

test.run <- run_model(inits = inits, func = SIRV.generic, times = times, params = params)
test.run$N <- apply(test.run[,c("S","E","I","R","V_total")], 1, sum)

# Check to make sure number in pop is stationary
summary(test.run$N)
# plot(apply(test.run[,c("S","E","I","R","V_total")], 1, sum), type = "l")
ggplot(test.run, aes(x = time/365, y = N)) + geom_line() + theme_bw() + xlab("Years") + ylab("N")

# Check to make sure enough vaccines were given (106,624)
test.run[nrow(test.run), "Vax"]

# Add the first few months pre-vaccination
time_first <- as.numeric(as.Date("2014-02-01"))

test.run$time <- test.run$time + time_start - time_first
added <- data.frame(matrix(nrow = time_start - time_first, ncol = ncol(test.run)))
names(added) <- names(test.run)
added$time <- 1:(time_start - time_first)
added[,c("S", "N")] <- 100000
added[,3:33] <- 0
added[,"Re"] <- params$beta/params$gamma
added[,"prob_outbreak_10"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 10)
added[,"prob_outbreak_50"] <- prob_outbreak_fcn(R = added[,"Re"], outbreak_size = 50)

test.run <- rbind(added, test.run)

test.run$time <- as.Date(test.run$time, origin = "2014-02-01")

# Store with different R
test.run$R0 <- params$beta/params$gamma

test.run[test.run$time > "2016-10-16","Re"][1]

# Plot results
ggplot(test.run, aes(x = time, y = Re/(params$beta/params$gamma))) + geom_line(alpha = 0.5) + 
  geom_ribbon(aes(x = time, ymin = Re/(params$beta/params$gamma), ymax = 1), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + xlab("Date") + ylab("X(t)") + ylim(0,1) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(panel.grid.minor = element_blank()) + theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(test.run, aes(x = time, y = Vax)) + geom_line() + theme_bw() + xlab("Date") + ylab("Number of Vaccine Courses Given") + scale_y_continuous(limits = c(0, max(test.run$Vax))) +  scale_x_date(date_labels = "%b '%y") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### Use Wallinga Teunis method to calculate Re using observed Bentiu case data ####
check.incid(df_cases$Cases, t=df_cases$Date)
plot(check.incid(df_cases$Cases, t=df_cases$Date)$incid, type = "b")

sa.GT(df_cases$Cases, GT.type="gamma", GT.mean=seq(3,10,1), GT.sd.range=1, est.method="EG")

#Time dependent method
mGT<-generation.time("gamma", c(5, 7.13))
plot(mGT,xlim = c(0, 30))
plot(cumsum(mGT$GT),xlim = c(1, 30), type = "b")
output <- est.R0.TD(df_cases$Cases, mGT,nsim=1000, begin=1, end=31)
plot(output)

df_cases$R <- output$R
df_cases$R_lower <- output$conf.int[,1]
df_cases$R_upper <- output$conf.int[,2]
df_cases$Date <- as.Date(as.character(df_cases$Date), format = "%m/%d/%y")

ggplot(df_cases, aes(x = Date)) +
  theme_bw() +
  geom_bar(aes(y = Cases/5), stat = "identity") +
  geom_hline(yintercept = 1, color = "red", lty = "longdash") +
  geom_ribbon(aes(ymin = R_lower, ymax = R_upper), fill = "pink", alpha = 0.5) +
  geom_line(aes( y = R), color = "darkred") +
  # scale_y_continuous(name = expression(R["t"])) +
  scale_y_continuous(name = "Daily Cholera Cases", breaks = c(0, 1, 2, 3, 4, 5), labels = c(0, 5, 10, 15, 20, 25), position = "right") +
  scale_x_date(date_labels = "%b %d") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "figures/Figure_GG_RtLeft.pdf", width = 4, height = 3, units = "in")
ggsave(file = "figures/Figure_GG_RtRight.pdf", width = 4, height = 3, units = "in")

#uses multiple methods
output2 <- estimate.R(df_cases$Cases, mGT, methods=c("EG", "ML", "TD", "AR", "SB"), pop.size=100000, begin=1, end=31, nsim=1000)
plot(output2)
output2

#### Save workspace ####
save.image(file = "src/Figure_GG.RData")
