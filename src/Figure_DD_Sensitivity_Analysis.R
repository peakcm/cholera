#### Header ####
# Sensitivity Analysis to see how optimal migration rate changes with input conditions
setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
source("src/vax_targeting.R")
source("src/prob_outbreak_fcn.R")
source("src/calculate_VE.R")
library(ggplot2)
library(lattice)
library(RColorBrewer) 
library(scales) # for muted colors

#### Load workspace ####
load(file = "src/Figure_DD_Sensitivity_Analysis.RData")

#### Define variable conditions ####
mig_rate_conditions <- unique(1/c(seq(1, 10, length.out = 20),seq(10, 40, length.out = 10)))
pop_size_conditions <- 10000 # I confirmed that the impact of N*avg_prob is compound and therefore only one needs to be varied
avg_prob_conditions <- 0.00001 * 10^seq(0, 2, length.out = 20) # One order of magnitude lower and higher than 0.0001 (or 1/10,000)
R0_conditions <- seq(0.75, 3, length.out = 10)

df_dd_sens <- data.frame(mig = rep(mig_rate_conditions, each = length(pop_size_conditions)*length(avg_prob_conditions)*length(R0_conditions)),
                         pop_size = rep(rep(pop_size_conditions, each = length(avg_prob_conditions)*length(R0_conditions)), times = length(mig_rate_conditions)),
                         avg_prob = rep(rep(avg_prob_conditions, each = length(R0_conditions), times = length(mig_rate_conditions)*length(pop_size_conditions))),
                         R0 = rep(R0_conditions, times = length(mig_rate_conditions)*length(pop_size_conditions)*length(avg_prob_conditions)))

df_dd_sens$prob_outbreak_novax <- NA
df_dd_sens$prob_outbreak_dukoral <- NA

nrow(df_dd_sens)

#### Define static conditions
horizon = 4 # years
time_step = 1 # day
seasonal_amp = 0 # 0 if no seasonality. 1 if doubles
vaccine_choices = c("None", "Dukoral")
outbreak_size = 10 # min cases for outbreak definition

#### Loop through conditions
for (i in 1:nrow(df_dd_sens)){
  df_dd_sens[i, "prob_outbreak_novax"] <- cum_prob_outbreak(years = horizon, mig_rate = df_dd_sens[i, "mig"], pop_size = df_dd_sens[i, "pop_size"], time_step = time_step, avg_prob = df_dd_sens[i, "avg_prob"], seasonal_amp = seasonal_amp, R = df_dd_sens[i, "R0"], outbreak_size = outbreak_size, vaccine_choice = "None")
  df_dd_sens[i, "prob_outbreak_dukoral"] <- cum_prob_outbreak(years = horizon, mig_rate = df_dd_sens[i, "mig"], pop_size = df_dd_sens[i, "pop_size"], time_step = time_step, avg_prob = df_dd_sens[i, "avg_prob"], seasonal_amp = seasonal_amp, R = df_dd_sens[i, "R0"], outbreak_size = outbreak_size, vaccine_choice = "Dukoral")
  cat(".")
  if (i %% 100 == 0){cat("|\n")}
  }

df_dd_sens$dukoral_impact <- df_dd_sens$prob_outbreak_novax - df_dd_sens$prob_outbreak_dukoral

#### Find optimal strategy for each setting ####
settings <-  data.frame(pop_size = rep(pop_size_conditions, each = length(avg_prob_conditions)*length(R0_conditions)),
                           avg_prob = rep(rep(avg_prob_conditions, each = length(R0_conditions), times = length(pop_size_conditions))),
                           R0 = rep(R0_conditions, times = length(pop_size_conditions)*length(avg_prob_conditions)))
settings$optimal_impact <- NA
settings$optimal_mig <- NA

for (i in 1:nrow(settings)){
  pop_size <- settings[i, "pop_size"]
  avg_prob <- settings[i, "avg_prob"]
  R0 <- settings[i, "R0"]
  
  settings[i, "optimal_impact"] <- max(df_dd_sens[df_dd_sens$pop_size == pop_size & df_dd_sens$R0 == R0 & df_dd_sens$avg_prob == avg_prob,]$dukoral_impact)
  settings[i, "optimal_mig"] <- df_dd_sens[df_dd_sens$pop_size == pop_size & df_dd_sens$R0 == R0 & df_dd_sens$avg_prob == avg_prob & df_dd_sens$dukoral_impact == settings[i, "optimal_impact"],]$mig[1]
  
}
settings

settings$mig_time <- 1/settings$optimal_mig
settings$pop_prob <- (settings$pop_size*settings$avg_prob)


hist(1/settings$optimal_mig, breaks = 1/mig_rate_conditions)

summary(factor(1/settings$optimal_mig))
sort(summary(factor(1/settings$optimal_mig)), decreasing = TRUE)[1] / nrow(settings)

plot(x = settings[settings$pop_prob == 1 & settings$R0 == unique(settings$R0)[3],"pop_size"], y = settings[settings$pop_prob == 1 & settings$R0 == unique(settings$R0)[3],"optimal_mig"])

#### Heatmap plot ####
ggplot(settings, aes(x = R0, y = avg_prob, fill = log(mig_time))) +
  geom_tile() +
  # geom_text(aes(label = round(mig_time,0)), col = "black") +
  scale_y_log10(expand = c(0,0), name = expression(paste("Fraction of migrants infected ",(mu)))) +
  scale_x_continuous(expand = c(0,0), name = expression(paste(R[0]))) +
  scale_fill_gradient2(mid="white", high = "blue", low="red", midpoint = 2) +
  guides(fill = FALSE) +
  coord_fixed() +

ggsave(file = "figures/Figure_DD_Sensitivity_Analysis.pdf", width = 6, height = 5, units = "in")


#### Save workspace ####
save.image(file = "src/Figure_DD_Sensitivity_Analysis.RData")
