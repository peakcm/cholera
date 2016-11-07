#### Influence of mobility on optimal vaccine targeting ####

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
source("src/prob_outbreak_fcn.R")
library(ggplot2)

# Sum from 0 to 5 years by day. for (A)*(B)*(C)
# (A) Expected number of immigrants at time t
  # Directly relate to the migration rate and assumed to be constant over time
# (B) Probability an immigrant is infected at time t
  # Can vary over time if immigrant is coming from a region with seasonality as well. Assume predominantly local travel, so seasonality is synchronized
# (C) Probability an infected immigrant will generate an outbreak at time t.
  # Depends on R_0(t) and X(t)

#### Function (A) ####
expected_migrants <- function(mig_rate, pop_size, time_step){
  # mig_rate is per year
  # time is days
  pop_size* (1-exp(-1*mig_rate*time_step/365))
}

#### Function (B) ####
prob_migrant_infected <- function(avg_prob, seasonal_amp=0.5, seasonal_time=0){
  # avg_prob is the yearly averge probability a migrant is infected
  # seasonal_amp is the amplitude of the seasonality, which ranges from 0 (no seasonality) to 1
  # seasonal_time is the current time in the season, assuming annual cosine seasonality
  avg_prob * (1+seasonal_amp * cos(2*pi*seasonal_time/365))
}

#### Function (C) ####
# See "prob_outbreak_fcn"

#### Calculate prob_outbreak ####
cum_prob_outbreak <- function(years, mig_rate, pop_size, time_step, avg_prob, seasonal_amp, R, outbreak_size, vaccine_choice){
  if (vaccine_choice == "None"){
    Xt <- rep(1, 365*years)
  } else {
    VE_t <- Create_VE(timesteps_per_month = 30, VE_shape = vaccine_choice, bound=TRUE, 12*years)
    
    Xt <- c(1 - (VE_t*exp(-1*(mig_rate/365)*(1:(length(VE_t))))), rep(1, times = (365*years - length(VE_t))))
  }
  
  df <- data.frame(days = 1:(365*years), seasonal_time = rep(seq_len(365), years), Xt = Xt, Nt = NA, Ct = NA, Rt = NA, prob_outbreak_t = NA)

  for (day in 1:nrow(df)){
    df[day, "Nt"] <- expected_migrants(mig_rate = mig_rate, pop_size = pop_size, time_step = time_step)
    df[day, "Ct"] <- df[day, "Nt"] * prob_migrant_infected(avg_prob = avg_prob, seasonal_amp = seasonal_amp, seasonal_time = df[day, "seasonal_time"])
    df[day, "Rt"] <- R * (1+seasonal_amp*cos(2*pi*df[day, "seasonal_time"]/365)) * df[day, "Xt"]
    df[day, "prob_outbreak_t"] <- 1-(1-prob_outbreak_fcn(R = df[day,"Rt"], outbreak_size = outbreak_size))^df[day, "Ct"]
  }

  return(1-prod((1-df[, "prob_outbreak_t"])))
}

#### Practice ####
years = 5
mig_rate = 1/200
pop_size = 1000
time_step = 1 # day
avg_prob = 0.01 # yearly average prob that immigrant is infected
seasonal_amp = 0.1 # 0 if no seasonality. 1 if doubles
R = 1.5 # yearly average
outbreak_size = 10 # min cases for outbreak definition
vaccine_choice = "None"
# vaccine_choice = "Shanchol"
# vaccine_choice = "Dukoral"

cum_prob_outbreak(years = years, mig_rate = mig_rate, pop_size = pop_size, time_step = time_step, avg_prob = avg_prob, seasonal_amp = seasonal_amp, R = R, outbreak_size = outbreak_size, vaccine_choice = vaccine_choice)

#### Loop ####
vaccine_choices = c("None", "Shanchol")
mig_rate_choices = seq(0, 1/5, length.out = 40)

df <- data.frame(vaccine_choice = rep(vaccine_choices, each = length(mig_rate_choices)), mig_rate = rep(mig_rate_choices, length(vaccine_choices)), prob_outbreak = NA)

for (i in 1:nrow(df)){
  df[i, "prob_outbreak"] <- cum_prob_outbreak(years = years, mig_rate = df[i, "mig_rate"], pop_size = pop_size, time_step = time_step, avg_prob = avg_prob, seasonal_amp = seasonal_amp, R = R, outbreak_size = outbreak_size, vaccine_choice = df[i, "vaccine_choice"])
  cat(".")
}

df_diff <- data.frame(mig_rate = mig_rate_choices, diff = NA)
for (i in 1:nrow(df_diff)){
  df_diff[i, "diff"] <- df[df$vaccine_choice == "None" & df$mig_rate == df_diff[i,"mig_rate"], "prob_outbreak"] - df[df$vaccine_choice == "Shanchol" & df$mig_rate == df_diff[i,"mig_rate"], "prob_outbreak"]
}

ggplot() + geom_bar(data = df_diff, aes(x = mig_rate, y = diff), stat = "identity", color = "lightgrey", fill = "grey") + geom_line(data = df, aes(x = mig_rate, y = prob_outbreak, color = vaccine_choice)) + theme_bw() + ylab("Probability of an Outbreak") + xlab("Migration Rate (per year)") + scale_color_discrete(name = "Vaccine Status")
