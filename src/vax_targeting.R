#### Influence of mobility on optimal vaccine targeting ####

# setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
# source("src/prob_outbreak_fcn.R")
# source("src/calculate_VE.R")

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
    VE_t <- Create_VE(timesteps_per_month = 30, VE_shape = vaccine_choice, bound=TRUE, max_V_months = 12*years)
    
    Xt <- c(1 - (VE_t*exp(-1*(mig_rate/365)*(1:(length(VE_t))))), rep(1, times = (365*years - length(VE_t))))
  }
  
  df <- data.frame(days = 1:(365*years), seasonal_time = rep(seq_len(365), years), Xt = Xt, Nt = NA, Ct = NA, Rt = NA, prob_outbreak_t_per_inf = NA, prob_outbreak_t = NA)

  for (day in 1:nrow(df)){
    df[day, "Nt"] <- expected_migrants(mig_rate = mig_rate, pop_size = pop_size, time_step = time_step)
    df[day, "Ct"] <- df[day, "Nt"] * prob_migrant_infected(avg_prob = avg_prob, seasonal_amp = seasonal_amp, seasonal_time = df[day, "seasonal_time"])
    df[day, "Rt"] <- R * (1+seasonal_amp*cos(2*pi*df[day, "seasonal_time"]/365)) * df[day, "Xt"]
    df[day, "prob_outbreak_t_per_inf"] <- prob_outbreak_fcn(R = df[day,"Rt"], outbreak_size = outbreak_size)
    df[day, "prob_outbreak_t"] <- 1-(1-prob_outbreak_fcn(R = df[day,"Rt"], outbreak_size = outbreak_size))^df[day, "Ct"]
  }

  return(1-prod((1-df[, "prob_outbreak_t"])))
}

#### Test to make sure this method works for fractional cases on each day ####
df_test <- data.frame(day = 1:6, Rt =1.1, Ct = c(1,1,0.5, 0.5, 0.5, 0.5), prob_outbreak_t = NA)
for (i in seq_len(nrow(df_test))){df_test[i, "prob_outbreak_t"] <- 1-(1-prob_outbreak_fcn(R = df_test[i,"Rt"], outbreak_size = outbreak_size))^df_test[i, "Ct"]}
df_test
1-prod(1-df_test[1:2, "prob_outbreak_t"])
1-prod(1-df_test[3:6, "prob_outbreak_t"])
0