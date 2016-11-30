#### Figure DD ####
setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
source("src/vax_targeting.R")
library(ggplot2)

#### Load workspace ####
load(file = "src/Figure_DD.RData")

#### Notes ####

# Higher R, lower optimal Migration Rate?
# Higher seasonal amplitude, lower optimal Migration Rate?
# Pop size and avg_prob are inversely associated?
# Higher min outbreak size increases the difference between vax and non-vax, but does not shift optimal migration rate?

#### Create data frame for migration rate call-outs ####
# Calculate the migration rate to satisfy the turnover in Bentiu. Between Feb 2016 and Oct 2016, an average of 2,000 entries and exits occured while the population remained static at around 104,000 people
N_bentiu <- 104000
(Bentiu_rate <- 1/(365*(log((N_bentiu-2000)/N_bentiu, base = exp(1))/30.5*-1)))

call_outs <- data.frame(labels = c("Dhaka\n[Qadri et al 2015]", "Calcutta\n[Sur et al 2011]", "Bentiu\n[IOM 2016]"),
                        x = c(1/2, 1/20, 1/Bentiu_rate),
                        y = c(1.075, 1.075, 1.075))

#### Practice ####
years = 5
mig_rate = 1/20
pop_size = 1000
time_step = 1 # day
avg_prob = 1/pop_size # yearly average prob that immigrant is infected
seasonal_amp = 0 # 0 if no seasonality. 1 if doubles
R = 1.5 # yearly average
outbreak_size = 10 # min cases for outbreak definition
vaccine_choice = "None"
# vaccine_choice = "Shanchol"
# vaccine_choice = "Dukoral"

cum_prob_outbreak(years = years, mig_rate = mig_rate, pop_size = pop_size, time_step = time_step, avg_prob = avg_prob, seasonal_amp = seasonal_amp, R = R, outbreak_size = outbreak_size, vaccine_choice = vaccine_choice)

#### Single Loop through migration rates ####
vaccine_choices = c("None", "Shanchol")
mig_rate_choices = seq(0, 1/1, length.out = 40)

df_single <- data.frame(vaccine_choice = rep(vaccine_choices, each = length(mig_rate_choices)), mig_rate = rep(mig_rate_choices, length(vaccine_choices)), prob_outbreak = NA)

for (i in 1:nrow(df_single)){
  df_single[i, "prob_outbreak"] <- cum_prob_outbreak(years = years, mig_rate = df_single[i, "mig_rate"], pop_size = pop_size, time_step = time_step, avg_prob = avg_prob, seasonal_amp = seasonal_amp, R = R, outbreak_size = outbreak_size, vaccine_choice = df_single[i, "vaccine_choice"])
  cat(".")
}

df_single_diff <- data.frame(mig_rate = mig_rate_choices, diff = NA)
for (i in 1:nrow(df_single_diff)){
  df_single_diff[i, "diff"] <- df_single[df_single$vaccine_choice == "None" & df_single$mig_rate == df_single_diff[i,"mig_rate"], "prob_outbreak"] - df_single[df_single$vaccine_choice == "Shanchol" & df_single$mig_rate == df_single_diff[i,"mig_rate"], "prob_outbreak"]
}

#### Plot single year-horizon ####
ggplot() + 
  # call-outs
  geom_vline(data = call_outs, aes(xintercept = x), col = "grey", lty = "dashed") + 
  geom_text(data = call_outs, aes(x = x, y = y, label = labels), color = "grey", angle = 90, nudge_x = 0, size = 2) +
  # bars
  geom_bar(data = df_single_diff, aes(x = mig_rate, y = diff, fill = "Vaccine\nImpact"), stat = "identity") +
  # lines
  geom_line(data = df_single, aes(x = mig_rate, y = prob_outbreak, color = vaccine_choice), alpha = 0.5) + 
  # Formatting
  theme_bw() + ylab("5-Year Probability of an Outbreak\nInitiated by an Imported Case") + scale_color_discrete(name = "Vaccine Status") +
  scale_x_continuous(breaks = c(0.05, 0.20, 0.333, 0.5, 1/1.5, 1), labels = c("1/20", "1/5", "1/3", "1/2", "1/1.5", "1/1"), name = "Migration Rate (per year)") + theme(text = element_text(size = 6), legend.text=element_text(size=4), legend.title=element_text(size=4), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) + theme(legend.position = c(.8, .4)) + scale_fill_manual(name = element_blank(), values = "darkgrey") + guides(color = guide_legend(order = 1)) + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.50", "0.75", "1"), limits = c(0,1.2)) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

ggsave(file = "figures/Figure_DD_single.pdf", width = 3, height = 3, units = "in")

#### Multiple Loops through migration rates ####
horizons <- c(2, 3, 4, 5, 6)

df <- NA
df_diff <- NA
for (years in horizons){
  vaccine_choices = c("None", "Shanchol")
  mig_rate_choices = seq(0, 1/1, length.out = 40)
  
  df_temp <- data.frame(vaccine_choice = rep(vaccine_choices, each = length(mig_rate_choices)), mig_rate = rep(mig_rate_choices, length(vaccine_choices)), prob_outbreak = NA)
  
  for (i in 1:nrow(df_temp)){
    df_temp[i, "prob_outbreak"] <- cum_prob_outbreak(years = years, mig_rate = df_temp[i, "mig_rate"], pop_size = pop_size, time_step = time_step, avg_prob = avg_prob, seasonal_amp = seasonal_amp, R = R, outbreak_size = outbreak_size, vaccine_choice = df_temp[i, "vaccine_choice"])
    cat(".")
  }
  
  df_temp_diff <- data.frame(mig_rate = mig_rate_choices, diff = NA)
  for (i in 1:nrow(df_temp_diff)){
    df_temp_diff[i, "diff"] <- df_temp[df_temp$vaccine_choice == "None" & df_temp$mig_rate == df_temp_diff[i,"mig_rate"], "prob_outbreak"] - df_temp[df_temp$vaccine_choice == "Shanchol" & df_temp$mig_rate == df_temp_diff[i,"mig_rate"], "prob_outbreak"]
  }
  
  df_temp$years <- years
  df_temp_diff$years <- years
  
  df <- rbind(df, df_temp)
  df_diff <- rbind(df_diff, df_temp_diff)
  cat(years, "\n")
}
df_diff <- df_diff[is.na(df_diff$years) == 0,]
df_diff$years <- as.factor(df_diff$years)

#### Dataset with just optimal migration rates ####
df_diff_optimal <- data.frame(years = unique(df_diff$years), optimal_mig_rate = NA, diff = NA)
for (i in 1:nrow(df_diff_optimal)){
  df_diff_optimal[i, c("diff")] <- max(df_diff[df_diff$years == df_diff_optimal[i,"years"], "diff"])
  df_diff_optimal[i, c("optimal_mig_rate")] <- df_diff[df_diff$years == df_diff_optimal[i,"years"] & df_diff$diff == df_diff_optimal[i,"diff"], "mig_rate"]
}

#### Plot multiple year-horizons ####
ggplot() + 
  # call-outs
  geom_vline(data = call_outs, aes(xintercept = x), col = "grey", lty = "dashed") + 
  geom_text(data = call_outs, aes(x = x, y = y-.6, label = labels), color = "grey", angle = 90, nudge_x = 0, size = 1.5) +
  geom_point(data = df_diff_optimal, aes(x = optimal_mig_rate, y = diff, color = years), pch = 17) +
  # bars
  geom_line(data = df_diff, aes(x = mig_rate, y = diff, group = years, color = years)) +
  # Formatting
  theme_bw() + ylab("Vaccine Impact\n(decrease in outbreak probability)") + 
  scale_color_brewer(type = "seq", palette = 8, name = "Time Horizon\n(years)") +
  scale_x_continuous(breaks = c(0.05, 0.20, 0.333, 0.5, 1/1.5, 1), labels = c("1/20", "1/5", "1/3", "1/2", "1/1.5", "1/1"), name = "Migration Rate (per year)") + theme(text = element_text(size = 6), legend.text=element_text(size=4), legend.title=element_text(size=4), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8)) + scale_y_continuous(breaks = c(0, 0.25, 0.5), labels = c("0", "0.25", "0.50"), limits = c(0,0.6)) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

ggsave(file = "figures/Figure_DD.pdf", width = 4, height = 2, units = "in")

#### Save workspace ####
save.image(file = "src/Figure_DD.RData")
