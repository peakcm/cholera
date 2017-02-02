#### Graphically show that all-or-none and leaky vaccines lead to the same X(t) profiles ####
# Show VE and V on a plot. For all-or-none, V wanes but not VE. For leaky, VE wanes but not V.
# Show X(t) for both below it

setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")
library(ggplot2)
library(reshape2)

#### Vaccine time series ####
VE_0 <- 0.75 #Initial VE
time_points <- 365*1 # 10 years

V_all_or_none <- exp(-5*(0:time_points)/time_points) * VE_0
VE_leaky <- V_all_or_none


V_leaky <- rep(1, time_points+1)
VE_all_or_none <- V_leaky

df <- data.frame(time = rep(0:time_points, 2), 
                 action = rep(c("all_or_none", "leaky"),each = time_points+1),
                 V = c(V_all_or_none, V_leaky),
                 VE = c(VE_all_or_none, VE_leaky))

df$action <- factor(df$action, labels = c("Failure in Take\n(All or Nothing)", "Failure in Degree\n(Leaky)"))

df$X <- (1-df$V) + df$V*(1-df$VE)  # 1 compartmental model.

df_melt <- melt(df, id.vars = c("time", "action"))
df_melt$variable <- factor(df_melt$variable, labels = c("V(t)", "VE(t)", "X(t)"))

#### Plot V and VE time series ####
ggplot(df_melt, aes(x = time/36.5, y = value, group = variable, lty = variable)) +
  facet_grid(.~action) +
  theme_bw() +
  geom_line() +
  scale_y_continuous(name = "") +
  scale_x_continuous(name = "Years") +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), name = "") +
  theme(text = element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8)) +
  

ggsave(file = "figures/Figure_HH.pdf", width = 5, height = 3, units = "in")

