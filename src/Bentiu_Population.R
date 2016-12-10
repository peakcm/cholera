# Bentiu case study #
setwd("/Users/peakcm/Dropbox/Cholera Amanda/cholera_waning")

#### Load Libraries ####
library(ggplot2)

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

#### Plot ####
ggplot(df) + 
  # geom_rect(data = df_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white", color = "black") +
  # geom_bar(data = df_cases, aes(x = Date, y = Cases_weekly*500), stat = "identity", width = 2) + 
  geom_bar(aes(x = Date, y = Cholera.cases*1000), stat = "identity", fill = "lightgrey") + 
  geom_text(aes(x = Date, y = (Cholera.cases*1000 + 3000), label = Cholera.cases), size = 2) + 
  geom_ribbon(data = df[is.na(df$influx)==0,], aes(x = Date, ymin = pop, ymax = influx), fill = "lightgreen") +
  geom_ribbon(data = df[is.na(df$outflux)==0,], aes(x = Date, ymin = outflux, ymax = pop), fill = "pink") +
  geom_line(aes(x = Date, y = pop), color = "darkblue") +
  geom_bar(data = df[is.na(df$Vaccine.Doses)==0,], aes(x = Date, y = Vaccine.Doses/2), stat = "identity", width = 6, color = "forestgreen") +
  scale_y_continuous(breaks = c(0, 5e4, 1e5, 1.5e5), labels = c("0", "50", "100", "150"), name = "Population (thousands)") +
  scale_x_date()
  theme_bw() +
  theme(text = element_text(size = 6), legend.text=element_text(size=4), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
  ggtitle("Bentiu PoC Camp")

ggsave(file = "figures/Figure_GG.pdf", width = 3, height = 3, units = "in")
