#### Header ####
# Analytic Solution
# Corey Peak 12/9/15

#### Load Libraries ####
library(ggplot2)

#### Define Initial Parms ####
parms_init = c(u = 1/(60*365),   # Death rate (1/days)
          b = 1/(60*365),   # Birth rate (1/days)
          B0 = 0.5,         # Transmission parameter for a non-vaccinated susceptible person contacting a non-vaccinated infectious person
          k1 = 0.0,         # Reduction in infectiousness of a once-vaccinated person
          k2 = 0.0,         # Reduction in infectiousness of a twice-vaccinated person
          vac1 = 0.5,       # Personal efficacy of one dose of vaccine
          vac2 = 0.85,      # Personal efficacy of two doses of vaccine
          r = 1/3,          # Recovery rate (1/days)
          e = 1/2,          # Incubation (1/days)
          ntot = 25000,     # Population size
          v1_day = 50,      # Day of first vaccine dose (If you don't want vaccination, then set this to 100 or whatever your t_final is)
          v2_day = 64,      # Day of second vaccine dose
          v1_count = 1000,  # number of vaccines intended for first dose
          v2_count = 1000,  # number of vaccines intended for second dose
          m = .10)          # annual turnover rate due to migration (same for immigration and emigration so constant pop size)

#### Define Initial States ####
states_init <- c(S_0=1000,
           S_1=0,
           S_2=0,
           E_0=0,
           E_1=0,
           E_2=0,
           I_0=0,
           I_1=0,
           I_2=0,
           R_0=0,
           R_1=0,
           R_2=0)

#### Define Vaccine Waning Function ####
VacWane <- function(parms, day, shape = "exponential", decay_parm){  
  if (shape == "exponential"){
    parms["vac1"] <- parms["vac1"] * (1 - decay_parm )
    parms["vac2"] <- parms["vac2"] * (1 - decay_parm )
  }
  return(parms)
}

#### Define Migrate function ####
Migrate <- function(states, parms){
  N <- sum(states)
  m <- parms["m"]
  S_0_in <- (m/365)*N #account for people coming in susceptible
  
  states <- states*(1 - m/365) #account for people leaving in any state proportionally
  states["S_0"] <- states["S_0"] + S_0_in 

  return(states)
}

#### Define Vaccinate function ####
Vaccinate <- function(states, parms){
  if (parms["v1_count"] > 0){
    if (parms["v1_count"] <= states["S_0"]){
      states["S_0"] <- states["S_0"] - parms["v1_count"]
      states["S_1"] <- parms["v1_count"]
      
    } else {cat("error v1_count")} # Haven't added all the possibilities yet
  }
  if (parms["v2_count"] > 0){
    if (parms["v2_count"] <= states["S_1"]){
      states["S_1"] <- states["S_1"] - parms["v2_count"]
      states["S_2"] <- parms["v2_count"]
    } else {cat("error v2_count")} # Haven't added all the possibilities yet
  }
  return(states)
}

#### Solve ####
Solve <- function(states, parms, years, decay_parm){
  output <- data.frame(matrix(rep(NA, length(states)*years*365), nrow = years*365))
  names(output) <- names(states)
  output$X <- NA  # X is the functional proportion suceptible
  
  states <- Vaccinate(states, parms)
  
  X <- ( states["S_0"] + (1-parms["vac1"])*states["S_1"] + (1-parms["vac2"])*states["S_2"] + sum(states[c("R_0", "R_1", "R_2")]) ) / sum(states)
  output[1,] <- c(states, X)
  
  for (day in 2:(365*years)){
    parms <- VacWane(parms, day, shape = "exponential", decay_parm)
    states <- Migrate(states, parms)
    
    X <- ( states["S_0"] + (1-parms["vac1"])*states["S_1"] + (1-parms["vac2"])*states["S_2"] + sum(states[c("R_0", "R_1", "R_2")]) ) / sum(states)
    output[day,] <- c(states, X)
  }
  
  return(output)
}

#### Test Solve function ####
years = 10
decay_parm = .00055 # Moore 2015 JRoyalSoc. 1/(365*5)
parms_init["m"] = 0.0

output <- Solve(states = states_init, parms_init, years, decay_parm)

#### Plot Solve Function Compartments ####
layout(mat=matrix(c(1,2,3), nrow = 3))
plot(output$S_0, main = "never vaccinated", xlab = "days", ylab = "count")
plot(output$S_2, main = "twice-vaccinated", xlab = "days", ylab = "count")
plot(output$X, main = "fraction susceptible", xlab = "days", ylim = c(0,1), ylab = "fraction")

#### Repeat for different levels of migration ####
layout(c(1))
m <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
dur <- rep(NA, length(m))

data <- data.frame(cbind(m, dur))

years = 10
decay_parm = .00055 # Moore 2015 JRoyalSoc. 1/(365*5)
R <- 2

for (i in 1:nrow(data)){
  parms_init["m"] <- data[i, "m"]
  output <- Solve(states = states_init, parms_init, years, decay_parm)
  herd_immunity_years <- which(1 > (R* output$X) )
  if (length(herd_immunity_years) == 0){
    data[i, "dur"] <- 0
  } else {
    data[i, "dur"] <- max(herd_immunity_years)
  }
  data$R <- R
  data$decay_parm <- decay_parm
}

ggplot(data, aes(x=m, y=dur/365)) +
  theme_bw()+
  geom_point() +
  xlab("Migration Rate (Turnover Fraction per Year)") +
  ylab("Years of Herd Immunity")

#### Repeat for different levels of migration AND R ####
layout(c(1))
m <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
R <- c(2, 3, 4, 5)

data <- data.frame(cbind(m=rep(m, length(R)), R=rep(R,each=length(m))))
data$dur <- NA

years = 10
decay_parm = .00055 # Moore 2015 JRoyalSoc. 1/(365*5)

for (i in 1:nrow(data)){
  parms_init["m"] <- data[i, "m"]
  R <- data[i, "R"]
  output <- Solve(states = states_init, parms_init, years, decay_parm)
  herd_immunity_years <- which(1 > (R* output$X) )
  if (length(herd_immunity_years) == 0){
    data[i, "dur"] <- 0
  } else {
    data[i, "dur"] <- max(herd_immunity_years)
  }
  data$decay_parm <- decay_parm
}

data$m <- factor(data$m)
data$R <- factor(data$R)

ggplot(data, aes(x=m, y=dur/365, group=R, col=R)) +
  theme_bw()+
  geom_line(size=2) +
  xlab("Migration Rate (Turnover Fraction per Year)") +
  ylab("Years of Herd Immunity")
