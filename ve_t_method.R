#### Practice proportional hazards testing of schoenfeld residuals ####
library(survival)

#### Create Fictional OCV data ####
years = 5
samples = 100000
VE = 0.4  # Best if VE < 0.5 because I'm assuming linear decay down to the background event probability
event_prob_untreated <- 0.5
event_prob_treated <- event_prob_untreated*(1-VE)
fake_ocv <- data.frame(ID = 1:samples,
                       treat = rep(c(0,1), each = samples/2),
                       time = round(runif(n=samples, min = 1, max = 365*years)),
                       event = c(rbinom(n = samples/2, size = 1, prob = event_prob_untreated), rbinom(n = samples/2, size = 1, prob = event_prob_treated))) # Event times occur randomly, but some are censoring and some are disease

# Add a time-varying effect that is strong at first and then gets weaker
fake_ocv <- fake_ocv[order(fake_ocv$treat, fake_ocv$time),] # Sort by treatment status, then time
fake_ocv$event_tv <- fake_ocv$event
fake_ocv[(samples/2+1):samples, "event_tv"] <- sapply(seq(min(1, VE*2), 0, length.out = samples/2), function(x) rbinom(n=1, size=1, prob = event_prob_untreated*(1-x)))

#### Confirm behavior of Fictional Data ####
# Calculate VE via 1-RR
VE
1 - mean(fake_ocv[fake_ocv$treat == 1, "event"]) / mean(fake_ocv[fake_ocv$treat == 0, "event"])
1 - mean(fake_ocv[fake_ocv$treat == 1, "event_tv"]) / mean(fake_ocv[fake_ocv$treat == 0, "event_tv"]) # Biased downwards b/c sample size larger earlier

# Calculate Hazard Ratio over all time
event_prob_treated/event_prob_untreated
summary(coxph(Surv(time,event) ~ treat, data=fake_ocv))[c("coefficients", "conf.int")]
summary(coxph(Surv(time,event_tv) ~ treat, data=fake_ocv))[c("coefficients", "conf.int")] 

# Compare Hazard Ratio in first half to second half
summary(coxph(Surv(time,event_tv) ~ treat, data=fake_ocv[fake_ocv$time < max(fake_ocv$time)/2,]))[c("coefficients", "conf.int")]
summary(coxph(Surv(time,event_tv) ~ treat, data=fake_ocv[fake_ocv$time > max(fake_ocv$time)/2,]))[c("coefficients", "conf.int")]

# Plot survival curves 
plot(survfit(formula = Surv(time, event)~ treat, data = fake_ocv, conf.type="none"), xlab="Time", ylab="Survival Probability", main = "Constant VE effect")
plot(survfit(formula = Surv(time, event_tv)~ treat, data = fake_ocv, conf.type="none"), xlab="Time", ylab="Survival Probability", main = "Time-waning VE effect")

#### Including Time Dependent Covariates in the Cox Model - static VE ####
time.dep <- coxph( Surv(time, event)~treat,data = fake_ocv,
                   method="breslow", na.action=na.exclude)
(time.dep.zph <- cox.zph(time.dep))
# plot(time.dep.zph, main = "Constant VE effect")

time.dep.zph.df <- data.frame(x = as.numeric(time.dep.zph$x), y = as.numeric(time.dep.zph$y))

# ggplot(time.dep.zph.df, aes(x=x, y=y)) + geom_point() + stat_smooth(method="lm", formula = y ~ splines::bs(x,4))

time.dep.zph.df.spline <- lm(time.dep.zph.df$y ~ splines::bs(time.dep.zph.df$x,4))
plot(time.dep.zph.df$x*years, 1-exp(time.dep.zph.df.spline$fitted.values), ylim = c(0,1), type = "l", ylab = "VE(t)", xlab = "Years")
abline(h = VE, lty = 3)

#### Including Time Dependent Covariates in the Cox Model - waning VE ####
time.dep.tv <- coxph( Surv(time, event_tv)~treat,data = fake_ocv,
                   method="breslow", na.action=na.exclude)
(time.dep.zph.tv <- cox.zph(time.dep.tv))
# plot(time.dep.zph.tv, main = "Time-waning VE effect", df=4) # Fits a spline (defualt df=4)

time.dep.zph.tv.df <- data.frame(x = as.numeric(time.dep.zph.tv$x), y = as.numeric(time.dep.zph.tv$y))
# ggplot(time.dep.zph.tv.df, aes(x=x, y=y)) + geom_point() + stat_smooth(method="lm", formula = y ~ splines::bs(x,4))

time.dep.zph.tv.df.spline <- lm(time.dep.zph.tv.df$y ~ splines::bs(time.dep.zph.tv.df$x,4))
plot(time.dep.zph.tv.df$x*years, 1-exp(time.dep.zph.tv.df.spline$fitted.values), ylim = c(-.1,1), type = "l", ylab = "VE(t)", xlab = "Years")
abline(h = VE, lty = 3)

