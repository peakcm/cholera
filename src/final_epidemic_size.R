# Adapted from [Black and Ross 2015 J Theoretical Biology]
# http://www.sciencedirect.com/science/article/pii/S0022519314006882

final_epidemic_size <- function(N, R0){
  
  beta = R0/(N-1)
  gamma = 1
  
  q = rep(0, N+1)
  q[2] = 1
  
  for (Z2 in seq(0,N)){
    for (Z1 in seq(Z2+1, N-1)){
      p1 = 1/(1 + gamma/(beta*(N-Z1)))
      
      q[Z1+2] = q[Z1+2] + q[Z1+1]*p1
      q[Z1+1] = q[Z1+1]*(1-p1)
    }
  }
  return(q)
}

# Test function behavior
N = 500
R = 1.25
q <- final_epidemic_size(N, R)
cat("The probability of <10 cases is",sum(q[1:10]))
cat("The probability of 10 to", round(N/3,0),"cases is", sum(q[11:(N/3)]))
cat("The probability of >", round(N/3,0),"cases is", sum(q[(N/3):N]))
plot(q*N, type = "h")
plot(q*N, type = "h", xlim = c(0,30))
q_hist <- hist(q*N*seq(0, length(q)-1), breaks = 10)
# 
# # Probability of a large outbreak (>1% AR) given you make it past "cut_point" cases
# cut_point <- 10
# sum(q[(N/100):N]) / sum(q[(cut_point+1):N])
# When N=20,000 and R=1.5, cut_point 10 is 89% predictive. cut_point 25 is 97% predictive
# When N=10,000 and R=1.5, cut_point 10 is 88% predictive. cut_point 25 is 97% predictive
# When N=10,000 and R=1.4, cut_point 10 is 82% predictive. cut_point 25 is 95% predictive
# When N=10,000 and R=1.1, cut_point 10 is 47% predictive. cut_point 25 is 68% predictive. cut_point 80 is 95% predictive....but not very interesting 

#### Sensitivity analysis of the use of 10 cases as outbreak threshold ####
library(data.table)
library(ggplot2)

N_conditions <- c(1000, 5000) # Plot is insensitive to N, but Sensitivity measurements gets accurate with 5,000
R0_conditions <- c(0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3)

df_sims <- data.table(N = rep(N_conditions, each = length(R0_conditions)), R0 = rep(R0_conditions, times = length(N_conditions)), q = NA)

for (row in 1:nrow(df_sims)){
  N <- df_sims$N[row]
  R <- df_sims$R[row]
  q <- final_epidemic_size(N, R)[2:N]
  df_sims$q[row] <- list(q)
  cat("\n", row)
}

plot((unlist(df_sims$q[4])), type = "h")

pdf("figures/final_epidemic_size.pdf", width=5, height=4)
plot((unlist(df_sims$q[1])), type = "l", xlim = c(0, 30), col = "red", ylab = "Density", xlab = "Cases", xaxt="n")
axis(1,at = 0:30)
lines((unlist(df_sims$q[4])), xlim = c(0, 30), col = "yellow")
lines((unlist(df_sims$q[7])), xlim = c(0, 30), col = "green")
lines((unlist(df_sims$q[10])), xlim = c(0, 30), col = "blue")
lines((unlist(df_sims$q[13])), xlim = c(0, 30), col = "black")
abline(v = 10, lty = "dashed")
dev.off()

pdf("figures/final_epidemic_size_subplot.pdf", width=5, height=4)
plot((unlist(df_sims$q[1])), type = "l", ylim = c(0, 0.05), col = "red", ylab = "Density", xlab = "Cases")
lines((unlist(df_sims$q[4])),  ylim = c(0, 0.05), col = "yellow")
lines((unlist(df_sims$q[7])),  ylim = c(0, 0.05),col = "green")
lines((unlist(df_sims$q[10])),  ylim = c(0, 0.05), col = "blue")
lines((unlist(df_sims$q[13])),  ylim = c(0, 0.05), col = "black")
dev.off()

plot(cumsum(unlist(df_sims$q[3])), type = "l")

df_sims$cut_10 <- NA
for (row in 1:nrow(df_sims)){
  cut_point = 10
  df_sims$cut_10[row] <- sum(unlist(df_sims$q[row])[100:(df_sims$N[row]-1)]) / sum(unlist(df_sims$q[row])[(cut_point+1):(df_sims$N[row]-1)])
}
df_sims$Population <- factor(df_sims$N)

ggplot(df_sims, aes(x = R0, y = cut_10, color = Population)) + geom_line() +
  geom_point() +
  scale_y_continuous(name = "Probability of 100+ cases\ngiven 10 cases already", limits = c(0, 1)) +
  scale_x_continuous(name = expression(paste(R[0]))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
