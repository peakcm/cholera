# Adapted from http://www.sciencedirect.com/science/article/pii/S0022519314006882

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
# q <- final_epidemic_size(1000, 1.1)
# cat("The probability of <10 cases is",sum(q[1:10]))
# cat("The probability of 10 to 50 cases is", sum(q[11:51]))
# cat("The probability of >50 cases is", sum(q[52:N]))
# hist(q*N)

