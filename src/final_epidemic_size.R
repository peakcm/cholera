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
# N = 20000
# R = 1.5
# q <- final_epidemic_size(N, R)
# cat("The probability of <10 cases is",sum(q[1:10]))
# cat("The probability of 10 to 50 cases is", sum(q[11:(N/3)]))
# cat("The probability of >50 cases is", sum(q[(N/3):N]))
# plot(q*N, type = "h")
# plot(q*N, type = "h", xlim = c(0,30))
# 
# # Probability of a large outbreak (>1% AR) given you make it past "cut_point" cases
# cut_point <- 10
# sum(q[(N/100):N]) / sum(q[(cut_point+1):N])
# When N=20,000 and R=1.5, cut_point 10 is 89% predictive. cut_point 25 is 97% predictive
# When N=10,000 and R=1.5, cut_point 10 is 88% predictive. cut_point 25 is 97% predictive
# When N=10,000 and R=1.4, cut_point 10 is 82% predictive. cut_point 25 is 95% predictive
# When N=10,000 and R=1.1, cut_point 10 is 47% predictive. cut_point 25 is 68% predictive. cut_point 80 is 95% predictive....but not very interesting 
