#### Estimate probabilty of an outbreak given introduction of a case ####

# See equation in "Modeling to Inform Infectious Disease Control" by Niels G. Becker on page 16
# http://lib.myilibrary.com.ezp-prod1.hul.harvard.edu/Open.aspx?id=772197

prob_outbreak_fcn <- function(R, outbreak_size=10, distribution = "Poisson"){
  if (distribution == "Poisson"){
    out = 0
    for (y in seq_len(outbreak_size)){
      out <- out + 1/factorial(y-1) * y^(y-2) * R^(y-1) * exp(-y*R)

    }
  }
  return(1-out)
}

# par(mar = c(6,6,6,6))
# plot(seq(0, 3, by=0.01), prob_outbreak_fcn(seq(0, 3, by=0.01), outbreak_size = 5), type="l", ylim = c(0,1), col = "darkgreen", xlab = "R", ylab = "Probability of an Outbreak Greater Than...\n5, 10, or 50 cases")
# lines(seq(0, 3, by=0.01), prob_outbreak_fcn(seq(0, 3, by=0.01), outbreak_size = 10), col = "orange")
# lines(seq(0, 3, by=0.01), prob_outbreak_fcn(seq(0, 3, by=0.01), outbreak_size = 50), col = "darkred")

