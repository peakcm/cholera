#### Function to perform revaccination ####

revaccination <- function(t, vac_freq, vac_frac){
  # Check to make sure vac_freq is greater than zero and that the current time is divisible by the vac_freq
  if(vac_frac > 0.99){vac_frac <- 0.99}
  if (t >= 1 && vac_freq > 0 && (round(t) %% vac_freq) == 0){ 
    return(-log(1-vac_frac, base = exp(1))) # Convert from the fraction that you want to vaccinate to the rate needed
  } else {return(0)}
}