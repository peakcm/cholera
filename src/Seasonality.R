#### Seasonal Forcing ####
beta_t_fcn <- function(t, shape, amplitude, phase_shift){
  # "times" should be a sequence of days starting from zero
  # "shape" should be "constant" or "sinusoidal"
  # "amplitude" should be 0 if no fluctuations and 1 if beta ranges from 0 to 2*beta
  # "phase shift" in days should be 0 if simulation starts at peak infectiousness 365/2 if at the trough
  
  if (shape == "constant"){
    beta_t <- 1
  } else if (shape == "sinusoidal"){
    beta_t <- 1+amplitude*cos(t/365*2*pi + phase_shift/365*2*pi)
  }
  return(beta_t)
}
