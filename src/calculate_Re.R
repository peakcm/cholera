# Helper function to calculate Re

calculate_Re <- function(simulation_output, params){
  X  <- sum(simulation_output["S"], sum(simulation_output[3:(params$n.comps.V+2)]*(1-params$VE))) / sum(simulation_output[2:(2+params$n.comps.V+3)])
  R0 <- params$beta / params$gamma * beta_t_fcn(simulation_output["time"], shape = params$beta_shape, amplitude = params$beta_amp, phase_shift = params$beta_phase_shift)
  Re <- X*R0
  return(Re)
}
