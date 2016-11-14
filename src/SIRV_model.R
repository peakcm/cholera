#### Header ####
# Create an SIRV model that can account for many V compartments
# Adapted from Azman et al. https://github.com/HopkinsIDD/singledose-ocv

#### Adapted function ####
SIRV.generic <- function(t,
                         x, 
                         params){
  
  ## get the names of each of the Vaccine compartments
  V.states <- paste0(rep(paste0("V."), params$n.comps.V), 1:params$n.comps.V)
  
  ## names of all the states
  state.vars.all <- c("S",
                      V.states,
                      "E",
                      "I",
                      "R",
                      "W",
                      "CI",
                      "Vax")
  
  if (length(state.vars.all) != length(x)) stop("input state length not correct")
  
  dV.1 <- numeric(params$n.comps.V) # Initialize with zeros

  ## assign our state variables within the environment
  for (y in seq_along(state.vars.all)) {assign(state.vars.all[y],x[y])}
  V.states.vector <- x[seq(1,params$n.comps.V)+1]
  names(V.states.vector) <- state.vars.all[seq(1,params$n.comps.V)+1]
  
  N  <- sum(x[1:(4+params$n.comps.V)]) # Total pop 

  ## force of infection (partitioned)
  lambda <- params$beta * I/N * beta_t_fcn(t, params$beta_shape, params$beta_amp, params$beta_phase_shift)
  
  ## revaccination
  revax <- N*revaccination(t, params$vac_freq, params$vac_frac) # Number of vaccines to be allocated today
  if ( Vax > params$max_vax){
    vax_remaining <- 0
  } else {
    vax_remaining <- 1
    revax <- min(revax, params$max_vax - Vax)
    }
  # print(c(t, revax))

  ## vaccine waning
  dV.states <- rep(0, params$n.comps.V)
  names(dV.states) <- paste0(rep(paste0("dV."), params$n.comps.V), 1:params$n.comps.V)
  
  if (length(dV.states) != length(params$VE)) stop("length of params$VE does not match number of V.states")
  
  for (y in seq_len(params$n.comps.V)){
    if (y==1){
      dV.states[1] <- -params$V_step*V.states.vector[1] -lambda*V.states.vector[1]*(1-params$VE[1]) -params$mig_out*V.states.vector[1] -params$birth_death_rate*V.states.vector[1] -revax*V.states.vector[1]/N*vax_remaining +revax*vax_remaining +params$mig_in*(N*(1-params$foreign_infection))*(params$vax_mig)*vax_remaining
    } else {
      dV.states[y] <- +params$V_step*V.states.vector[y-1] -params$V_step*V.states.vector[y] -lambda*V.states.vector[y]*(1-params$VE[y])  -params$mig_out*V.states.vector[y] -params$birth_death_rate*V.states.vector[y] -revax*V.states.vector[y]/N*vax_remaining
    }
  }
  
  ## susceptibles
  dS <- params$nat_wane*R + as.numeric(params$V_step*V.states.vector[params$n.comps.V]) +params$mig_in*(N*(1-params$foreign_infection))*(1-params$vax_mig) +params$birth_death_rate*N -lambda*S -params$mig_out*S -params$birth_death_rate*S -revax*S/N*vax_remaining
  
  ## latent
  dE  <- - params$sigma*E + lambda*S + sum(lambda*V.states.vector*(1-params$VE)) - params$mig_out*E -params$birth_death_rate*E -revax*E/N*vax_remaining
  
  ## infectious 
  dI  <-  params$sigma*E + params$mig_in*(N*params$foreign_infection) - params$gamma*I -params$mig_out*I -params$birth_death_rate*I -revax*I/N*vax_remaining
  
  ## removed 
  dR  <- params$gamma*I -params$nat_wane*R -params$mig_out*R -params$birth_death_rate*R -revax*R/N*vax_remaining
  
  ## Water 
  dW <- 0
  
  ## Cumulative Incidence 
  dCI  <- params$sigma*E #lambda*S #CI in unvaccinated class
  
  ## Number of vaccine courses
  dVax <- revax*vax_remaining + params$mig_in*(N*(1-params$foreign_infection))*(params$vax_mig)*vax_remaining
  
  out  <- c(dS = dS,
            dV.states,
            dE = dE,
            dI = dI,
            dR = dR,
            dW = dW,
            dCI = dCI,
            dVax = dVax)
  
  return(list(out))
}


