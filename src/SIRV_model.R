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
  lambda <- as.numeric(params$beta * I/N * beta_t_fcn(t, params$beta_shape, params$beta_amp, params$beta_phase_shift))
  
  ## Update migration rates
  if (params$mig_rates_constant == TRUE){
    mig_in <- params$mig_in
    mig_out <- params$mig_out
  } else{
    mig_in <- migration_rate_calculator(t = t, rate_vector = params$mig_in)
    mig_out <- migration_rate_calculator(t = t, rate_vector = params$mig_out)
  }
  
  ## revaccination
  if ("revaccination_log_trans" %in% names(params)){ # Added this feature much later, so this check will make sure it only applies when I define it.
    revax_out <- revaccination(t = t, N = N, S = S, V.states.vector =  V.states.vector, Vax =  Vax, vac_routine_count = params$vac_routine_count, vac_mass_freq =  params$vac_mass_freq, vac_mass_frac = params$vac_mass_frac, vac_mig_frac = params$vac_mig_frac, vac_recip = params$vac_recip, vac_max = params$vac_max, birth_rate = params$birth_death_rate, mig_in = mig_in, vac_birth_frac = params$vac_birth_frac, vac_stopper = params$vac_stopper, revaccination_log_trans = params$revaccination_log_trans)
  } else {
    revax_out <- revaccination(t = t, N = N, S = S, V.states.vector =  V.states.vector, Vax =  Vax, vac_routine_count = params$vac_routine_count, vac_mass_freq =  params$vac_mass_freq, vac_mass_frac = params$vac_mass_frac, vac_mig_frac = params$vac_mig_frac, vac_recip = params$vac_recip, vac_max = params$vac_max, birth_rate = params$birth_death_rate, mig_in = mig_in, vac_birth_frac = params$vac_birth_frac, vac_stopper = params$vac_stopper)
    }
  
  ## vaccine waning
  dV.states <- rep(0, params$n.comps.V)
  names(dV.states) <- paste0(rep(paste0("dV."), params$n.comps.V), 1:params$n.comps.V)
  
  if (length(dV.states) != length(params$VE)) stop("length of params$VE does not match number of V.states")
  
    # Special case for the first vaccine state because it receives new vaccinees
  dV.states[1] <- -params$V_step*V.states.vector[1] -lambda*V.states.vector[1]*(1-params$VE[1]) -mig_out*V.states.vector[1] -params$birth_death_rate*V.states.vector[1] -revax_out$mass_V_rate[1] +sum(revax_out$routine_S_rate, revax_out$routine_V_rate, revax_out$mass_S_rate, revax_out$mass_V_rate, revax_out$birth, revax_out$migrant)
  
    # For other vaccine states, need to create a temporary data frame that lines up the current V members and the lagged members so we can use an apply function to advance people instead of using a loop. See next section with the looped version (commented out)
  temp.V.states.df <- data.frame(orig = as.numeric(V.states.vector[2:params$n.comps.V]), lag_1 = as.numeric(V.states.vector[1:(params$n.comps.V-1)]), VE = params$VE[2:params$n.comps.V], mass_V_rate = revax_out$mass_V_rate[2:params$n.comps.V], routine_V_rate = revax_out$routine_V_rate[2:params$n.comps.V])
  dV.states[2:params$n.comps.V] <- apply(temp.V.states.df, 1, function(x) params$V_step*x["lag_1"] -params$V_step*x["orig"] -lambda*x["orig"]*(1-x["VE"]) -mig_out*x["orig"] -params$birth_death_rate*x["orig"] -x["routine_V_rate"] -x["mass_V_rate"])
  rm("temp.V.states.df")
  
    # Slower version of updating states using loop instead of apply
  # for (y in seq_len(params$n.comps.V)){
  #   if (y==1){
  #     dV.states[1] <- -params$V_step*V.states.vector[1] -lambda*V.states.vector[1]*(1-params$VE[1]) -mig_out*V.states.vector[1] -params$birth_death_rate*V.states.vector[1] -revax_out$mass_V_rate[1] +sum(revax_out$mass_S_rate, revax_out$mass_V_rate, revax_out$birth, revax_out$migrant)
  #   } else {
  #     dV.states[y] <- +params$V_step*V.states.vector[y-1] -params$V_step*V.states.vector[y] -lambda*V.states.vector[y]*(1-params$VE[y])  -mig_out*V.states.vector[y] -params$birth_death_rate*V.states.vector[y] -revax_out$mass_V_rate[y]
  #   }
  # }
  
  ## susceptibles
  dS <- as.numeric(params$nat_wane*R + as.numeric(params$V_step*V.states.vector[params$n.comps.V]) +mig_in*(N*(1-params$foreign_infection)) +params$birth_death_rate*N -lambda*S -mig_out*S -params$birth_death_rate*S -revax_out$mass_S_rate -revax_out$birth - revax_out$migrant -revax_out$routine_S_rate)
  
  ## latent
  dE  <- as.numeric(- params$sigma*E + lambda*S + sum(lambda*V.states.vector*(1-params$VE)) - mig_out*E -params$birth_death_rate*E)
  
  ## infectious 
  dI  <-  as.numeric(params$sigma*E + mig_in*(N*params$foreign_infection) - params$gamma*I -mig_out*I -params$birth_death_rate*I)
  
  ## removed 
  dR  <- as.numeric(params$gamma*I -params$nat_wane*R -mig_out*R -params$birth_death_rate*R)
  
  ## Water 
  dW <- 0
  
  ## Cumulative Incidence 
  dCI  <- +as.numeric(params$sigma*E) #lambda*S #CI in unvaccinated class
  
  ## Number of vaccine courses
  dVax <- +sum(revax_out$routine_S_rate, revax_out$routine_V_rate, revax_out$mass_S_rate, revax_out$mass_V_rate, revax_out$birth, revax_out$migrant)
  # cat("dVax:", dVax, "\n")
  # cat("mass_S_rate:",revax_out$mass_S_rate,"\n")
  # cat("dV.states:", dV.states, "\n")

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


