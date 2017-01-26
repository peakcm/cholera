#### Function to perform revaccination ####

revaccination <- function(t, N, S, V.states.vector, Vax, vac_routine_count, vac_mass_freq, vac_mass_frac, vac_mig_frac, vac_recip, vac_max, birth_rate, mig_in, vac_birth_frac, vac_stopper=1e20, revaccination_log_trans = TRUE){
  # In the ODE, the entire compartment cannot move at once, so the vaccine fraction needs to be reduced
  if(vac_mass_frac > 0.99){vac_mass_frac <- 0.99}
  
  # cat("Vax:", Vax, "\n") # Troubleshoot
  vax_rem <- vac_max - Vax
  if (vac_routine_count > vax_rem){
    vac_routine_count <- vax_rem
  }
  
  # Check if vac_stopper < t, then set vax_rem to zero. Useful only for the Bentiu case study
  if (vac_stopper < t){
    vax_rem <- 0
  }

  # Initalize
  birth     <- 0
  migrant   <- 0
  routine_S   <- 0
  routine_V <- rep(0, length(V.states.vector))
  mass_S <- 0
  mass_V <- rep(0, length(V.states.vector))

  if ( as.numeric(vax_rem) > 0 ){

    # Give first priority to births
    if ("birth" %in% vac_recip){
      birth   <- min(vac_routine_count, as.numeric(min(vax_rem, birth_rate*N*vac_birth_frac)))
      vax_rem <- vax_rem - birth
      vac_routine_count <- vac_routine_count - birth
    }

    # Give second priority to immigrants
    if (vax_rem > 0 && "migrant" %in% vac_recip){
      migrant <- min(vac_routine_count,as.numeric(min(vax_rem, mig_in*N*vac_mig_frac)))
      vax_rem <- vax_rem - migrant
      vac_routine_count <- vac_routine_count - migrant
    }
    
    # Give third priority to routine vaccination of susceptible people
    if (vax_rem > 0 && ("routine_all" %in% vac_recip | "routine_S" %in% vac_recip) && S > 1){
      routine_S <- min(max(0,S-1), as.numeric(vac_routine_count)) # Remove one S so that the fraction is never 1
      vax_rem <- vax_rem - routine_S
      vac_routine_count <- vac_routine_count - routine_S
    }
    
    # Give fourth priority to routine vaccination of those vaccinated the longest time ago
    if (vax_rem > 0 && "routine_all" %in% vac_recip && sum(V.states.vector)>0){
      V.states.vector_reverse_cumsum <- as.numeric(rev(cumsum(rev(V.states.vector)))) # cumsum from the back of the vector
      V_to_vaccinate <- V.states.vector_reverse_cumsum < as.numeric(vac_routine_count)
      
      routine_V <- as.numeric(V.states.vector * V_to_vaccinate) # Results in a little bit of vaccine wastage
      
      vax_rem <- vax_rem - sum(routine_V)
      vac_routine_count <- vac_routine_count - sum(routine_V)
    }

    # Give fifth priority to mass vaccination of susceptible people
    if (vax_rem > 0 && ("mass_all" %in% vac_recip | "mass_S" %in% vac_recip)){
      if (t >= 1 && vac_mass_freq*vac_mass_frac > 0 && (round(t) %% vac_mass_freq) == 0){ # If today is the day of a mass revaccination
        mass_S <- as.numeric(min(vax_rem, (S-routine_S)*vac_mass_frac))
        vax_rem    <- vax_rem - mass_S
      }
    }

    # Give sixth priority to mass vaccination of those vaccinated the longest time ago
    if (vax_rem > 0 && "mass_all" %in% vac_recip){
      if (t >= 1 && vac_mass_freq*vac_mass_frac > 0 && (round(t) %% vac_mass_freq) == 0){ # If today is the day of a mass revaccination
        V.states.vector_reverse_cumsum <- rev(cumsum(rev((V.states.vector-routine_V)*vac_mass_frac))) # cumsum from the back of the vector
        mass_V <- as.numeric(V.states.vector*vac_mass_frac* (V.states.vector_reverse_cumsum < as.numeric(vax_rem))) # Results in a little bit of vaccine wastage
      }
    }
  }
  
  # Convert from the desired number to vaccinate to the desired fraction of each compartment
    # Note that the fraction cannot be 1 because the rate would be infinite. This was adjusted for at the beginning of the function
  if (S > 0){
    routine_S_frac <- as.numeric(routine_S / S)
    mass_S_frac <- as.numeric(mass_S / S)
  } else {
    routine_S_frac <- 0
    mass_S_frac <- 0
  }
  routine_V_frac <- as.numeric(routine_V/V.states.vector)
  routine_V_frac[which(is.nan(routine_V_frac))] <- 0
  routine_V_frac[which(routine_V_frac > 0.99)] <- 0.99
  mass_V_frac <- as.numeric(mass_V/(V.states.vector-routine_V))
  mass_V_frac[which(is.nan(mass_V_frac))] <- 0
  
  # Convert from the desired fraction to the desired transition rate for the ode
  routine_S_rate = as.numeric(S * -log(1-routine_S_frac, base = exp(1)))
  mass_S_rate = as.numeric((S-routine_S) * -log(1-mass_S_frac, base = exp(1)))
  routine_V_rate = as.numeric(V.states.vector * sapply(routine_V_frac, function(x) -log(1-x, base = exp(1))))
  mass_V_rate = as.numeric((V.states.vector-routine_V) * sapply(mass_V_frac, function(x) -log(1-x, base = exp(1))))
  
  if (revaccination_log_trans == FALSE){ # Remove the log transformation of the rates
    routine_S_rate <- routine_S
    routine_V_rate <- routine_V
    mass_S_rate <- mass_S
    mass_V_rate <- mass_V
  }
  
  # cat(routine_V_rate, "\n") # Troubleshoot
  return(list(routine_S = routine_S,
              routine_S_rate = routine_S_rate,
              routine_V = routine_V,
              routine_V_rate = routine_V_rate,
              mass_S = mass_S,
              mass_S_rate = mass_S_rate,
              mass_V = mass_V,
              mass_V_rate = mass_V_rate,
              birth     = birth,
              migrant   = migrant))
  }
   

## Simple, old function
# revaccination <- function(t, vac_mass_freq, vac_mass_frac){
#   # Check to make sure vac_freq is greater than zero and that the current time is divisible by the vac_mass_freq
#   if(vac_mass_frac > 0.99){vac_mass_frac <- 0.99}
#   if (t >= 1 && vac_mass_freq > 0 && (round(t) %% vac_mass_freq) == 0){ 
#     return(-log(1-vac_mass_frac, base = exp(1))) # Convert from the fraction that you want to vaccinate to the rate needed
#   } else {return(0)}
# }

## Practice
# revaccination(t = 100,N =  10000,S = 1000,V.states.vector =  V.states.vector,Vax =  100,vac_mass_freq =  365, vac_mass_frac = 1, vax_mig = 1, vac_recip = c("all", "migrant", "birth"), vac_max = 10000, birth_rate = 1/(40*365), mig_in = 1/(20*365))