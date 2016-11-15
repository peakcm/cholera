#### Function to perform revaccination ####

revaccination <- function(t, N, S, V.states.vector, Vax, vac_routine_freq, vac_routine_frac, vac_mig_frac, vac_recip, vac_max, birth_rate, mig_in, vac_birth_frac){
  # In the ODE, the entire compartment cannot move at once, so the vaccine fraction needs to be reduced
  if(vac_routine_frac > 0.99){vac_routine_frac <- 0.99}
  
  vax_rem <- vac_max - Vax

  # Initalize
  routine_S <- 0
  routine_V <- rep(0, length(V.states.vector))
  birth     <- 0
  migrant   <- 0

  if ( as.numeric(vax_rem) > 0 ){

    # Give first priority to births
    if ("birth" %in% vac_recip){
      birth   <- as.numeric(min(vax_rem, birth_rate*N*vac_birth_frac))
      vax_rem <- vax_rem - birth
    }

    # Give second priority to immigrants
    if (vax_rem > 0 && "migrant" %in% vac_recip){
      migrant <- as.numeric(min(vax_rem, mig_in*N*vac_mig_frac))
      vax_rem <- vax_rem - migrant
    }

    # Give third priority to routine vaccination of susceptible people
    if (vax_rem > 0 && ("all" %in% vac_recip | "S" %in% vac_recip)){
      if (t >= 1 && vac_routine_freq*vac_routine_frac > 0 && (round(t) %% vac_routine_freq) == 0){ # If today is the day of a routine revaccination
        routine_S <- as.numeric(min(vax_rem, S*vac_routine_frac))
        vax_rem    <- vax_rem - routine_S
      }
    }

    # Give final priority to routine vaccination of those vaccinated the longest time ago
    if (vax_rem > 0 && "all" %in% vac_recip){
      if (t >= 1 && vac_routine_freq*vac_routine_frac > 0 && (round(t) %% vac_routine_freq) == 0){ # If today is the day of a routine revaccination
        V.states.vector_reverse_cumsum <- rev(cumsum(rev(V.states.vector))) # cumsum from the back of the vector
        routine_V <- as.numeric(V.states.vector* (V.states.vector_reverse_cumsum <= vax_rem)) # Results in a little bit of vaccine wastage
      }
    }
  }
  
  # Convert from the desired number to vaccinate to the desired fraction of each compartment
    # Note that the fraction cannot be 1 because the rate would be infinite. This was adjusted for at the beginning of the function
  if (S > 0){
    routine_S_frac <- as.numeric(routine_S / S)
  } else {routine_S_frac <- 0}
  routine_V_frac <- as.numeric(routine_V/V.states.vector)
  routine_V_frac[which(is.nan(routine_V_frac))] <- 0
  
  # Convert from the desired fraction to the desired transition rate for the ode
  routine_S_rate = as.numeric(S * -log(1-routine_S_frac, base = exp(1)))
  routine_V_rate = as.numeric(V.states.vector * sapply(routine_V_frac, function(x) -log(1-x, base = exp(1))))
  
  return(list(routine_S = routine_S,
              routine_S_rate = routine_S_rate,
              routine_V = routine_V,
              routine_V_rate = routine_V_rate,
              birth     = birth,
              migrant   = migrant))
  }
   

## Simple, old function
# revaccination <- function(t, vac_freq, vac_frac){
#   # Check to make sure vac_freq is greater than zero and that the current time is divisible by the vac_freq
#   if(vac_frac > 0.99){vac_frac <- 0.99}
#   if (t >= 1 && vac_freq > 0 && (round(t) %% vac_freq) == 0){ 
#     return(-log(1-vac_frac, base = exp(1))) # Convert from the fraction that you want to vaccinate to the rate needed
#   } else {return(0)}
# }

## Practice
# revaccination(t = 100,N =  10000,S = 1000,V.states.vector =  V.states.vector,Vax =  100,vac_freq =  365, vac_frac = 1, vax_mig = 1, vac_recip = c("all", "migrant", "birth"), vac_max = 10000, birth_rate = 1/(40*365), mig_in = 1/(20*365))