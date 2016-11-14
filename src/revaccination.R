#### Function to perform revaccination ####

revaccination <- function(t, N, S, V.states.vector, Vax, vac_routine_freq, vac_routine_frac, vac_mig_frac, vac_recip, vac_max, birth_rate, mig_in, vac_birth_frac){
  # # Check to make sure vac_freq is greater than zero and that the current time is divisible by the vac_freq
  # if(vac_frac > 0.999){vac_frac <- 0.999}
  vax_rem <- vac_max - Vax
  
  # Initalize
  routine_S <- 0
  routine_V <- rep(0, length(V.states.vector))
  birth     <- 0
  migrant   <- 0
  
  if ( vax_rem > 0 ){
    
    # Give first priority to births
    if ("birth" %in% vac_recip){
      birth   <- min(vax_rem, birth_rate*N*vac_birth_frac)
      vax_rem <- vax_rem - birth
    }
    
    # Give second priority to immigrants
    if (vax_rem > 0 && "migrant" %in% vac_recip){
      migrant <- min(vax_rem, mig_in*N*vac_mig_frac)
      vax_rem <- vax_rem - migrant
    }
    
    # Give third priority to routine vaccination of susceptible people
    if (vax_rem > 0 && ("all" %in% vac_recip | "S" %in% vac_recip)){
      if (t >= 1 && vac_routine_freq*vac_routine_frac > 0 && (round(t) %% vac_routine_freq) == 0){ # If today is the day of a routine revaccination
        routine_S <- min(vax_rem, S*vac_routine_frac)
        vax_rem    <- vax_rem - routine_S
      }
    }
    
    # Give final priority to routine vaccination of those vaccinated the longest time ago
    if (vax_rem > 0 && "all" %in% vac_recip){
      if (t >= 1 && vac_routine_freq*vac_routine_frac > 0 && (round(t) %% vac_routine_freq) == 0){ # If today is the day of a routine revaccination
        V.states.vector_reverse_cumsum <- rev(cumsum(rev(V.states.vector))) # cumsum from the back of the vector
        routine_V <- V.states.vector* (V.states.vector_reverse_cumsum <= vax_rem) # Results in a little bit of vaccine wastage
      }
    }
  }
  
  return(list(routine_S = routine_S,
              routine_V = routine_V,
              birth     = birth,
              migrant   = migrant))
  }
    
    # -log(1-vac_frac, base = exp(1)) # Convert from the fraction that you want to vaccinate to the rate needed

# revaccination(t = 100,N =  10000,S = 1000,V.states.vector =  V.states.vector,Vax =  100,vac_freq =  365, vac_frac = 1, vax_mig = 1, vac_recip = c("all", "migrant", "birth"), vac_max = 10000, birth_rate = 1/(40*365), mig_in = 1/(20*365))