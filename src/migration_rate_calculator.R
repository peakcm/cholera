#### Function to update migration rates ####
migration_rate_calculator <- function(t, rate_vector){
  t_floor <- max(1,floor(t))
  rate_vector[t_floor]
}

