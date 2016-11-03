#### Function to run SIRV model ####
run_model <- function(func, times, params){
  require(deSolve)
  
  run <- ode(y=x,
             times = times,
             func = func,
             parms = params)
  run <- data.frame(run)
  names(run) <- c("time", "S", paste0("V.", 1:n.comps.V), "E", "I", "R", "W", "CI")
  run$V_total <- apply(run[,3:(2+n.comps.V)], 1, sum)
  run$Re <- apply(run, 1, function(x) calculate_Re(x, params))
  run$prob_outbreak_10 <- prob_outbreak_fcn(R = run$Re, outbreak_size = 10)
  
  return(run)
}

