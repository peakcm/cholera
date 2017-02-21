################
# Input VE Fcn #
################

Create_VE <- function(timesteps_per_month,VE_shape = "Shanchol",bound=TRUE, max_V_months, custom_input){
  if (VE_shape == "Shanchol"){
    input <- matrix(c(0,6,12,18,24,30,36,42,48,.430,.525,.579,.583,.538,.438,.245,-.073,-.590),nrow=9)
  } else if (VE_shape == "Dukoral"){
    input <- matrix(c(0,6,12,18,24,30,36,42,48,54,60,.713,.650,.572,.476,.374,.280,.202,.141,.092,0.046,0),nrow=11)
  } else if (VE_shape == "Linear"){
    input <- matrix(c(0, max_V_months, .713, 0), nrow = 2)
  } else if (VE_shape == "Custom"){
    # input should be a matrix with first column of months and a second column of VE at that month
    # first entry time is at month 0
    input <- custom_input
  } else if (VE_shape == "Perfect"){
    input <- matrix(c(0,max_V_months, 1, 1), nrow=2)
  }
  
  last_month = max(max(input[,1]), max_V_months)
  VE <- matrix(c(seq(0, last_month, by = 1/timesteps_per_month),rep(0,1+last_month*timesteps_per_month)),nrow=1+last_month*timesteps_per_month)
  VE[VE[,1] %in% input[,1],2] <- input[,2]
  
  # estimate VE for other months (linear)
  for (i in seq_len(nrow(input)-1)){
    diff <- input[i+1,2] - input[i,2]
    VE_rows_to_fill <- which(VE[,1] > input[i,1] & VE[,1] < input[i+1,1])
    for (j in VE_rows_to_fill){
      VE[j,2] <- VE[j-1,2]+diff/(1+length(VE_rows_to_fill))
    }
  }
  
  if (nrow(VE) > max_V_months){
    VE <- VE[VE[,1] <= max_V_months,]
  }
  
  # set negative results to zero if bound=TRUE
  if (!bound){
    print("Some created VE estimates are negative")
  } else{
    VE[VE[,2]<0,2] <- 0
  }
  return(VE[2:nrow(VE),2])
}

# Create Shancol and Dukoral Estimates
# Shanchol <- matrix(c(0,6,12,18,24,30,36,42,48,.430,.525,.579,.583,.538,.438,.245,-.073,-.590),nrow=9)
# Dukoral  <- matrix(c(0,6,12,18,24,30,36,42,48,54,60,.713,.650,.572,.476,.374,.280,.202,.141,.092,0.046,0),nrow=11)
# VE_Shanchol <- Create_VE(timesteps_per_month = 1,VE_shape = "Shanchol",bound = TRUE, max_V_months = 60)
# VE_Dukoral  <- Create_VE(timesteps_per_month = 1,VE_shape = "Dukoral",bound = TRUE, max_V_months = 60)
