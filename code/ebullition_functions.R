# specify the FOA coefficients for the measurement variance extremes

# LAKE MEAN FOA EQUATIONS

FOA_mean_lake <- function(x){
  predicted_ebu_rate = 91.269685 * 1.0269642^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_lake <-  function(x){
  predicted_ebu_rate = 47.042139 * 0.9977923^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_lake <-  function(x){
  predicted_ebu_rate = 201.892220 * 1.0496897^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_lake <-  function(x){
  predicted_ebu_rate = 3.126347 * 1.3298923^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_lake <-  function(x){
  predicted_ebu_rate = 265.229264 * 1.0309585^(x-20)
  return(predicted_ebu_rate)
}


# RESERVOIR MEAN FOA EQUATIONS

FOA_mean_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_res <-  function(x){
  predicted_ebu_rate = 87.821582 * 1.173435^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_res <-  function(x){
  predicted_ebu_rate = 1.629564 * 2.517342^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_res <-  function(x){
  predicted_ebu_rate = 25.753039 * 1.444188^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_res <-  function(x){
  predicted_ebu_rate = 222.604025 * 1.270625^(x-20)
  return(predicted_ebu_rate)
}


# specify the model error functions to scale up
FOA_mean_error_low_lake <- function(x){
  predicted_ebu_rate = 91.269685 * 1.0269642^(x-20) - 125
  return(predicted_ebu_rate)
}

FOA_mean_error_high_lake <- function(x){
  predicted_ebu_rate = (91.269685 * 1.0269642^(x-20)) + 125
  return(predicted_ebu_rate)
}


FOA_mean_error_low_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20) - 438.3
  return(predicted_ebu_rate)
}

FOA_mean_error_high_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20) + 438.3
  return(predicted_ebu_rate)
}


# specify the coefficient error around the global predictions

FOA_coefficient_low_lake <- function(x){
  predicted_ebu_rate = (79.81968 * 1.00997^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_lake <- function(x){
  predicted_ebu_rate = (102.7197 * 1.04395^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_ebu_rate = (9.9645 * 1.1524^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_ebu_rate = (182.2323 * 1.5516^(x-20))
  return(predicted_ebu_rate)
}
