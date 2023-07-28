
# specify the FOA coefficients for the measurement variance extremes

# LAKE MEAN FOA EQUATIONS

FOA_mean_lake <- function(x){
  predicted_ebu_rate = 25.774853 * 1.065859^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_lake <-  function(x){
  predicted_ebu_rate = 8.022998 * 1.075325^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_lake <-  function(x){
  predicted_ebu_rate = 60.092906 * 1.064871^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_lake <-  function(x){
  predicted_ebu_rate = 3.036065 * 1.176036^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_lake <-  function(x){
  predicted_ebu_rate = 71.765966 * 1.034914^(x-20)
  return(predicted_ebu_rate)
}


# RESERVOIR MEAN FOA EQUATIONS

FOA_mean_res <- function(x){
  predicted_ebu_rate = 25.929964 * 1.029465^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_res <-  function(x){
  predicted_ebu_rate = 3.374613 * 1.032785^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_res <-  function(x){
  predicted_ebu_rate = 67.520539 * 1.026307^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_res <-  function(x){
  predicted_ebu_rate = 5.209963 * 1.068294^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_res <-  function(x){
  predicted_ebu_rate = 64.685000 * 1.022218^(x-20)
  return(predicted_ebu_rate)
}


# specify the model error functions to scale up
FOA_mean_error_low_lake <- function(x){
  predicted_ebu_rate = 91.269685 * 1.0269642^(x-20) - 69.44
  return(predicted_ebu_rate)
}

FOA_mean_error_high_lake <- function(x){
  predicted_ebu_rate = (91.269685 * 1.0269642^(x-20)) + 69.44
  return(predicted_ebu_rate)
}


FOA_mean_error_low_res <- function(x){
  predicted_ebu_rate = 25.92996 * 1.02947^(x-20) - 97.1
  return(predicted_ebu_rate)
}

FOA_mean_error_high_res <- function(x){
  predicted_ebu_rate = 25.92996 * 1.02947^(x-20) + 97.1
  return(predicted_ebu_rate)
}


# specify the coefficient error around the global predictions

FOA_coefficient_low_lake <- function(x){
  predicted_ebu_rate = (21.90539 * 1.04649^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_lake <- function(x){
  predicted_ebu_rate = (29.64431 * 1.08523^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_ebu_rate = (14.05996 * 0.96092^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_ebu_rate = (37.79996 * 1.09802^(x-20))
  return(predicted_ebu_rate)
}
