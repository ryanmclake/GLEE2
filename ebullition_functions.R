# specify the FOA coefficients for the measurement variance extremes

FOA_mean_lake <- function(x){
  predicted_ebu_rate = 91.269685 * 1.0269642^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_lake <-  function(x){
  predicted_ebu_rate = 47.042248 * 0.9977927^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_lake <-  function(x){
  predicted_ebu_rate = 981.252156 * 1.0324378^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_lake <-  function(x){
  predicted_ebu_rate = 3.126442 * 1.3298893^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_lake <-  function(x){
  predicted_ebu_rate = 821.953696 * 1.0283409^(x-20)
  return(predicted_ebu_rate)
}

SOA_mean_lake <- function(x, p){
  predicted_ebu_rate = 115.37323 * exp((-0.0021017625*x) * (p/20179.458+p))
  return(predicted_ebu_rate)
}

SOA_time_low_lake <- function(x, p){
  predicted_ebu_rate = 50.82891 * exp((-0.0005940116*x) * (p/5798.549+p))
  return(predicted_ebu_rate)
}

SOA_time_high_lake <- function(x, p){
  predicted_ebu_rate = 889.82795 * exp((0.0005763015*x) * (p/8777186.414+p))
  return(predicted_ebu_rate)
}

SOA_space_low_lake <- function(x, p){
  predicted_ebu_rate = 20.68000 * exp((0.0006338250*x) * (p/-3065.443+p))
  return(predicted_ebu_rate)
}

SOA_space_high_lake <- function(x, p){
  predicted_ebu_rate = 775.78974 * exp((0.0002532683*x) * (p/29899305.937+p))
  return(predicted_ebu_rate)
}

FOA_mean_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_low_res <-  function(x){
  predicted_ebu_rate = 87.82158 * 1.173435^(x-20)
  return(predicted_ebu_rate)
}

FOA_time_high_res <-  function(x){
  predicted_ebu_rate = 50.11549 * 1.632997^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_low_res <-  function(x){
  predicted_ebu_rate = 25.75302 * 1.444188^(x-20)
  return(predicted_ebu_rate)
}

FOA_space_high_res <-  function(x){
  predicted_ebu_rate = 312.80327 * 1.227426^(x-20)
  return(predicted_ebu_rate)
}


SOA_mean_res <- function(x, p){
  predicted_ebu_rate = 927.0834 * exp((-0.03023548*x) * (p/-1.687342e+03+p))
  return(predicted_ebu_rate)
}

SOA_time_low_res <- function(x, p){
  predicted_ebu_rate = 381.2238 * exp((-0.08918259*x) * (p/-1.008395e+04+p))
  return(predicted_ebu_rate)
}

SOA_time_high_res <- function(x, p){
  predicted_ebu_rate = 1925.6869 * exp((-0.01859822*x) * (p/-6.364054e+08+p))
  return(predicted_ebu_rate)
}

SOA_space_low_res <- function(x, p){
  predicted_ebu_rate = 473.8978 * exp((-0.09406905*x) * (p/-3.336566e+03+p))
  return(predicted_ebu_rate)
}

SOA_space_high_res <- function(x, p){
  predicted_ebu_rate = 1443.8990 * exp((-0.01658993*x) * (p/-3.551146e+06+p))
  return(predicted_ebu_rate)
}


# specify the model error functions to scale up
FOA_mean_error_low_lake <- function(x){
  predicted_ebu_rate = (91.269685 * 1.0269642^(x-20)) - 792
  return(predicted_ebu_rate)
}

FOA_mean_error_high_lake <- function(x){
  predicted_ebu_rate = (91.269685 * 1.0269642^(x-20)) + 792
  return(predicted_ebu_rate)
}

SOA_error_low_lake <- function(x, p){
  predicted_ebu_rate = 927.0834 * exp((-0.03023548*x) * (p/-1.687342e+03+p)) - 884
  return(predicted_ebu_rate)
}

SOA_error_high_lake <- function(x, p){
  predicted_ebu_rate = 927.0834 * exp((-0.03023548*x) * (p/-1.687342e+03+p)) + 884
  return(predicted_ebu_rate)
}


FOA_mean_error_low_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20) - 2061
  return(predicted_ebu_rate)
}

FOA_mean_error_high_res <- function(x){
  predicted_ebu_rate = 86.13395 * 1.352015^(x-20) + 2061
  return(predicted_ebu_rate)
}

SOA_error_low_res <- function(x, p){
  predicted_ebu_rate = 927.0834 * exp((-0.03023548*x) * (p/-1.687342e+03+p)) - 1816
  return(predicted_ebu_rate)
}

SOA_error_high_res <- function(x, p){
  predicted_ebu_rate = 927.0834 * exp((-0.03023548*x) * (p/-1.687342e+03+p)) + 1816
  return(predicted_ebu_rate)
}



# specify the coefficient error around the global predictions

FOA_coefficient_low_lake <- function(x){
  predicted_ebu_rate = (49.76578 * 1.068475^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_lake <- function(x){
  predicted_ebu_rate = (68.81703 * 1.105989^(x-20))
  return(predicted_ebu_rate)
}

SOA_coefficient_low_lake <- function(x, p){
  predicted_ebu_rate = 54.89363 * exp((-18.65061*x) * (p/-18075200441+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_high_lake <- function(x, p){
  predicted_ebu_rate = 100.6445 * exp((18.64799*x) * (p/18075202691+p))
  return(predicted_ebu_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_ebu_rate = (-41.79118 * 1.081439^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_ebu_rate = (151.2742 * 1.826618^(x-20))
  return(predicted_ebu_rate)
}

SOA_coefficient_low_res <- function(x, p){
  predicted_ebu_rate = 772.9693 * exp((-67340.63*x) * (p/-5833387514+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_high_res <- function(x, p){
  predicted_ebu_rate = 1058.831 * exp((67340.58*x) * (p/5833387612+p))
  return(predicted_ebu_rate)
}
