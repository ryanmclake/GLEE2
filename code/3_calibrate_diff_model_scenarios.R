rm(list=ls())
gc()

library(minpack.lm)
library(tidyverse)

predictNLS <- function(
    object, 
    newdata,
    level = 0.95, 
    nsim = 10000,
    ...
)
{
  require(MASS, quietly = TRUE)
  
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## coefficients
  COEF <- coef(object)
  
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  
  ## check that 'newdata' has same name as predVAR
  if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")
  
  ## get parameter coefficients
  COEF <- coef(object)
  
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
  
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
  
  ## define counter function
  counter <- function (i) 
  {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  
  outMAT <- NULL  
  
  for (i in 1:NR) {
    counter(i)
    
    ## get predictor values and optional errors
    predVAL <- newdata[i, 1]
    if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
    
    ## create mean vector for 'mvrnorm'
    MU <- c(COEF, predVAL)
    
    ## create variance-covariance matrix for 'mvrnorm'
    ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
    
    ## create MC simulation matrix
    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
    
    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
    
    ## collect statistics
    PRED <- data.frame(predVAL)
    colnames(PRED) <- predNAME   
    FITTED <- predict(object, newdata = data.frame(PRED))
    MEAN.sim <- mean(EVAL, na.rm = TRUE)
    SD.sim <- sd(EVAL, na.rm = TRUE)
    MEDIAN.sim <- median(EVAL, na.rm = TRUE)
    MAD.sim <- mad(EVAL, na.rm = TRUE)
    QUANT <- quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
    RES <- c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
    outMAT <- rbind(outMAT, RES)
  }
  
  colnames(outMAT) <- c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
  rownames(outMAT) <- NULL
  
  cat("\n")
  
  return(outMAT)  
}

set.seed(1098)

#How does observation, coefficient, and model variability contribute to total global emissions uncertainty?
### EXPLORE ERROR IN THE OBSERVATIONS FOR THE Q10 TEMP MODELS ###
# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)

diff_base_temp_lake <- vroom::vroom("./output/filtered_GLEE_diffusion_w_GLCP_link_and_ERROR.csv") %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, temp_for_model_C) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.) %>%
  filter(waterbody_type == "lake")

diff_base_temp_res <- vroom::vroom("./output/filtered_GLEE_diffusion_w_GLCP_link_and_ERROR.csv") %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, temp_for_model_C) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "reservoir")


# RES BASELINE
mean_diff_MFOA = nlsLM(ch4_diff ~ A * a^(temp_for_model_C-20),
                      start = list(A = 100, a = 1.1),
                      data = diff_base_temp_lake,
                      control = nls.lm.control(maxiter=1000))

summary(mean_diff_MFOA)

baseline_lake_diff <- as.data.frame(unlist(predictNLS(mean_diff_MFOA, newdata = data.frame(temp_for_model_C = seq(-5,35,0.1)))))

baseline_lake_diff_fit <- as.data.frame(as.numeric(unlist(baseline_lake_diff[,1]))) %>% rename(fit = `as.numeric(unlist(baseline_lake_diff[, 1]))`)
baseline_lake_diff_sd <- as.data.frame(as.numeric(unlist(baseline_lake_diff[,3]))) %>% rename(sd = `as.numeric(unlist(baseline_lake_diff[, 3]))`)
baseline_lake_diff_2.5 <- as.data.frame(as.numeric(unlist(baseline_lake_diff[,6]))) %>% rename(low_CI = `as.numeric(unlist(baseline_lake_diff[, 6]))`)
baseline_lake_diff_97.5 <- as.data.frame(as.numeric(unlist(baseline_lake_diff[,7]))) %>% rename(high_CI = `as.numeric(unlist(baseline_lake_diff[, 7]))`)
newdata = data.frame(temp_for_model_C = seq(-5,35,0.1))

baseline_lake_diff <- cbind(baseline_lake_diff_fit, baseline_lake_diff_sd, baseline_lake_diff_2.5, baseline_lake_diff_97.5, newdata)



# LAKE TEMPORAL 
time_low_diff_MFOA = nlsLM(diff_time_low ~ A * a^(temp_for_model_C-20),
                          start = list(A = 100, a = 1.1),
                          data = diff_base_temp_lake,
                          control = nls.lm.control(maxiter=1000))

summary(time_low_diff_MFOA)

time_low_lake_diff <- as.data.frame(unlist(predictNLS(time_low_diff_MFOA, newdata = data.frame(temp_for_model_C = seq(-5,35,0.1)))))

time_low_lake_diff_fit <- as.data.frame(as.numeric(unlist(time_low_lake_diff[,1]))) %>% rename(fit = `as.numeric(unlist(time_low_lake_diff[, 1]))`)
time_low_lake_diff_sd <- as.data.frame(as.numeric(unlist(time_low_lake_diff[,3]))) %>% rename(sd = `as.numeric(unlist(time_low_lake_diff[, 3]))`)
time_low_lake_diff_2.5 <- as.data.frame(as.numeric(unlist(time_low_lake_diff[,6]))) %>% rename(low_CI = `as.numeric(unlist(time_low_lake_diff[, 6]))`)
time_low_lake_diff_97.5 <- as.data.frame(as.numeric(unlist(time_low_lake_diff[,7]))) %>% rename(high_CI = `as.numeric(unlist(time_low_lake_diff[, 7]))`)
newdata = data.frame(temp_for_model_C = seq(-5,35,0.1))

time_low_lake_diff <- cbind(time_low_lake_diff_fit, time_low_lake_diff_sd, time_low_lake_diff_2.5, time_low_lake_diff_97.5, newdata)
time_low_lake_diff$scenario <- "time_low"

time_high_diff_MFOA = nlsLM(diff_time_high ~ A * a^(temp_for_model_C-20),
                           start = list(A = 100, a = 1.1),
                           data = diff_base_temp_lake,
                           control = nls.lm.control(maxiter=1000))
summary(time_high_diff_MFOA)
time_high_lake_diff <- as.data.frame(unlist(predictNLS(time_high_diff_MFOA, newdata = data.frame(temp_for_model_C = seq(-5,35,0.1)))))

time_high_lake_diff_fit <- as.data.frame(as.numeric(unlist(time_high_lake_diff[,1]))) %>% rename(fit = `as.numeric(unlist(time_high_lake_diff[, 1]))`)
time_high_lake_diff_sd <- as.data.frame(as.numeric(unlist(time_high_lake_diff[,3]))) %>% rename(sd = `as.numeric(unlist(time_high_lake_diff[, 3]))`)
time_high_lake_diff_2.5 <- as.data.frame(as.numeric(unlist(time_high_lake_diff[,6]))) %>% rename(low_CI = `as.numeric(unlist(time_high_lake_diff[, 6]))`)
time_high_lake_diff_97.5 <- as.data.frame(as.numeric(unlist(time_high_lake_diff[,7]))) %>% rename(high_CI = `as.numeric(unlist(time_high_lake_diff[, 7]))`)
newdata = data.frame(temp_for_model_C = seq(-5,35,0.1))

time_high_lake_diff <- cbind(time_high_lake_diff_fit, time_high_lake_diff_sd, time_high_lake_diff_2.5, time_high_lake_diff_97.5, newdata)
time_high_lake_diff$scenario <- "time_high"

time_lake_diff <- rbind(time_low_lake_diff, time_high_lake_diff)


space_low_diff_MFOA = nlsLM(diff_space_low ~ A * a^(temp_for_model_C-20),
                           start = list(A = 100, a = 1.1),
                           data = diff_base_temp_lake,
                           control = nls.lm.control(maxiter=1000))

space_low_lake_diff <- as.data.frame(unlist(predictNLS(space_low_diff_MFOA, newdata = data.frame(temp_for_model_C = seq(-5,35,0.1)))))

space_low_lake_diff_fit <- as.data.frame(as.numeric(unlist(space_low_lake_diff[,1]))) %>% rename(fit = `as.numeric(unlist(space_low_lake_diff[, 1]))`)
space_low_lake_diff_sd <- as.data.frame(as.numeric(unlist(space_low_lake_diff[,3]))) %>% rename(sd = `as.numeric(unlist(space_low_lake_diff[, 3]))`)
space_low_lake_diff_2.5 <- as.data.frame(as.numeric(unlist(space_low_lake_diff[,6]))) %>% rename(low_CI = `as.numeric(unlist(space_low_lake_diff[, 6]))`)
space_low_lake_diff_97.5 <- as.data.frame(as.numeric(unlist(space_low_lake_diff[,7]))) %>% rename(high_CI = `as.numeric(unlist(space_low_lake_diff[, 7]))`)
newdata = data.frame(temp_for_model_C = seq(-5,35,0.1))

space_low_lake_diff <- cbind(space_low_lake_diff_fit, space_low_lake_diff_sd, space_low_lake_diff_2.5, space_low_lake_diff_97.5, newdata)
space_low_lake_diff$scenario <- "space_low"

space_high_diff_MFOA = nlsLM(diff_space_high ~ A * a^(temp_for_model_C-20),
                            start = list(A = 100, a = 1.1),
                            data = diff_base_temp_lake,
                            control = nls.lm.control(maxiter=1000))

space_high_lake_diff <- as.data.frame(unlist(predictNLS(space_high_diff_MFOA, newdata = data.frame(temp_for_model_C = seq(-5,35,0.1)))))

space_high_lake_diff_fit <- as.data.frame(as.numeric(unlist(space_high_lake_diff[,1]))) %>% rename(fit = `as.numeric(unlist(space_high_lake_diff[, 1]))`)
space_high_lake_diff_sd <- as.data.frame(as.numeric(unlist(space_high_lake_diff[,3]))) %>% rename(sd = `as.numeric(unlist(space_high_lake_diff[, 3]))`)
space_high_lake_diff_2.5 <- as.data.frame(as.numeric(unlist(space_high_lake_diff[,6]))) %>% rename(low_CI = `as.numeric(unlist(space_high_lake_diff[, 6]))`)
space_high_lake_diff_97.5 <- as.data.frame(as.numeric(unlist(space_high_lake_diff[,7]))) %>% rename(high_CI = `as.numeric(unlist(space_high_lake_diff[, 7]))`)
newdata = data.frame(temp_for_model_C = seq(-5,35,0.1))

space_high_lake_diff <- cbind(space_high_lake_diff_fit, space_high_lake_diff_sd, space_high_lake_diff_2.5, space_high_lake_diff_97.5, newdata)
space_high_lake_diff$scenario <- "space_high"

space_lake_diff <- rbind(space_low_lake_diff, space_high_lake_diff)


coefficients_a <- rbind(mean_diff_MFOA$m$getPars(),time_low_diff_MFOA$m$getPars(),
                        time_high_diff_MFOA$m$getPars(),space_low_diff_MFOA$m$getPars(),
                        space_high_diff_MFOA$m$getPars())

# RESERVOIRS
mean_diff_MFOA_res = nlsLM(ch4_diff ~ A * a^(temp_for_model_C-20),
                          start = list(A = 100, a = 1.1),
                          data = diff_base_temp_res,
                          control = nls.lm.control(maxiter=1000))

summary(mean_diff_MFOA_res)

time_low_diff_MFOA_res = nlsLM(diff_time_low ~ A * a^(temp_for_model_C-20),
                              start = list(A = 100, a = 1.1),
                              data = diff_base_temp_res,
                              control = nls.lm.control(maxiter=1000))

time_high_diff_MFOA_res = nlsLM(diff_time_high ~ A * a^(temp_for_model_C-20),
                               start = list(A = 100, a = 1.1),
                               data = diff_base_temp_res,
                               control = nls.lm.control(maxiter=1000))

space_low_diff_MFOA_res = nlsLM(diff_space_low ~ A * a^(temp_for_model_C-20),
                               start = list(A = 100, a = 1.1),
                               data = diff_base_temp_res,
                               control = nls.lm.control(maxiter=1000))

space_high_diff_MFOA_res = nlsLM(diff_space_high ~ A * a^(temp_for_model_C-20),
                                start = list(A = 100, a = 1.1),
                                data = diff_base_temp_res,
                                control = nls.lm.control(maxiter=1000))

coefficients_d <- rbind(mean_diff_MFOA_res$m$getPars(),time_low_diff_MFOA_res$m$getPars(),
                        time_high_diff_MFOA_res$m$getPars(),space_low_diff_MFOA_res$m$getPars(),
                        space_high_diff_MFOA_res$m$getPars())
