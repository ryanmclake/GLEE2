library(dplyr)
library(readr)
library(vroom)


vroom::vroom("/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_BASE.csv") %>%
  dplyr::group_by(centr_lat, centr_lon) %>%
  dplyr::summarise(`Baseline Ebullition Predictions` = (sum(FOA_mean_ebullition))/5021) %>%
  readr::write_csv(., "/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_BASE_mean.csv")

vroom::vroom("/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_TIME.csv") %>%
  dplyr::group_by(centr_lat, centr_lon) %>%
  dplyr::summarise(`Time-low Ebu Predictions` = (sum(FOA_time_low_ebu))/5021,
            `Time-high Ebu Predictions` = (sum(FOA_time_high_ebu))/5021) %>%
  readr::write_csv(., "/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_TIME_mean.csv")

vroom::vroom("/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_SPACE.csv") %>%
  dplyr::group_by(centr_lat, centr_lon) %>%
  dplyr::summarise(`Space-low Ebu Predictions` = (sum(FOA_space_low_ebu))/5021,
            `Space-high Ebu Predictions` = (sum(FOA_space_high_ebu))/5021) %>%
  readr::write_csv(., "/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_SPACE_mean.csv")

vroom::vroom("/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_MODEL.csv") %>%
  dplyr::group_by(centr_lat, centr_lon) %>%
  dplyr::summarise(`Model-low Ebu Predictions` = (sum(FOA_mean_error_low_ebu))/5021,
            `Model-high Ebu Predictions` = (sum(FOA_mean_error_high_ebu))/5021) %>%
  readr::write_csv(., "/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_MODEL_mean.csv")

vroom::vroom("/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_PARAM.csv") %>%
  dplyr::group_by(centr_lat, centr_lon) %>%
  dplyr::summarise(`Param-low Ebu Predictions` = (sum(FOA_coefficient_low_ebu))/5021,
            `Param-high Ebu Predictions` = (sum(FOA_coefficient_high_ebu))/5021) %>%
  readr::write_csv(., "/central/groups/carnegie_poc/rmcclure/project-GLEE/output/Global_ebullition_variability_PARAM_mean.csv")

