# Bootstrap model data to estimate confidence intervals ---------------------------
#' Bootstrap model data
#' 
#' This function is designed to bootstrap the train and test data separately.
#' Each dataset is sampled with replacement based on a seed set by bootstrap iteration `b`.
#' MICE datasets are resampled identically.
#' 
#' @param seed Seed to set for sampling with replacement. In practice, `seed`=`b`, 
#' where `b`=bootstrap iteration.
#' @param data Data frame of data. It should at least contain the columns
#' `patientid` and `indication` as well as `mice.imp` if data were multiply imputed by MICE
#' @param assignments A table of patientids and assignments == "Training" or "Testing"
#' 
#' @return A data frame with the same columns as `data_mice`, but containing the bootstrap
#' resampled data. Some records will be duplicated. 
#' @export

bootstrap_data <- function(seed, data, assignments){
  
  # Get train and test patient ids to resample
  train_patientid <- assignments %>%
    dplyr::filter(assignment == "Training") %>%
    dplyr::pull(patientid)
  
  test_patientid <- assignments %>%
    dplyr::filter(assignment == "Testing") %>%
    dplyr::pull(patientid)
  
  # Need to resample by mice dataset
  all_samples <- list()
  for(j in unique(data$mice.imp)){
    # Subset to train/test
    data_train <- data %>% filter(patientid%in%train_patientid) %>% filter(mice.imp==j)
    data_test <- data %>% filter(patientid%in%test_patientid) %>% filter(mice.imp==j)
    
    set.seed(seed)
    data_train <- data_train[sample(nrow(data_train), replace=TRUE),]
    
    set.seed(seed)
    data_test <- data_test[sample(nrow(data_test), replace=TRUE),]
    
    # Recombine to produce single bootstrapped dataset (train + test)
    all_samples[[j]] = rbind(data_train, data_test)
  }
  
  data_resample <- do.call(rbind, all_samples)
  return(data_resample)
}