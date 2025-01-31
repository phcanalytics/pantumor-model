
split_data <- function(assignments){
  # Split the data into training and testing based on the assignment
  train_set <- assignments %>% filter(assignment == "Training")
  test_set <- assignments %>% filter(assignment == "Testing")
  
  # Assign cross-validation folds to the training set
  train_set <- train_set %>%
    assign_cv_folds() %>%
    select(patientid, indication, foldset, fold_id)
  
  return(list(train = train_set, test = test_set))
}

assign_cv_folds <- function(assignments){
  set.seed(7)
  
  # Assign cross-validation folds
  folds <- rsample::vfold_cv(assignments, strata = "indication", v = 5)
  
  # Get fold id per patient for the leave-out fold
  foldids <- 
    purrr::map2_dfr(folds$splits, 
                    folds$id, 
                    function(x, y){
                      rsample::testing(x) %>%
                        dplyr::mutate(Dataset = 'CrossValidatedTraining', Fold = y, FoldSet = 'Validation') %>%
                        dplyr::select(patientid, indication, Dataset, Fold, FoldSet)
                    }
    ) %>%
    dplyr::rename_all(tolower) %>%
    dplyr::filter(foldset == "Validation") %>% 
    dplyr::mutate(fold_id = as.numeric(substr(fold, 5, nchar(fold))))
  
  return(foldids)
}


#' Make model-ready data and CV folds
#' 
#' Make finalized, model-ready data and corresponding CV folds
#' Briefly, this function prepares an object consisting of a data frame of
#' all input features + outcomes, as well as separate train and test object used
#' for future model training and evaluation.
#' 
#' @param data Data frame of data. It should at least contain the columns
#' `patientid` and `indication` as well as `mice.imp` if data were multiply imputed by MICE
#' @param assignments A table of patientids and assignments == "Training" or "Testing"
#' 
#' @return An object of class [`model_data`] with the following:
#' \describe{
#' \item{data}{A data frame of the model-ready dataset.}
#' \item{train_set}{A data frame of the train set containing `patientid`,`indication`, and `fold_id`.}
#'  \item{test_set}{A data frame of the test set containing `patientid` and `indication`.}
#' }
#' 
#' @export

make_model_data <- function(data, assignments){
  
  sets <- split_data(assignments)
  
  # Here you can add additional functions to process `data` as needed
  # If `data` contains multiply-imputed datasets, then you can apply to these on the train and test separately by:
  # df_train <- data %>% 
  #   dplyr::filter(patientid %in% assignments$patientid[assignments$assignment=="Training"]) %>%
  #   dplyr::nest_by(mice.imp) %>%
  #   dplyr::mutate(final_df = 
  #                   list(data %>%
  #                          your_function_here(.))) %>%
  #   dplyr::select(mice.imp, final_df) %>%
  #   tidyr::unnest(cols = c(final_df)) %>%
  #   dplyr::ungroup()
  # 
  # df_test <- data %>% 
  #   dplyr::filter(patientid %in% assignments$patientid[assignments$assignment=="Testing"]) %>%
  #   dplyr::nest_by(mice.imp) %>%
  #   dplyr::mutate(final_df = 
  #                   list(data %>%
  #                          your_function_here(.))) %>%
  #   dplyr::select(mice.imp, final_df) %>%
  #   tidyr::unnest(cols = c(final_df)) %>%
  #   dplyr::ungroup()
  # 
  # data <- bind_rows(df_train, df_test) 
    
  res <- list(
    data = data,
    train_set = sets$train,
    test_set = sets$test
  )
  
  class(res) <- "model_data"
  
  return(res)
  
}
