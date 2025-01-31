#' Survival modeling pipeline described in McGough et al. 
#' 
#' Run a single instance of a survival modeling pipeline by sequentially running
#'  [preprocess()], [fit()], [predict()], and [evaluate()]. The pipeline is 
#'  run according to the scenario described by `config`.
#'
#' @param config A configuration file describing all modeling scenarios of interest.
#' @param m The 'm' dataset selected from multiple imputation.
#' @param row The row from the configuration file used when running the pipeline.
#' @param model_data Model-ready data containing train and test splits. Must be
#' an object of class `model_data`, see `make_model_data()`.
#' @param return_fits Logical, return the fitted model objects from the run? Default TRUE. Can result in large object sizes.
#' @param verbose Integer value for level of verboseness of logging. Default to 0 - none
#' 1 - ERROR, 2 - WARN, 3 - INFO, 4 - DEBUG
#'
#' @return A [`dplyr::tibble`] with one row and the following columns:
#' \describe{
#' \item{scenario}{An integer valued ID for the modeling scenario.}
#' \item{fit}{The fitted survival model from [fit()].}
#' \item{prediction}{Predictions from [predict_surv()].}
#' \item{evaluation}{Model evaluations from [evaluate()].}
#' }
#' 
#' @export
run_pipeline <- function(config, m, row, model_data, return_fits = TRUE, verbose = 0){
  log_lvl <- (verbose + 1)*100
  class(log_lvl) <- "loglevel"
  logger::log_threshold(log_lvl)
  config_json <- jsonlite::toJSON(config[row,])
  
  xy <- preprocess(model_data = model_data, m = m,
                   features=config$features[row],
                   indication_method=config$indication_method[row])
  
  logger::log_info("Starting fit")
  mod <- fit(data = xy$train, 
             algorithm=config$algorithm[row],
             indication_method=config$indication_method[row],
             univariate_screen=config$univariate_screen[row])
  logger::log_info("Finished fit")
  
  logger::log_info("Starting predict")
  predictions <- predict_surv(object = mod, 
                              newdata = xy$test,
                              s = "lambda.min")
  logger::log_info("Finished predict")
  
  # We need to run a basic separate model to allow the 'factor' models to evaluate on the same time points t
  if(config$indication_method[row]=="factor") {
    ref_xy <- preprocess(model_data = model_data, m = m,
                         features="benchmark_1",
                         indication_method="separate")
    
    ref_mod <- fit(data = ref_xy$train, 
                   algorithm="coxph",
                   indication_method="separate",
                   univariate_screen=FALSE)
    
    ref_predictions <- predict_surv(object = ref_mod, 
                                    newdata = ref_xy$test,
                                    s = "lambda.min")
    
    reference_object <- ref_predictions
  } else{
    reference_object <- NULL
  }
  
  concordance <- cindex(object = predictions, newdata = xy$test) 
  concordance_var <- cindex_var(object = predictions, newdata = xy$test)
  deviance <- deviance(object = predictions, newdata = xy$test)
  
  calibration <- calibrate(object = predictions, newdata = xy$test)
  
  ibs <- ibs(object = predictions,
             newdata = xy$test,
             reference_object = if(is.null(reference_object)) NULL else reference_object %>% filter(purrr::map_int(surv,nrow) > 1))
  
  isSep <- config$indication_method[row]=="separate"
  
  if (isSep) {
    warnings <- lapply(mod, function(x) attr(x, "warning"))
  } else {
    warnings <- attr(mod,"warning")
  }
  
  if(return_fits){
    return(
      tibble(
        scenario = row,
        mice_m = m,
        xy_train = list(xy$train),
        xy_test = list(xy$test),
        fit = list(mod),
        prediction = list(predictions),
        concordance = list(concordance),
        concordance_var = list(concordance_var),
        deviance = list(deviance),
        calibration = list(calibration),
        ibs = list(ibs),
        config = config_json,
        warnings = list(warnings)
      )
    )
  } else {
    return(
      tibble(
        scenario = row,
        mice_m = m,
        concordance = list(concordance),
        concordance_var = list(concordance_var),
        deviance = list(deviance),
        calibration = list(calibration),
        ibs = list(ibs),
        config = config_json,
        warnings = list(warnings)
      )
    )
  }
  
}