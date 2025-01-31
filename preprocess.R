standardize <- function(x, center = NULL, scale = NULL){
  mx <- if(is.null(center)) colMeans(x) else center
  sx <- if(is.null(scale)) sqrt(apply(x, 2, stats::var)) else scale
  x <- scale(x, mx, sx)
  #x[, sx == 0] <- 0 # We should have already done this
  return(x)
}

get_zero_variance <- function(train){
  
  # Store zero variance numeric predictors
  zero_var <- 
    train %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select_if(~var(., na.rm = TRUE) == 0) %>%
    colnames # store just the predictor names
  
  # Store strictly zero variance factor predictors
  zero_var_factors <- 
    train %>%
    dplyr::select(where(is.factor) | where(is.character)) %>%
    dplyr::select_if(~(length(unique(.[!is.na(.)])) == 1L) | all(is.na(.))) %>% # potentially remove the all(is.na(.)) part once the "assume zero" imputation is fixed
    colnames # store just the predictor names
  
  # Store zero variance biomarkers- biomarkers with var == 0 when rounded to nearest whole number
  zero_var_biomarkers <-
    train %>% 
    dplyr::select(contains("genomics_binary")) %>%
    dplyr::mutate_all(function(x) round(x, 0)) %>%
    dplyr::select_if(~var(., na.rm = TRUE) == 0) %>%
    colnames # store just the predictor names
  
  # Return zero variance & low-frequency/high class imbalance predictors
  return(c(zero_var, zero_var_factors, zero_var_biomarkers)) #, low_prop, low_prop_bin))
}

get_uncommon_dummy <- function(data) {
  # Select binary variables and look for rare counts
  data_binary <-  data %>%
    dplyr::select_if(function(x) is.binary(x) | length(unique(x))==1L)
  low_0 <- data_binary[,which(colSums(data_binary)>=nrow(data_binary)/2),drop=FALSE] %>% select_if(apply(., 2, function(x) sum(x==0)<20)) %>% colnames
  low_1 <- data_binary[,which(colSums(data_binary)<nrow(data_binary)/2),drop=FALSE] %>% select_if(apply(., 2, function(x) sum(x==1)<20)) %>% colnames
  
  # Here we will also get uncommon biomarkers imputed on continuous scale
  data_biomarkers <- data %>% 
    dplyr::select(contains("genomics_binary")) %>%
    dplyr::mutate_all(function(x) round(x, 0)) %>%
    dplyr::select_if(function(x) is.binary(x) | length(unique(x))==1L)
  low_0_bm <- data_biomarkers[,which(colSums(data_biomarkers)>=nrow(data_biomarkers)/2),drop=FALSE] %>% select_if(apply(., 2, function(x) sum(x==0)<20)) %>% colnames
  low_1_bm <- data_biomarkers[,which(colSums(data_biomarkers)<nrow(data_biomarkers)/2),drop=FALSE] %>% select_if(apply(., 2, function(x) sum(x==1)<20)) %>% colnames
  
  return(c(low_0, low_1, low_0_bm, low_1_bm))
}

select_features <- function(df, 
                            features = c("benchmark_1", "ropro_2","full_3"),...){
  
  # Define outcomes data  
  outcomes <- c("entry_days_l1", "os_days_l1", "died")
  
  if (features =="benchmark_1"){
    model_features = c("age", "race", "gender", "ml_smoking_status","time_dx_to_l1", "entry_days_l1_predictor", "baseline_l1_ecog", "groupstage") 
    cols <- c(model_features)
    data <- df[, names(df) %in% unique(cols),drop=FALSE]
  }  else if (features =="ropro_2"){
    model_features = c("age", "gender", "ml_smoking_status", "baseline_l1_ecog", "groupstage",
                       "last_value_body_weight", 
                       "last_value_body_height", 
                       "last_value_heart_rate", 
                       "last_value_hemoglobin_whole_blood", 
                       "last_value_systolic_blood_pressure", 
                       "last_value_diastolic_blood_pressure",
                       "last_value_urea_nitrogen", 
                       "last_value_alkaline_phosphatase_alp", 
                       "last_value_alanine_aminotransferase_alt_or_sgpt",
                       "last_value_aspartate_aminotransferase_ast_or_sgot",
                       "last_value_calcium_serum", 
                       "last_value_creatinine_serum", 
                       "last_value_proteintotal_serum", 
                       "last_value_bilirubin_total_serum", 
                       "last_value_albumin_serum",
                       "last_value_hematocrit", 
                       "last_value_glucose", 
                       "last_value_platelet_count",
                       "last_value_lymphocyte_count_absolute", 
                       "last_value_monocyte_count", 
                       "last_value_neutrophil_count_absolute",
                       "BMI_cat", 
                       "entry_days_l1_predictor")
    cols <- c(model_features) 
    data <- df[, names(df) %in% unique(cols),drop=FALSE] 
  } else {
    # For the full model, take all variables, but exclude outcomes and other metadata
    data = df %>%
      dplyr::select(-outcomes) %>%
      dplyr::select(-mice.imp)%>%
      dplyr::select(-indication)%>% 
      dplyr::select(-patientid) 
  }
  # Save the feature names for future reference
  return(tolower(names(data)))
}




# survival_xy object -----------------------------------------------------------
#' Survival input and response data
#' 
#' Preprocess a data frame to create a `survival_xy` object, 
#' which stores an input matrix, a survival object storing response variables in 
#' a survival model, and other data attributes. 
#' @param data A data frame containing the variables for modeling. May either
#' be the training or the test data.
#' @param m The multiple imputation dataset to select for preprocessing. This should be 
#' an integer between 1 and 5 stored in the column `mice.imp`.
#' @param object A `survival_xy` object. This is useful if `data` is the test data
#' and there is a need to incorporate information from training data when 
#' preprocessing the test dataset. 
#' @param features A character vector containing the names of the features in `data`
#' to include.
#' 
#' @return A `survival_xy` object which is a list containing:
#' \describe{
#' \item{x}{An input matrix where each row is an observation and each column is
#' a model term.}
#' \item{y}{A [survival::Surv()] object containing the response variables 
#' for survival modeling.}
#' \item{m}{The multiple imputation dataset number.}
#' \item{patientid}{The patient IDs corresponding to rows in `x`.}
#' \item{indication}{The cancer indication corresponding to rows in `x`.}
#' \item{fold_id}{If the object contains training data, the cross-validation fold
#' in which rows in `x` are in the validation set.}
#' \item{features}{Character vector of the features using in `x`. These are not
#' the same as the columns in `x` because categorical variables are treated as
#' a single feature, rather than as separate dummy variables as in `x`.}
#' }
#' @export
survival_xy <- function(data, object = NULL, features) {
  
  
  # (1) Create "y" output
  y <- survival::Surv(data$entry_days_l1,
                      data$os_days_l1,
                      data$died)
  
  # (2) Create input "x" matrix
  ## (2a) Subset columns
  ### If test data, then just use the same features as train
  if (!is.null(object)) {
    x <- data[, names(data) %in% object$features, drop=FALSE]
  } else { ### Otherwise...
    ### First subset to relevant columns
    x <- data[, names(data) %in% features, drop=FALSE]
    
    ### Then remove features with zero variance 
    zero_vars <- get_zero_variance(x)
    features <- setdiff(features, zero_vars)
    x <- x[, names(x) %in% features, drop=FALSE]
    
  } # End else for feature subsetting
  
  ## (2b) Convert categorical variables to dummies via one-hot encoding
  if (any(sapply(x, is.factor))) {
    if (is.null(object)) { # training data
      x <- fastDummies::dummy_cols(x,
                                   remove_most_frequent_dummy = TRUE, # remove reference level (most frequent)
                                   ignore_na = TRUE, # assume all NAs have been dealt with
                                   remove_selected_columns = TRUE) # don't keep original data
    } else { # make choice of ref level align with training by NOT dropping any levels here and relying on later step to filter
      x <- fastDummies::dummy_cols(x,
                                   ignore_na = TRUE, # assume all NAs have been dealt with
                                   remove_selected_columns = TRUE) # don't keep original data
      
    }
  }
  
  # While fastDummies was told to not create a separate level for NA entries,
  # it did not fill them with 0 for the levels that entry was not, so we do it now
  
  x <- x %>% mutate(across(where(is.numeric), ~tidyr::replace_na(.,0)))
  
  ## (2b.2) Now eliminate dummy variables that have counts of the least frequent level < 20  
  # ONLY IF it's the train set! Otherwise, you might remove variables that are used in train but are zero variance in test
  if (is.null(object)) {
    zero_vars <- get_uncommon_dummy(x)
    z <- setdiff(names(x), zero_vars)
    x <- x[, names(x) %in% z, drop=FALSE]
  } else { # Now make sure that for the test set, we subset only to the dummy variables retained in the train set
    x <- x[, names(x) %in% colnames(object$x), drop=FALSE]
  }
  
  ##  (2c) Standardize all variables
  if (!is.null(object)) { ### Use training information if "data" is the test object
    x_scale <- attr(object$x, "scaled:scale")
    x_center <- attr(object$x, "scaled:center")
  } else{
    x_scale <- x_center <- NULL # Done by default in standardize() function
  }
  x <- standardize(as.matrix(x), center = x_center, scale = x_scale)
  
  ## (2d) Truncate outliers +/- 3 z-scores away
  truncate <- function(x){
    min <- min( x[x>=-3])
    max <- max(x[x<=3])
    x[x <= -3] <- min
    x[x >= 3] <- max
    x
  }
  
  # ecog and groupstage are exceptions- treated as continuous but result in high z-values, plan to keep all
  ecog <- c(as.vector(which(colnames(x)=="baseline_l1_ecog")))
  gpstg <- c(as.vector(which(colnames(x)=="groupstage")))
  binary_colnums <- c(as.vector(which(apply(x, 2, is.binary))),ecog,gpstg)
  nonbinary_colnums <- c(as.vector(which((apply(x, 2, is.binary)==FALSE))))[-c(which(c(as.vector(which((apply(x, 2, is.binary)==FALSE))))==ecog),which(c(as.vector(which((apply(x, 2, is.binary)==FALSE))))==gpstg))]
  
  nonbinary_trunc <- apply(x[,nonbinary_colnums, drop=FALSE], 2, truncate) %>%
    {
      if (nrow(x) == 1){
        t(.)
      }
      else {
        (.)
      }
    }
  
  x_trunc <- cbind(x[,binary_colnums, drop=FALSE], # do not apply the truncation to binary variables!
                   nonbinary_trunc
  )
  colnames(x_trunc) <- c(colnames(x)[binary_colnums, drop=FALSE], 
                         colnames(x)[nonbinary_colnums, drop=FALSE])
  
  x_trunc <- x_trunc[,colnames(x),drop=FALSE]
  
  attributes(x_trunc)$`scaled:center` <- attributes(x)$`scaled:center`
  attributes(x_trunc)$`scaled:scale`  <- attributes(x)$`scaled:scale`
  
  x <- x_trunc
  
  # (3) Check for any NAs
  if(any(is.na(x))) {
    stop("x matrix has NA elements.")
  }
  
  # (4) Return
  fold_id <- if (!is.null(data[["fold_id"]])) data[["fold_id"]] else NULL
  res <- list(x = x,
              y = y,
              indication = as.character(data$indication),
              patientid = data$patientid,
              fold_id = fold_id,
              features = features)
  class(res) <- "survival_xy"
  return(res)
}

#' Get response data
#' 
#' Get the `$y` element of a `survival_xy` object or a list of `survival_xy`
#' objects. If a list, then the elements of `y ` are combined into a single 
#' [`survival::Surv`] object.
#' @param object An object of class `survival_xy` or a list of `survival_xy`
#' objects.
#' @return A [`survival::Surv`] object  In cases where `object` is a list,
#' then attributes `indication` and `patientid` are added, which are vectors 
#' denoting the indication and patient ID associated each row of the `Surv` 
#' object.
#' @export
get_y <- function(object) {
  # (1) First case is we have one pantumor model. Then we just return $y
  if(inherits(object, "survival_xy")) return(object$y)
  
  # (2) Second case is we have a list - separate cancer models (1 element per cancer)
  ## First check the elements of the list are correct
  if (!inherits(object, "list")) {
    stop("'object' must either be a 'survival_xy' object or a list.")
  }
  
  ## We now combine the elements of the list into a single `Surv` object
  y_df <- lapply(object, function (z) {
    if (!inherits(z, "survival_xy")) {
      stop("If 'object' is a list, each element must be of class 'survival_xy'.")
    }
    dplyr::tibble(indication = z$indication,
                  patientid = z$patientid,
                  start = z$y[, "start"],
                  stop = z$y[, "stop"],
                  status = z$y[, "status"])
  }) %>%
    dplyr::bind_rows()
  y <- survival::Surv(time = y_df$start, 
                      time2 = y_df$stop,
                      event = y_df$status)
  attr(y, "patientid") <- y_df$patientid
  attr(y, "indication") <- y_df$indication
  return(y)
}

# Main preprocess() function ---------------------------------------------------
#' Preprocess data
#' 
#' Preprocess the data required to train and test the survival models into a format
#' that is suitable for modeling.
#'
#' @inheritParams fit
#' @param model_data An object of class `model_data` containing: 
#' `data`: a data frame containing patients in both the training and test set;
#' It should contain at least the columns `patientid`, `indication`, `mice.imp`.
#' `train_set`: a data frame containing information on the training set.
#' It should contain the columns `patientid`, `indication`, `mice.imp`,
#' and `fold_id`, where `fold_id` corresponds to use in [cv.glmnet()`] and 
#' is the fold in which a given row is used as validation data;
#' `test_set`: a data frame containing information on the test set.
#' It should contain the columns `patientid` and `indication`.
#' @param m The multiple imputation dataset to select for preprocessing. This should be 
#' an integer between 1 and 5 stored in the column `mice.imp`.
#' @param features A character vector containing the names of the features in `data`
#' to include.
#'
#' @return If a single model is going to be fit across all indications 
#' (`indication != "separate"`), then an object of class `survival_xy`, which is
#' a list with the following elements:
#' \describe{
#' \item{x}{The input matrix for survival modeling.}
#' \item{y}{A [`survival::Surv`] object used as the response variable when modeling.}
#' }
#' In cases where separate models are fit for each indication (`indication == "separate"`),
#' then a list of `survival_xy` objects is returned.
#' @export
preprocess <- function(model_data, 
                       m,
                       features, 
                       indication_method = c("factor", "stratify", "separate")){
  
  if (!inherits(model_data, "model_data")){
    stop("'model_data' must be an object of class 'model_data'. See 'make_model_data()'.")
  }
  
  data <- model_data$data
  
  features = select_features(data, features)
  
  indication_method <- match.arg(indication_method)
  if (indication_method == "factor") features <- c(features, "indication")
  
  # Use all lowercase names
  colnames(data) <- tolower(colnames(data))
  
  # Filter to multiple imputation dataset
  data <- data %>%
    dplyr::filter(mice.imp==m)
  
  # Keep columns we need
  outcomes <- c("entry_days_l1", "os_days_l1", "died")
  cols <- c(features, outcomes, "patientid", "indication")
  data <- data[, names(data) %in% unique(cols),drop=FALSE]
  
  # Death variable should be numeric and entry days should be non-negative, and NAs should be 0
  data <- data %>%
    dplyr::mutate(
      entry_days_l1 = pmax(entry_days_l1, 0), # Set negative entry days to zero
      died = as.integer(died)# Death indicator should be numeric 
    ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric),~tidyr::replace_na(.,0))) # Replace NAs with 0 - for block missing variables (variables which have been imputed for their relevant cancer subgroups and are 0 for all else)
  
  
  # Convert character variables (other than ID variables) to factors. Note that
  # indication will be converted to a factor below.
  character_cols <- names(which(sapply(
    data[, !names(data) %in% c("patientid", "indication"),drop=FALSE], 
    is.character
  )))
  data <- data %>% 
    dplyr::mutate( dplyr::across(any_of(character_cols), as.factor))
  
  # For factor variables, collapse small levels into a single category. That way they are not eliminated in the uncommon dummy step
  data <- data %>%
    dplyr::mutate_if(is.factor, ~forcats::fct_lump_min(., min = 20))
  
  # If separate models are required for each indication, then we need to convert
  # our data to a list with one element for each indication. In cases where we
  # have one pantumor model, we still use a list (of length 1) so that our data
  # always has a standard format
  make_data_list <- function(data, set, indication_method) {
    
    # Subset data training/test set patients
    colnames(set) <- tolower(colnames(set))
    d <- data %>%
      dplyr::inner_join(set, by = c("patientid", "indication")) %>%
      dplyr::mutate(indication = as.factor(indication)) # Needs to be a factor if used in "x" matrix
      
      # Now create a list of datasets. The list is of length equal to the number
      # of indications if indication_method == "separate" and of length 1 otherwise
      if (indication_method != "separate") {
        d_list <- list(d)
      } else {
        d_list <- dplyr::group_split(d, by = indication)
      }
    return(d_list)
  }
  train_list <- make_data_list(data, model_data$train_set, indication_method)
  test_list <- make_data_list(data, model_data$test_set, indication_method)
  
  # Create a list of survival_xy objects
  map_survival_xy <- function(data_list, object_list = NULL) {
    res <- vector(mode = "list", length = length(data_list))
    for (i in 1:length(res)) {
      res[[i]] <- survival_xy(data_list[[i]], object_list[[i]],
                              features = features)
      names(res)[i] <- res[[i]]$indication[1]
    }
    return(res)
  }
  xy_train <- map_survival_xy(train_list)
  xy_test <- map_survival_xy(test_list, xy_train)
  
  # Return
  return(list(
    train = if (length(xy_train) == 1) xy_train[[1]] else xy_train,
    test = if (length(xy_test) == 1) xy_test[[1]] else xy_test
  ))
}
