check_survival_xy <- function(data) {
  
  # Check if data is a list
  if (inherits(data, "list")) {
    is_survival_xy <- sapply(data, function (x) inherits(x, "survival_xy"))
    is_survival_xy <- all(is_survival_xy)
  } else if (inherits(data, "survival_xy")) {
    is_survival_xy <- TRUE
  } else {
    is_survival_xy <- FALSE
  }
  
  if (!is_survival_xy) {
    stop(paste0("'data' must be an object of class 'survival_xy' or ", 
                "a list of objects of class 'survival_xy'."))
  }
}

#' Fit a survival model
#' 
#' This is a general function for fitting a survival model using different
#' possible algorithms. 
#' 
#' @param data Data used to fit the model. Must be an object of class 
#' [`survival_xy`].
#' @param algorithm The algorithm used to fit the model. Either an unpenalized 
#' ([survival::coxph()]) or cross-validated penalized ([glmnet::cv.glmnet()]) 
#' Cox proportional hazards model.
#' @param indication_method Whether to account for different tumor indications
#' with a factor (i.e., dummy) variable (`"factor"`), 
#' by stratifying the baseline hazard (`"stratify"`),
#' or by fitting separate models for each indication (`"separate"`).
#' @param univariate_screen Logical. If `TRUE`, then univariate screening is used
#' to determine which features should be included in the model within k-fold CV.
#' @param ... Additional arguments to pass to underlying fitting functions
#' specified in `algorithm`.
#' 
#' @return If a single model is fit, then an object of class `model_fit`, 
#' in addition to a subsequent class related to the fitted model (e.g., `_coxph`).
#' If separate models are fit, then a list of `model_fit` objects.
#' @importFrom survival strata
#' @export
fit <- function(data, algorithm = c("coxph", "cv.glmnet"), 
                indication_method = c("factor", "stratify", "separate"),
                univariate_screen = FALSE,
                ...) {
  
  algorithm <- match.arg(algorithm)
  indication_method <- match.arg(indication_method)
  check_survival_xy(data)
  
  if(!univariate_screen){
    if (algorithm  == "coxph" && indication_method  != "separate") {
      fit <- fit_coxph(data = data, indication_method = indication_method, ...)
    } else if (algorithm  == "cv.glmnet" && indication_method  != "separate") {
      fit <- fit_cv.glmnet(data = data, indication_method = indication_method, ...)
    } else {
      fit <- fit_list(data = data, algorithm = algorithm, 
                      indication_method = indication_method, ...)
    }
    return(fit)
  } else {
    if (algorithm  == "cv.glmnet" && indication_method  != "separate") {
      fit <- fit_cv.glmnet(data = data, indication_method = indication_method, exclude=prescreen_genomic_vars, ...)
    } else if (algorithm  == "cv.glmnet" && indication_method  == "separate") {
      fit <- fit_list(data = data, algorithm = algorithm, 
                      indication_method = indication_method, exclude=prescreen_genomic_vars, ...)
    } else {
      stop("Univariate screen is only implemented when algorithm == cv.glmnet.")
    }
    return(fit)
  }
  
}

fit_coxph <- function(data, indication_method, ...) {
  
  message_postfix <- ifelse(
    "indication" %in% names(attributes(data)),
    paste0(", indication: ",attributes(data)$indication),
    "")
  
  warn <- NULL
  
  fit <- tryCatch(suppressWarnings(withCallingHandlers(
    {
      if (indication_method == "stratify"){
        survival::coxph(data$y ~ . -indication + strata(indication), data = data.frame(data$x, indication = data$indication,
                                                                                       check.names = FALSE))
      } else {
        survival::coxph(data$y ~ ., data = data.frame(data$x, check.names = FALSE))
      }
    },
    warning = function(w) {
      logger::log_warn("coxph: {w$message}{message_postfix}")
      warn <<- w
    },
    error = function(e) {
      logger::log_error("coxph: {e$message}{message_postfix}")
    })),
    error = function(e) e
  )
  
  
  object <- list(fit = fit,
                 indication_method = indication_method)
  class(object) <- c("model_fit", "_coxph")
  if(!is.null(warn)) {
    attr(object, "warning") <- warn
  }
  return(object)
}

fit_cv.glmnet <- function(data, indication_method, ...) {
  
  # Set glmnet control
  glmnet::glmnet.control(eps = 1e-04, mxitnr=10000)
  
  message_postfix <- ifelse(
    "indication" %in% names(attributes(data)),
    paste0(", indication: ",attributes(data)$indication),
    "")
  
  if(indication_method=="stratify"){
    data$y <- glmnet::stratifySurv(data$y, strata = data$indication)
  }
  
  warn <- NULL
  fit <- tryCatch(suppressWarnings(withCallingHandlers(
    {
      glmnet::cv.glmnet(x = data$x, y = data$y, family = "cox", 
                        foldid = data$fold_id, standardize = FALSE, type.measure = "C", ...) # automatically fits stratified model if y is a stratifySurv object
    },
    warning = function(w) {
      logger::log_warn("glmnet::cv.glmnet {w$message}{message_postfix}")
      warn <<- w
    },
    error = function(e) {
      logger::log_error("glmnet::cv.glmnet {e$message}{message_postfix}")
    })),
    error = function(e) e
  )
  object <- list(
    fit = fit,
    x = data$x,
    y = data$y,
    indication_method = indication_method
  )
  if(!is.null(warn)) {
    attr(object, "warning") <- warn
  }
  class(object) <- c("model_fit", "_cv.glmnet")
  return(object)
}

fit_list <- function(data, algorithm, indication_method, ...) {
  fun <- paste0("fit_", algorithm)
  args <- c(list(indication_method = indication_method), 
            list(...))
  #curr_log_format <- logger::log_layout()
  fits <- list()
  for(i in seq_along(data)) {
    x <- data[[i]]
    attributes(x)$indication <- names(data)[[i]]
    #format_with_ind <- paste0(as.character(curr_log_format$format),", indication:",names(data)[[i]])
    #logger::log_layout(format_with_ind)
    args_i <- c(list(data = x), args)
    fit <- do.call(fun, args_i)
    fits[[i]] <- fit
  } 
  names(fits) <- names(data)
  #logger::log_layout(curr_log_format)
  return(fits)
}
