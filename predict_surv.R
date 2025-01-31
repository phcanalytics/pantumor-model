#' @importFrom dplyr as_tibble
#' @export
dplyr::as_tibble

# Generic function for linear predictors ---------------------------------------
predict_lp <- function(object, ...) {
  UseMethod("predict_lp", object)
}

predict_lp._coxph <- function(object, newdata, ...) {
  newx <- make_coxph_newx(newdata, object)
  stats::predict(object$fit, newdata = newx)
}

predict_lp._cv.glmnet <- function(object, newdata, s, ...) {
  stats::predict(object$fit, newx = newdata$x, s = s)[, 1]
}

#' @importFrom foreach %do%
predict_lp.list <- function(object, newdata, ...) {
  lp <- foreach::foreach(i = 1:length(object), .errorhandling = "pass") %do%
    {
      predict_lp(object = object[[i]], newdata = newdata[[i]], ...)
    }
  # If any failed, fill with NA
  err_loc <- purrr::map_lgl(lp, ~is(.,"error"))
  lp[err_loc] <- purrr::map(newdata[err_loc], ~rep(NA_real_, nrow(.$x)))
  return(unname(unlist(lp)))
}

# Generic function for survival distributions ----------------------------------
#' Coerce survival object to data frame
#' 
#' Coerce a `survival::survfit.object` into a tidy `dplyr::tibble`.
#' @param x A `survival::survfit.object` object.
#' @param row_id Optional row ID for each subject that a survival curve is
#' generated for. If `NULL`, then rows IDs are integers ordered from 1 to
#' the number of subjects.
#' @export
as_tibble.survfit <-  function(x, row_id = NULL) {
  
  if (length(x$surv) == 0) {
    stop("survfit() generated a survival curve of length 0.")
  }
  
  if (is.matrix(x$surv)) { # Non-stratified model
    # In this case surv is a matrix. We therefore need to flatten survival
    # curves into a vector and align the times with the flattened curves
    n_obs <- ncol(x$surv)
    n_times <- length(x$time)
    if (is.null(row_id)) row_id <- 1:n_obs
    
    x$surv <- c(x$surv)
    x$time <- rep(x$time, n_obs)
    row_id <- rep(row_id, each = n_times)
  } else { # Stratified model
    # In this case surv and times are already in "long" format, so we just 
    # need to get the IDs right
    n_obs <- length(x$strata)
    row_id <- rep(1:n_obs, times = x$strata)
  }
  
  # Make table
  surv_tbl <- dplyr::tibble(
    row_id = row_id,
    time = x$time,
    surv = x$surv
  )
}

predict_survfit <- function(object, ...) {
  UseMethod("predict_survfit", object)
}

predict_survfit._coxph <- function(object, newdata, ...) {
  newx <- make_coxph_newx(newdata, object)
  survival::survfit(object$fit, newdata = newx,
                    se.fit = FALSE)
}

predict_survfit._cv.glmnet <- function(object, newdata, ...) {
  if (!object$indication_method == "stratify") {
    newstrata <- NULL
  } else{
    newstrata <- newdata$indication
  }
  x <- survival::survfit(object$fit, newx = newdata$x, x = object$x, y = object$y,
                         se.fit = FALSE, newstrata = newstrata)
}

predict_survdistr <- function(object, ...) {
  UseMethod("predict_survdistr", object)
}

predict_survdistr.default <- function(object, newdata, ...) {
  surv <- predict_survfit(object, newdata = newdata, ...)
  as_tibble(surv)
}

#' @importFrom foreach %do%
predict_survdistr.list <- function(object, newdata, ...) {
  Ns <- purrr::map_int(newdata, ~nrow(.$x))
  ids_start <- head(cumsum(c(1,Ns)),n = -1)
  ids_end <- cumsum(Ns)
  surv <- foreach::foreach(i = 1:length(object),
                           id_start = ids_start,
                           id_end = ids_end,
                           .errorhandling = "pass") %do% {
                             survobj <- predict_survfit(object[[i]], newdata = newdata[[i]], ...)
                             as_tibble(survobj, row_id = id_start:id_end)
                           }
  errs <- purrr::map_lgl(surv, ~is(.,"error"))                         
  for (failed_idx in which(errs)) {
    logger::log_error(sprintf("Prediction failed for %s with message %s; relevant ids filled with NA time and surv",
                              names(object)[failed_idx],
                              as.character(surv[[failed_idx]])))
    surv[[failed_idx]] = tibble(row_id = ids_start[failed_idx] : ids_end[failed_idx],
                                time = NA_real_,
                                surv = NA_real_)
  }
  #surv <- surv[!errs]
  dplyr::bind_rows(surv)
}

# Data organization functions --------------------------------------------------
get_predict_id <- function(newdata) {
  if (!inherits(newdata, "list")) newdata <- list(newdata)
  return(tibble(
    patientid = unlist(lapply(newdata, function (z) z$patientid)),
    indication = unlist(lapply(newdata, function (z) z$indication))
  ))
}


predict_surv_tibble <- function(id, lp, surv) {
  # List of survival objects
  surv_list <- dplyr::group_split(surv, row_id)
  
  # Return
  dplyr::tibble(
    id,
    lp = lp,
    surv = surv_list
  )
}

make_coxph_newx <- function(newdata, fit) {
  
  newx <- as.data.frame(newdata$x, check.names = FALSE)
  if (fit$indication_method == "stratify") {
    newx <- data.frame(newx, indication = newdata$indication, check.names = FALSE)
  } 
  return(newx)
}

# Main function: generic predictions from survival models  ---------------------
#' Predictions from survival models
#' 
#' Generate predictions from survival models fit using [fit()] in a standardized format.
#' @param object A survival model fit using [`fit()`] for which prediction is desired.
#' @param newdata Data used for prediction. Must be an object of class 
#' [`survival_xy`].
#' @param s The value of the penalty parameter. Passed to [glmnet::predict.cv.glmnet()].
#' @return An object of class `surv_prediction` that inherits from [dplyr::tibble]. 
#' It contains the following columns:
#' \describe{
#' \item{patientid}{The patient ID.}
#' \item{indication}{The cancer indication.}
#' \item{lp}{The linear predictor.}
#' \item{surv}{A list of tibbles where each row is itself a tibble containing 
#' survival probabilities. The columns in each tibble are `row_id` (the row number
#' that a survival probability corresponds to), `time` (the time at which 
#' a survival prediction is made) and `surv` (the predicted survival probability).}
#' }
#' @export 
predict_surv <- function(object, newdata, s = "lambda.min") {
  if (inherits(object, "list") & !inherits(newdata, "list")) {
    stop("If 'object' is a list, then 'newdata' must be a list.")
  }
  
  lp <- tryCatch(
    suppressWarnings(withCallingHandlers(
      predict_lp(object, newdata = newdata, s = s),
      error = function(e) {
        logger::log_error("predict_lp: {e$message}")
      },
      warning = function(w) {
        logger::log_warn("predict_lp: {w$message}")
      }
    )), error = function(e) e
  )
  if (is(lp, "error")){
    return(lp)
  }
  
  surv <- tryCatch(
    suppressWarnings(withCallingHandlers(
      predict_survdistr(object, newdata = newdata),
      error = function(e) {
        logger::log_error("predict_survdistr: {e$message}")
      },
      warning = function(w) {
        logger::log_warn("predict_survdistr: {w$message}")
      }
    )), error = function(e) e
  )
  if (is(surv, "error")){
    return(surv)
  }
  out <- predict_surv_tibble(id = get_predict_id(newdata),
                             lp = lp, surv = surv)
  class(out) <- c("predict_surv", class(out))
  return(out)
}


# Convenient post-processing of predictions  -----------------------------------
#' Bind survival predictions
#' 
#' Bind a list of survival curves generated from [`predict_surv`] into a single 
#'  object. `rbindlist_surv()` returns a `data.table` object and `bind_surv()`
#'  returns a [`dplyr::tibble`]. 
#' @param object An object of class [`predict_surv`].
#' @export
rbindlist_surv <- function(object) {
  surv <- data.table::rbindlist(object$surv)
  object <- object %>%
    mutate(row_id = unlist(lapply(object$surv, function(x) unique(x$row_id)))) # max(row_id) may not = nrow() if an indication has been removed from the middle of the data
  surv$indication <- object$indication[match(surv$row_id, object$row_id)]
  return(surv)
}

#' @export
#' @rdname rbindlist_surv
bind_surv <- function(object) {
  dplyr::as_tibble(rbindlist_surv(object))
}