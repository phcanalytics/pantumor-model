# Evaluate statistic given predictions on a test set ---------------------------
#' Evaluate statistic
#' 
#' This is a general function to evaluate a statistic using predictions from a survival model.
#' It is used internally within the functions [cindex()] and [deviance.predict_surv()].
#' @param object An object of class [`predict_surv`], which contains predictions made on a test set.
#' @param newdata An object of class [`survival_xy`] (or a list of objects of class [`survival_xy`]) 
#' containing new values of `$y` at which predictions should be evaluated. 
#' @param f A function that evaluates a statistic given the values of `object`
#' and `newy`, where `newy` is the result of calling [get_y()] on `newdata`. `f` 
#' must return a single value.
#' 
#' @return A [`dplyr::tibble`] with two columns: `indication` for the indication and 
#' `estimate` for an estimate computed using `f`. There is also a "Pan tumor" estimate
#' which evaluates the statistic across all indications.
#' @export
evaluate_statistic <- function(object, newdata, f, reference_object=NULL) {
  
  newy <- get_y(newdata)
  
  # Pan tumor model
  pan_estimate <- do.call(f, args = list(newy = newy, object = object))
  names(pan_estimate) <- "Pan tumor"
  
  # Indication specific models
  object_list <- split(object[, ], object$indication)
  newy_list <- split(newy, object$indication)
  if(!is.null(reference_object)){
    ref_list <- split(reference_object, reference_object$indication)
    shared_indications <- intersect(names(object_list), names(ref_list))
    
    indication_estimates <- mapply(function(x, y, z) {
      do.call(f, args = list(object = x, newy = y, reference_object = z))
    }, object_list[shared_indications], newy_list[shared_indications], ref_list[shared_indications]
    )
  } else{
    indication_estimates <- mapply(function(x, y) {
      do.call(f, args = list(object = x, newy = y))
    }, object_list, newy_list
    )
  }
  
  estimates <- c(pan_estimate, indication_estimates)
  return(tibble(
    indication = names(estimates),
    estimate = estimates
  ))
}

# C-index ----------------------------------------------------------------------
#' C-index
#'
#' Evaluate concordance from a survival model.
#'
#' @inheritParams evaluate_statistic
#' @param ... Additional arguments affecting the evaluation of concordance. Currently unused. 
#'
#' @return A [`dplyr::tibble`] with two columns: `indication` for the indication and 
#' `estimate` for the estimate of the C-index.
#' @export
cindex <- function(object, newdata, ...){
  UseMethod("cindex", object)
}

#' @rdname cindex
#' @export
cindex.predict_surv <- function(object, newdata, ...){
  cindex_fun <- function(object, newy) {
    if (nrow(newy) > 1) {
      survival::concordance(newy ~ object$lp, reverse = TRUE)$concordance
    } else {
      NA_real_
    }
  }
  evaluate_statistic(object, newdata = newdata, f = cindex_fun)
}

# C-index variance ----------------------------------------------------------------------
#' C-index variance
#'
#' Obtain variance of concordance estimate from a survival model.
#'
#' @inheritParams evaluate_statistic
#' @param ... Additional arguments affecting the evaluation of concordance. Currently unused. 
#'
#' @return A [`dplyr::tibble`] with two columns: `indication` for the indication and 
#' `estimate` for the estimate of the C-index variance.
#' @export
cindex_var <- function(object, newdata, ...){
  UseMethod("cindex_var", object)
}

#' @rdname cindex_var
#' @export
cindex_var.predict_surv <- function(object, newdata, ...){
  cindex_var_fun <- function(object, newy) {
    if (nrow(newy) > 1) {
      survival::concordance(newy ~ object$lp, reverse = TRUE)$var
    } else {
      NA_real_
    }
  }
  evaluate_statistic(object, newdata = newdata, f = cindex_var_fun)
}

# Deviance ---------------------------------------------------------------------
#' Deviance
#' 
#' Evaluate deviance from a survival model.
#'
#' @inheritParams cindex
#' @return A [`dplyr::tibble`] with two columns: `indication` for the indication and 
#' `estimate` for the estimate of the deviance.
#' @export
deviance.predict_surv <- function(object, newdata, ...){
  deviance_fun <- function(object, newy) {
    glmnet::coxnet.deviance(object$lp, newy)
  }
  evaluate_statistic(object, newdata = newdata, f = deviance_fun)
} 


# Integrated Brier Score --------------------------------------------------
#' Integrated Brier Score
#' 
#' Estimate Integrated Brier Score from a survival model.
#' The IBS is scaled by max_t in the data.
#' Uses left truncation adjustment from `LTRCforests`.
#'
#' @inheritParams evaluate_statistic
#' @return A [`dplyr::tibble`] with two columns: `indication` for the indication and 
#' `estimate` for the estimate of the Integrated Brier Score, scaled by max_t.
#' @export
ibs <- function(object, newdata, ...){
  UseMethod("ibs", object)
}

#' @rdname ibs
#' @export

ibs.predict_surv <- function(object, newdata, reference_object=NULL, ...){
  
  ibs_fun <- function(object, newy, reference_object=NULL){
    
    newy_mat = as.matrix(newy)
    df = rbindlist_surv(object) 
    names(df)[3] <- "pred"
    n_times <- df[, .N, by = c("row_id")] # Account for potential differences in followup 
    df <- cbind(df, newy_mat[rep(1:nrow(newy_mat), times = n_times$N), ])
    
    # Define the event times at which the IBS score will be computed
    # Use the time points from referenceobject
    if(!is.null(reference_object)){
      if(length(unique(df$indication))>1){
        times <- unname(stats::quantile(df$time, probs = seq(0, 1, by=0.01), na.rm = TRUE))
      } else {
        times <- sort(unique(reference_object$surv[[1]]$time))
      }
    } else{
      if(length(unique(df$indication))>1){
        times <- unname(stats::quantile(df$time, probs = seq(0, 1, by=0.01), na.rm = TRUE))
      } else {
        times <- sort(unique(df$time))
      }
    }
    
    # Interpolates survival probabilities for a common set of time points, for all patients and all indications.
    req_times = tidyr::expand_grid(row_id = unique(df$row_id), time = times)
    data.table::setDT(req_times)
    data.table::setkey(df, row_id, time)
    df <- df[req_times,roll=Inf, rollends=c(T,T)]
    
    # Now define the inputs to "pred": probabilities, evaluation times, max time, and n
    probs <- lapply(split(df , f = df$row_id ), function(x) as.vector(unlist(x[,3])))
    unique_times = times # this returns the unique vector of times for the indication or pan-tumor
    tau = unname(quantile(df$time, probs=0.9, na.rm = TRUE)) # set upper limit of evaluation times = 90th percentile
    n = length(unique(df$row_id))
    
    pred = list(
      survival.probs = lapply(probs, function(x) c(1,x)), # start at t=0, prob=1
      survival.times = c(0,as.vector(unique_times)), # start at t=0, prob=1
      survival.tau = rep(tau, n),
      survival.obj = newy,
      survival.id = 1:n
    )
    
    ibs <- unname(LTRCforests::sbrier_ltrc(obj = newy, id= 1:n, pred = pred, type = "IBS"))
    ibs # returns the average IBS score = IBS / tau, which is ALREADY computed by the function
  }
  
  evaluate_statistic(object, newdata = newdata, f = ibs_fun, reference_object = reference_object)
}

# Calibration ------------------------------------------------------------------
#' Subset survival predictions
#' 
#' A predicted survival curves at specific times, which may not times at
#' which the original curves were evaluated. The step-function nature of the
#' predictions are preserved so that if `times` lies within `[t1, t2)`, the 
#' survival curves evaluated at `t1` will be used. An efficient approach using
#' `data.table` is employed.
#' @param x Predicted survival curves as generated by [predict_surv()]. 
#' @param times Times at which to subset survival curves to.
#' @param extend Whether to extend survival predictions beyond maximum followup
#' time.
subset_surv <- function(x, times, extend = FALSE) {
  x <- data.table::as.data.table(x)
  time_index <- findInterval(x$time, times, left.open = TRUE) + 1L 
  x$u <- times[time_index]
  x2 <- x[x[, .I[.N], by = c("row_id", "u")]$V1] # Get last row for each u by ID
  
  # Remove "times" if they are beyond what is observed in the survival curves
  # from "x"
  if (extend == FALSE) {
    x2[, max_time := max(time), by = "indication"]
    x2 <- x2[u <= max_time]
  }
  
  # Return
  return(x2)
}


#' Calibrate survival model
#'
#' Calibrate survival models by comparing predicted survival
#' probabilities with pseudo-observed probabilities at selected time points.
#' 
#' @inheritParams evaluate_statistic
#' @param n_groups Number of groups to break predicted survival probabilities
#' into. 
#' @param times A list of selected time points for which the calibration should be calculated.
#'
#' @return An object of class `calibrate` that inherits from [dplyr::tibble] 
#' comparing predicted survival probabilities to pseudo-observed probabilities at 
#' specified time points. The table contains the following columns:
#' \describe{
#' \item{time}{The time point at which a survival probability is evaluated.}
#' \item{interval}{An identifier of a group whose predicted survival probabilities lie 
#' within a specific interval. The number of intervals is determined by
#'  `n_groups`.}
#' \item{pred}{The average predicted survival probability among patients in an 
#' interval.}
#' \item{n_patients}{The number of patients withn an interval.}
#' \item{obs}{The observed survival probability among patients in an interval. 
#' Computed with the Kaplan-Meier estimator using the `$y` element of `newdata`.}
#' \item{indication}{The cancer indication.}
#' 
#' }
#' @importFrom data.table ':='
#' @export
calibrate <- function(object, newdata, n_groups = 10,
                      times = NULL){
  
  # (1) Combine the predicted survival probabilities with the observed
  ## Row bind predicted survival probabilities
  if (!inherits(object, "predict_surv")) {
    stop("'object' must be of class 'predict_surv'.")
  }
  surv <- rbindlist_surv(object)
  
  ## Restrict number of times 
  if(is.null(times)) {
    times <- stats::quantile(surv$time, probs = seq(0, 1, length.out = 7)[-1], na.rm = TRUE)
  }
  surv <- subset_surv(surv, times)
  
  ## Add Observed survival outcomes
  ## Note that if baseline hazards vary by indication, then not all rows of y will
  ## contain all times because of differences in followup. 
  newy <- get_y(newdata);  ymat <- as.matrix(newy)
  data.table::setnames(surv, c("surv"), c("pred"))
  n_times <- surv[, .N, by = c("row_id")] # Account for potential differences in followup 
  surv <- cbind(surv,
                ymat[rep(1:nrow(ymat), times = n_times$N), ])
  
  # (2) Calibrate
  ## Calibrate pan tumor model
  cal_pan <- calibrate_surv(surv, newy = newy, n_groups = n_groups,
                            times = times)
  cal_pan[, indication := "Pan tumor"]
  
  ## Indication specific models
  surv_list <- split(surv, by = "indication")
  newy_list <- split(newy, object$indication)
  
  cal_indication <- mapply(function(x, y) {
    calibrate_surv(surv = x, newy = y, n_groups = n_groups,
                   times = times)
  }, surv_list, newy_list, SIMPLIFY = FALSE
  )
  cal_indication <- data.table::rbindlist(cal_indication, idcol = "indication")
  
  ## Combine
  res <- rbind(cal_pan, cal_indication)
  res <- dplyr::as_tibble(res)
  class(res) <- c("calibrate", class(res))
  return(res)
}

calibrate_surv <- function(surv, newy, n_groups = NULL, times = NULL,
                           ...) {
  
  # (1) Cutpoints of predicted S(u|x)
  ntile <- function(x, n){
    cut(x, 
        breaks = quantile(x, 
                          probs = seq(0, 1, length = n + 1), 
                          na.rm = TRUE,
                          type = 2),
        include.lowest = TRUE,
        labels = FALSE)
  }
  if (n_groups > 1){
    #NOTE: ntile() and dplyr::ntile() should be equivalent but ntile() seems to cause more errors
    # when n + 1 quantiles cannot be recreated.
    surv[, interval := dplyr::ntile(pred, n = n_groups), by = "u"] 
  } else{
    surv[, interval := 1]
  }
  
  # (2) Average predicted S(u|x) in each interval
  surv_mean <- surv[, .(pred = mean(pred),
                        n_patients = .N),
                    by = c("u", "interval")]
  
  # (3) Compute KM for patients in each interval
  if(ncol(newy) == 3){
    kmfit <- survival::survfit(survival::Surv(start, stop, status) ~ survival::strata(interval), 
                               data = surv)
  } else{
    kmfit <- survival::survfit(survival::Surv(time, status) ~ survival::strata(interval), 
                               data = surv)
  }
  kmfit_summary <- summary(kmfit, times = times, extend = TRUE)
  if (n_groups > 1){
    strata <- as.integer(gsub("survival::strata(interval)=interval=", "",
                              kmfit_summary$strata, fixed = TRUE))
  } else{
    strata <- 1
  }
  kmfit_df <- data.table::data.table(
    interval = strata,
    time = kmfit_summary$time,
    obs = kmfit_summary$surv
  )
  data.table::setnames(surv_mean, "u", "time")
  surv_mean <- merge(surv_mean, kmfit_df, 
                     by = c("time", "interval"), 
                     all.x = TRUE)
  
  # Return
  return(surv_mean)
}

autoplot.calibrate <- function(object, colour = NULL){
  object$f_time <- factor(object$time,
                          levels = object$time,
                          labels = paste0("Time = ", object$time))
  p <- ggplot(object)
  if (is.null(colour)){
    aes <- aes_string(x = "pred", y = "obs", label = "interval")
  } else{
    aes <- aes_string(x = "pred", y = "obs", col = colour, label = "interval")
  }
  p <- p +
    aes +
    geom_point(size = 2) +
    geom_abline(slope = 1) +
    facet_wrap(~f_time) +
    scale_shape_discrete(name = "Model") +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    xlab("Predicted survival probability") + 
    ylab("Observed survival probability") 
  return(p)
}
