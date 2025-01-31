#' Detect binary variable
#' 
#' Detect a binary variable, for purposes of factor conversion.
#' 
#' @param v The variable of interest.
#' @export

is.binary <- function(v) { 
  x <- unique(v) 
  length(x) - sum(is.na(x)) == 2L 
}

#' Detect single level variable
#' 
#' Detect a single level variable, for purposes of factor conversion.
#' 
#' @param v The variable of interest.
#' @export

is.singular <- function(v) { 
  x <- unique(v) 
  length(x) - sum(is.na(x)) == 1L 
}