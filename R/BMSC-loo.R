#' @title loo and waic.
#'
#' @description bmscstan wrapper for computing approximate leave-one-out cross-validation
#' (loo) and Watanabe-Akaike Information Criterion or
#' Widely Applicable Information Criterion (WAIC)
#' using PSIS-LOO
#' for the single case and the control group
#'
#'
#' @param x An object of class \code{BMSC}, resulting from the
#'        \link{BMSC} function.
#' @param cores The number of cores for the
#'        `loo::relative_eff` function
#'
#' @param ... for `BMSC_loo` and `BMSC_waic`
#'        further arguments passed to the `loo::extract_log_lik` function.
#'        for `print` and `plot` methods further arguments to be passed to the
#'        `print` or `plot` functions
#'
#' @return for `BMSC_loo` a list with the log likelihood of the single case and the
#'         control group,
#'         the MCMC effective sample size divided by the total sample size,
#'         and the leave-one-out cross-validation.
#'         For `BMSC_waic` a list with the log likelihood of the single
#'         case and the control group, and the waic scores.
#' @export
BMSC_loo = function( x, cores = 1, ...){

  if(!inherits(x, "BMSC"))
    stop("Not a valid BMSC object.")

  log_lik_sc <- loo::extract_log_lik(x[[2]],
                                parameter_name = "log_lik_pt",
                                merge_chains = FALSE,
                                ...)

  r_eff_sc <- loo::relative_eff(exp(log_lik_sc), cores = cores)

  loo_sc <- loo::loo(log_lik_sc, r_eff = r_eff_sc, cores = cores)

  log_lik_ct <- loo::extract_log_lik(x[[2]],
                                parameter_name = "log_lik_ct",
                                merge_chains = FALSE,
                                ...)

  r_eff_ct <- loo::relative_eff(exp(log_lik_ct), cores = cores)

  loo_ct <- loo::loo(log_lik_ct, r_eff = r_eff_ct, cores = cores)

  out <- list(
    log_lik_sc,
    r_eff_sc,
    loo_sc,
    log_lik_ct,
    r_eff_ct,
    loo_ct
  )

  class(out) <- append( class(out), "loo_BMSC")

  return( out )
}

#' @rdname BMSC_loo
#' @export
BMSC_waic = function( x, ...){

  if(!inherits(x, "BMSC"))
    stop("Not a valid BMSC object.")

  log_lik_sc <- loo::extract_log_lik(x[[2]],
                                     parameter_name = "log_lik_pt",
                                     merge_chains = FALSE,
                                     ...)

  waic_sc <- loo::waic( log_lik_sc )

  log_lik_ct <- loo::extract_log_lik(x[[2]],
                                     parameter_name = "log_lik_ct",
                                     merge_chains = FALSE,
                                     ...)

  waic_ct <- loo::waic( log_lik_ct )

  out <- list(
    log_lik_sc,
    waic_sc,
    log_lik_ct,
    waic_ct
  )

  class(out) <- append( class(out), "waic_BMSC")

  return( out )
}

#' @rdname BMSC_loo
#' @method plot loo_BMSC
#' @exportS3Method plot loo_BMSC
plot.loo_BMSC = function( x, ... ){
  op <- graphics::par( no.readonly = TRUE )

  graphics::par( mfrow = c( 1 , 2 ) )
  plot( x[[3]], main = "PSIS diagnostic plot\nSingle Case", ...)

  plot( x[[3]], main = "PSIS diagnostic plot\nControl group", ...)

  graphics::par( op )
}

#' @rdname BMSC_loo
#' @method print waic_BMSC
#' @exportS3Method print waic_BMSC
print.waic_BMSC = function( x, ... ){
  cat("\nWidely applicable information criterion (WAIC)
      for the single case\n\n")

  print( x[[2]], ... )

  cat("\nWidely applicable information criterion (WAIC)
      for the control group\n\n")

  print( x[[4]], ... )
}

#' @rdname BMSC_loo
#' @method print loo_BMSC
#' @exportS3Method print loo_BMSC
print.loo_BMSC = function( x, ... ){
  cat("\nLeave-One-Out Cross-Validation using PSIS-LOO for the single case\n\n")

  print( x[[3]], ... )

  cat("\nLeave-One-Out Cross-Validation using PSIS-LOO for the control group\n\n")

  print( x[[6]], ... )
}

#' @title bmscstan wrapper for model comparison.
#'
#' @param x A list of \code{loo_BMSC} or  \code{waic_BMSC} objects.
#' @param simplify For the print method only, should only the essential columns
#'       of the summary matrix be printed? The entire matrix is always returned,
#'       but by default only the most important columns are printed.
#' @param ... further arguments passed to the function.
#'
#' @return a list with the log likelihood of the single case and the
#'         control group,
#'         the MCMC effective sample size divided by the total sample size,
#'         and the leave-one-out cross-validation.
#'
#' @export
BMSC_loo_compare = function( x, ... ){

  x_ct <- x_sc <- list()

  class_output <- NULL

  for(i in 1:length(x)){

    if ( !inherits( x[[i]], "loo_BMSC") && !inherits( x[[i]], "waic_BMSC") )
      stop(paste0(i, "-th is not a valid loo_BMSC or waic_BMSC object."))

    x_sc[[ i ]] <- x[[i]][[3]]
    x_ct[[ i ]] <- x[[i]][[3]]
  }

  if ( inherits( x[[i]], "loo_BMSC") ){
    class_output <- "loo_compare_BMSC"
  } else {
    class_output <- "waic_compare_BMSC"
  }

  comp_sc <- loo::loo_compare( x_sc, ... )
  comp_ct <- loo::loo_compare( x_ct, ... )

  tmp <- list(
    comp_sc,
    comp_ct
  )

  class( tmp ) <- append( class( tmp ), class_output )

  return( tmp )
}

#' @rdname BMSC_loo_compare
#' @method print loo_compare_BMSC
#' @exportS3Method print loo_compare_BMSC
print.loo_compare_BMSC = function( x , simplify = TRUE, ... ){
  cat("\nLeave-One-Out Cross-Validation model comparison for the single case\n\n")

  print( x[[1]] , simplify = simplify, ... )

  cat("\nLeave-One-Out Cross-Validation model comparison for the control group\n\n")

  print( x[[2]] , simplify = simplify, ... )
}

#' @rdname BMSC_loo_compare
#' @method print waic_compare_BMSC
#' @exportS3Method print waic_compare_BMSC
print.waic_compare_BMSC = function( x , simplify = TRUE, ... ){
  cat("\nWidely applicable information criterion (WAIC)
      for the single case\n\n")

  print( x[[1]] , simplify = simplify, ... )

  cat("\nWidely applicable information criterion (WAIC)
      for the control group\n\n")

  print( x[[2]] , simplify = simplify, ...)
}

#' @title bmscstan wrapper for diagnostics for Pareto smoothed importance
#' sampling (PSIS)
#'
#' @param x An object of class \code{loo_BMSC}
#' @param threshold for the `pareto_k_ids` method is the minimum $k$ value to flag.
#'        for the `mcse_loo` method all the $k$ values greater than the
#'        `threshold` will be returned as NA.
#' @param ... further arguments passed to the `print` function.
#'
#' @return
#' \itemize{
#'   \item{pareto_k_table}{ returns an object of class "pareto_k_table_BMSC",
#'      which is a matrix with columns "Count", "Proportion", and "Min. n_eff"
#'      }
#'   \item{pareto_k_ids}{ returns an integer vector indicating which observations
#'      have Pareto k estimates above threshold
#'      }
#'   \item{mcse_loo}{ returns the Monte Carlo standard error (MCSE) estimate for
#'      PSIS-LOO.
#'      MCSE will be NA if any Pareto kk values are above threshold.
#'      }
#'   \item{pareto_k_values}{ returns a vector of the estimated Pareto k parameters.
#'      These represent the reliability of sampling.
#'      }
#'   \item{pareto_k_influence_values}{ returns a vector of the estimated Pareto
#'      k parameters.
#'      These represent influence of the observations on the model
#'      posterior distribution.
#'      }
#'   \item{psis_k_influence_table}{ returns a vector of the estimated PSIS
#'      effective sample sizes.
#'      }
#' }
#'
#' @export
BMSC_pareto_k_table = function( x ){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::pareto_k_table( x[[3]] ),
    "Control group" = loo::pareto_k_table( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_table_BMSC")

  return( out )

}

#' @rdname BMSC_pareto_k_table
#' @export
BMSC_pareto_k_ids = function( x , threshold = 0.5){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::pareto_k_ids( x[[3]] , threshold = threshold ),
    "Control group" = loo::pareto_k_ids( x[[6]] , threshold = threshold )
  )

  class( out ) <- append( class( out ) , "pareto_k_ids_BMSC")

  return( out )
}

#' @rdname BMSC_pareto_k_table
#' @export
BMSC_mcse_loo = function( x, threshold = 0.7 ){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::mcse_loo( x[[3]], threshold = threshold ),
    "Control group" = loo::mcse_loo( x[[6]], threshold = threshold )
  )

  class( out ) <- append( class( out ) , "mcse_loo_BMSC")

  return( out )
}

#' @rdname BMSC_pareto_k_table
#' @export
BMSC_pareto_k_values = function( x ){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::pareto_k_values( x[[3]] ),
    "Control group" = loo::pareto_k_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_values_BMSC")

  return( out )
}

#' @rdname BMSC_pareto_k_table
#' @export
BMSC_pareto_k_influence_values = function( x ){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::pareto_k_influence_values( x[[3]] ),
    "Control group" = loo::pareto_k_influence_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_influence_values_BMSC")

  return( out )
}

#' @rdname BMSC_pareto_k_table
#' @export
BMSC_psis_n_eff_values = function( x ){

  if(!inherits(x, "loo_BMSC"))
    stop("Not a valid loo_BMSC object.")

  out <- list(
    "Single Case" = loo::psis_n_eff_values( x[[3]] ),
    "Control group" = loo::psis_n_eff_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "psis_n_eff_values_BMSC")

  return( out )
}

#' @rdname BMSC_pareto_k_table
#' @method print pareto_k_table_BMSC
#' @export
#' @exportS3Method print pareto_k_table_BMSC
print.pareto_k_table_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}

#' @rdname BMSC_pareto_k_table
#' @method print pareto_k_ids_BMSC
#' @export
#' @exportS3Method print pareto_k_ids_BMSC
print.pareto_k_ids_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}

#' @rdname BMSC_pareto_k_table
#' @method print pareto_k_values_BMSC
#' @export
#' @exportS3Method print pareto_k_values_BMSC
print.pareto_k_values_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}

#' @rdname BMSC_pareto_k_table
#' @method print pareto_k_influence_values_BMSC
#' @export
#' @exportS3Method print pareto_k_influence_values_BMSC
print.pareto_k_influence_values_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}

#' @rdname BMSC_pareto_k_table
#' @method print psis_n_eff_values_BMSC
#' @export
#' @exportS3Method print psis_n_eff_values_BMSC
print.psis_n_eff_values_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}

#' @rdname BMSC_pareto_k_table
#' @method print mcse_loo_BMSC
#' @export
#' @exportS3Method print mcse_loo_BMSC
print.mcse_loo_BMSC = function( x, ... ){

  cat("\nSingle case\n\n")

  print( x[[1]], ... )

  cat("\nControl group\n\n")

  print( x[[2]], ... )
}
