#' bmscstan wrapper for computing approximate leave-one-out cross-validation
#' using PSIS-LOO
#' for the single case and the control group
#'
#'
#' @param mdl An object of class \code{BMSC}, resulting from the
#'        \link{BMSC} function.
#' @param cores The number of cores fo the \link{loo::relative_eff} function
#'
#' @param ... further arguments passed to the `loo::extract_log_lik` function.
#'
#' @method loo BMSC
#' @return a list with the log likelihood of the single case and the
#'         control group,
#'         the MCMC effective sample size divided by the total sample size,
#'         and the leave-one-out cross-validation.
#' @export
loo.BMSC = function( mdl, cores = 1, ...){

  log_lik_sc <- loo::extract_log_lik(mdl[[2]],
                                parameter_name = "log_lik_pt",
                                merge_chains = FALSE,
                                ...)

  r_eff_sc <- loo::relative_eff(exp(log_lik_sc), cores = cores)

  loo_sc <- loo::loo(log_lik_sc, r_eff = r_eff_sc, cores = cores)

  log_lik_ct <- loo::extract_log_lik(mdl[[2]],
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

#' @export
loo_compare.loo_BMSC = function( mdl1, mdl2 ){

  if ( class( mdl1 )[2] != "loo_BMSC" )
    stop("mdl1 is not a valid loo_BMSC object.")
  if ( class( mdl2 )[2] != "loo_BMSC" )
    stop("mdl2 is not a valid loo_BMSC object.")

  comp_sc <- loo::loo_compare( mdl1[[3]], mdl2[[3]] )
  comp_ct <- loo::loo_compare( mdl1[[6]], mdl2[[6]] )

  tmp <- list(
    comp_sc,
    comp_ct
  )

  class( tmp ) <- append( class( tmp ), "loo_compare_BMSC" )

  return( tmp )
}

#' @export
print.loo_BMSC = function( x ){
  cat("\nLeave-One-Out Cross-Validation using PSIS-LOO for the single case\n\n")

  print( x[[3]] )

  cat("\nLeave-One-Out Cross-Validation using PSIS-LOO for the control group\n\n")

  print( x[[6]] )
}

#' @export
print.loo_compare_BMSC = function( x , simplify = TRUE ){
  cat("\nLeave-One-Out Cross-Validation model comparison for the single case\n\n")

  print( x[[1]] , simplify = simplify )

  cat("\nLeave-One-Out Cross-Validation model comparison for the control group\n\n")

  print( x[[2]] , simplify = simplify )
}

#' @export
pareto_k_table.loo_BMSC = function( x ){

  out <- list(
    "Single Case" = loo::pareto_k_table( x[[3]] ),
    "Control group" = loo::pareto_k_table( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_table_BMSC")

  return( out )

}

#' @export
print.pareto_k_table_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
pareto_k_ids.loo_BMSC = function( x , threshold = 0.5){

  out <- list(
    "Single Case" = loo::pareto_k_ids( x[[3]] , threshold = threshold ),
    "Control group" = loo::pareto_k_ids( x[[6]] , threshold = threshold )
  )

  class( out ) <- append( class( out ) , "pareto_k_ids_BMSC")

  return( out )
}

#' @export
print.pareto_k_ids_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
pareto_k_values.loo_BMSC = function( x ){

  out <- list(
    "Single Case" = loo::pareto_k_values( x[[3]] ),
    "Control group" = loo::pareto_k_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_values_BMSC")

  return( out )
}

#' @export
print.pareto_k_values_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
pareto_k_influence_values.loo_BMSC = function( x ){

  out <- list(
    "Single Case" = loo::pareto_k_influence_values( x[[3]] ),
    "Control group" = loo::pareto_k_influence_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "pareto_k_influence_values_BMSC")

  return( out )
}

#' @export
print.pareto_k_influence_values_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
psis_n_eff_values.loo_BMSC = function( x ){

  out <- list(
    "Single Case" = loo::psis_n_eff_values( x[[3]] ),
    "Control group" = loo::psis_n_eff_values( x[[6]] )
  )

  class( out ) <- append( class( out ) , "psis_n_eff_values_BMSC")

  return( out )
}

#' @export
print.psis_n_eff_values_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
mcse_loo.loo_BMSC = function( x, threshold = 0.7 ){

  out <- list(
    "Single Case" = loo::mcse_loo( x[[3]], threshold = threshold ),
    "Control group" = loo::mcse_loo( x[[6]], threshold = threshold )
  )

  class( out ) <- append( class( out ) , "mcse_loo_BMSC")

  return( out )
}

#' @export
print.mcse_loo_BMSC = function( x ){

  cat("\nSingle case\n\n")

  print( x[[1]] )

  cat("\nControl group\n\n")

  print( x[[2]] )
}

#' @export
plot.loo_BMSC = function( x ){
  op <- par( no.readonly = TRUE )

  par( mfrow = c( 1 , 2 ) )
  plot( x[[3]], main = "PSIS diagnostic plot\nSingle Case")

  plot( x[[3]], main = "PSIS diagnostic plot\nControl group")

  par( op )
}

#' bmscstan wrapper to compute the Widely applicable information criterion
#' (WAIC)
#' for the single case and the control group
#'
#'
#' @param mdl An object of class \code{BMSC}, resulting from the
#'            \link{BMSC} function.
#'
#' @param ... further arguments passed to the `loo::extract_log_lik` function.
#'
#' @method waic BMSC
#' @return a list with the log likelihood of the single case
#'         and the control group,
#'         and the waic scores.
#' @export
waic.BMSC = function( mdl, ...){

  log_lik_sc <- loo::extract_log_lik(mdl[[2]],
                                     parameter_name = "log_lik_pt",
                                     merge_chains = FALSE,
                                     ...)

  waic_sc <- loo::waic( log_lik_sc )

  log_lik_ct <- loo::extract_log_lik(mdl[[2]],
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

#' @export
print.waic_BMSC = function( x ){
  cat("\nWidely applicable information criterion (WAIC)
      for the single case\n\n")

  print( x[[2]] )

  cat("\nWidely applicable information criterion (WAIC)
      for the control group\n\n")

  print( x[[4]] )
}

#' @export
loo_compare.waic_BMSC = function( mdl1, mdl2 ){

  if ( class( mdl1 )[2] != "waic_BMSC" )
    stop("mdl1 is not a valid waic_BMSC object.")
  if ( class( mdl2 )[2] != "waic_BMSC" )
    stop("mdl2 is not a valid waic_BMSC object.")

  comp_sc <- loo::loo_compare( mdl1[[3]], mdl2[[3]] )
  comp_ct <- loo::loo_compare( mdl1[[6]], mdl2[[6]] )

  tmp <- list(
    comp_sc,
    comp_ct
  )

  class( tmp ) <- append( class( tmp ), "waic_compare_BMSC" )

  return( tmp )
}

#' @export
print.waic_compare_BMSC = function( x , simplify = TRUE ){
  cat("\nWidely applicable information criterion (WAIC)
      for the single case\n\n")

  print( x[[1]] , simplify = simplify )

  cat("\nWidely applicable information criterion (WAIC)
      for the control group\n\n")

  print( x[[2]] , simplify = simplify )
}
