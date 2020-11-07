#' Summarizing Bayesian Multilevel Single Case objects
#'
#' summary method for class "BMSC".
#'
#'
#' @param object An object of class \code{BMSC}, resulting from the \link{BMSC} function.
#'
#' @param ... other arguments are ignored.
#' @method summary BMSC
#' @return a \code{summary.BMSC} object
#' @export
summary.BMSC = function(object, ...) {

    if (class(object)[2] != "BMSC")
        stop("Not a valid BMSC object.")

    se <- function(object) {
        sd(object)/sqrt(length(object))
    }

    if (object[[7]] == "normal") {
        d0 <- dnorm(0, 0, object[[8]])
    } else if (object[[7]] == "cauchy") {
        d0 <- dcauchy(0, 0, object[[8]])
    } else if (object[[7]] == "student") {
        d0 <- LaplacesDemon::dst(0, object[[8]], 3)
    }

    delta = extract(object[[2]], pars = "b_Delta")
    delta_logspl = apply(delta$b_Delta, 2, .suppresslogspline)
    BF10_delta = lapply(delta_logspl, FUN = function(x) {
        d0/.suppressdlogspline(0, x)
    })

    beta = extract(object[[2]], pars = "b_Ctrl")
    beta_logspl = apply(beta$b_Ctrl, 2, .suppresslogspline)
    BF10_beta = lapply(beta_logspl, FUN = function(x) {
        d0/.suppressdlogspline(0, x)
    })

    pts = beta$b_Ctrl + delta$b_Delta
    pts_logspl = apply(pts, 2, .suppresslogspline)
    BF10_pts = lapply(pts_logspl, FUN = function(x) {
        d0/.suppressdlogspline(0, x)
    })

    sum01 = as.data.frame(summary(object[[2]], pars = "b_Ctrl")[[1]])
    sum02 = as.data.frame(summary(object[[2]], pars = "sigmaC")[[1]])

    sum03 = as.data.frame(summary(object[[2]], pars = "b_Delta")[[1]])
    sum04 = as.data.frame(summary(object[[2]], pars = "sigmaP")[[1]])

    sum05 = as.data.frame(cbind(apply(pts, 2, mean), apply(pts, 2, se), apply(pts, 2, sd), apply(pts, 2, quantile, probs = 2.5/100), apply(pts, 2, quantile,
        probs = 25/100), apply(pts, 2, quantile, probs = 50/100), apply(pts, 2, quantile, probs = 75/100), apply(pts, 2, quantile, probs = 97.5/100)))

    colnames(sum05) = c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")

    rownames(sum01) <- colnames(object[[5]]$XF_Ctrl)
    rownames(sum03) <- colnames(object[[5]]$XF_Pts)
    rownames(sum05) <- colnames(object[[5]]$XF_Pts)

    sum01$BF10 <- BF10_beta

    sum03$BF10 <- BF10_delta

    sum05$BF10 <- BF10_pts

    out = list(sum01, sum02, sum03, sum04, object, sum05, object[[7]], object[[8]])

    class(out) = append(class(out), "summary.BMSC")

    return(out)
}

#' Print summaries of Bayesian Multilevel Single Case objects
#'
#'
#' @param x An object of class \code{summary.BMSC}, resulting from the \link{summary.BMSC} function.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @method print summary.BMSC
#' @export
print.summary.BMSC = function(x, ...) {
    cat("\nBayesian Multilevel Single Case model\n\n")

    print(x[[5]][[1]], ...)
    cat("\n")
    print(paste0("Priors for the regression coefficients: ", x[[7]],
                " distribution; Dispersion parameter (scale or sigma): ", x[[8]]))

    cat("\n\n  Fixed Effects for the Control Group\n\n")

    print(x[[1]], ...)
    cat("\n")
    print(x[[2]], ...)

    cat("\n\n  Fixed Effects for the Patient\n\n")

    print(x[[6]], ...)

    cat("\n\n  Fixed Effects for the difference between the Patient and the Control Group\n\n")

    print(x[[3]], ...)
    cat("\n")
    print(x[[4]], ...)
}
