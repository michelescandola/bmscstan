#' @export
summary.BMSC = function(x) {
    
    if (class(x)[2] != "BMSC") 
        stop("Not a valid BMSC object.")
    
    se <- function(x) {
        sd(x)/sqrt(length(x))
    }
    
    if (x[[7]] == "normal") {
        d0 <- dnorm(0, 0, 10)
    } else if (x[[7]] == "cauchy") {
        d0 <- dcauchy(0, 0, sqrt(2)/2)
    } else if (x[[7]] == "student") {
        d0 <- LaplacesDemon::dst(0, 10, 3)
    }
    
    delta = extract(x[[2]], pars = "b_Delta")
    delta_logspl = apply(delta$b_Delta, 2, logspline)
    BF10_delta = lapply(delta_logspl, FUN = function(x) {
        d0/dlogspline(0, x)
    })
    
    beta = extract(x[[2]], pars = "b_Ctrl")
    beta_logspl = apply(beta$b_Ctrl, 2, logspline)
    BF10_beta = lapply(beta_logspl, FUN = function(x) {
        d0/dlogspline(0, x)
    })
    
    pts = beta$b_Ctrl + delta$b_Delta
    pts_logspl = apply(pts, 2, logspline)
    BF10_pts = lapply(pts_logspl, FUN = function(x) {
        d0/dlogspline(0, x)
    })
    
    sum01 = as.data.frame(summary(x[[2]], pars = "b_Ctrl")[[1]])
    sum02 = as.data.frame(summary(x[[2]], pars = "sigmaC")[[1]])
    
    sum03 = as.data.frame(summary(x[[2]], pars = "b_Delta")[[1]])
    sum04 = as.data.frame(summary(x[[2]], pars = "sigmaP")[[1]])
    
    sum05 = as.data.frame(cbind(apply(pts, 2, mean), apply(pts, 2, se), apply(pts, 2, sd), apply(pts, 2, quantile, probs = 2.5/100), apply(pts, 2, quantile, 
        probs = 25/100), apply(pts, 2, quantile, probs = 50/100), apply(pts, 2, quantile, probs = 75/100), apply(pts, 2, quantile, probs = 97.5/100)))
    
    colnames(sum05) = c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
    
    rownames(sum01) <- colnames(x[[5]]$XF_Ctrl)
    rownames(sum03) <- colnames(x[[5]]$XF_Pts)
    rownames(sum05) <- colnames(x[[5]]$XF_Pts)
    
    sum01$BF10 <- BF10_beta
    
    sum03$BF10 <- BF10_delta
    
    sum05$BF10 <- BF10_pts
    
    out = list(sum01, sum02, sum03, sum04, x, sum05, x[[7]])
    
    class(out) = append(class(out), "summary.BMSC")
    
    return(out)
}

#' @export
print.summary.BMSC = function(x, ...) {
    cat("\nBayesian Multilevel Single Case model\n\n")
    
    print(x[[5]][[1]], ...)
    cat("\n")
    print(paste("Prior:", x[[7]]))
    
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