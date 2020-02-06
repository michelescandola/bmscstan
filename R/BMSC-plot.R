#' Posterior predictive checks for BMSC objects
#'
#' \code{pp_check()} plots the posterior predictive check for BMSC objects.
#'
#' @param object a BMSC object
#' @param type
#' \describe{
#'         \item{dens} {density overlay plot}
#'         \item{hist} {histogram plot}
#'         \item{mean} {the distribution of the mean statistic, over the simulated datasets, compared to the mean of the real data}
#'     }
#' @param limited logical. TRUE if the output should be limited within 4 standard deviations, FALSE it should not. default TRUE
#' @return a ggplot2 object
#' @export

pp_check.BMSC = function(object, type = "dens", limited = TRUE) {

    if (class(object)[2] != "BMSC")
        stop("Not a valid BMSC object.")

    y_ct_rep = extract(object[[2]], pars = "y_ct_rep")
    y_pt_rep = extract(object[[2]], pars = "y_pt_rep")

    y_ct = object[[4]][, as.character(object[[1]][2])]
    y_pt = object[[3]][, as.character(object[[1]][2])]

    ans = list()

    if ((requireNamespace("bayesplot", quietly = TRUE)) && (requireNamespace("gridExtra", quietly = TRUE))) {
        if (type == "hist") {
            ct = bayesplot::ppc_hist(y_ct, y_ct_rep[[1]][1:8, ]) + ggtitle("Control Group")
            pt = bayesplot::ppc_hist(y_pt, y_pt_rep[[1]][1:8, ]) + ggtitle("Patient")

            if (limited) {
                ct = ct + xlim(c(-sd(y_ct) * 4, sd(y_ct) * 4))
                pt = pt + xlim(c(-sd(y_pt) * 4, sd(y_pt) * 4))
            }

            ans = gridExtra::grid.arrange(ct, pt)
        } else if (type == "mean") {
            ct = bayesplot::ppc_stat(y_ct, y_ct_rep[[1]]) + ggtitle("Control Group")
            pt = bayesplot::ppc_stat(y_pt, y_pt_rep[[1]]) + ggtitle("Patient")

            if (limited) {
                ct = ct + xlim(c(-sd(y_ct) * 4, sd(y_ct) * 4))
                pt = pt + xlim(c(-sd(y_pt) * 4, sd(y_pt) * 4))
            }

            ans = gridExtra::grid.arrange(ct, pt)
        } else {
            ct = bayesplot::ppc_dens_overlay(y_ct, y_ct_rep[[1]][1:200, ]) + ggtitle("Control Group")
            pt = bayesplot::ppc_dens_overlay(y_pt, y_pt_rep[[1]][1:200, ]) + ggtitle("Patient")

            if (limited) {
                ct = ct + xlim(c(-sd(y_ct) * 4, sd(y_ct) * 4))
                pt = pt + xlim(c(-sd(y_pt) * 4, sd(y_pt) * 4))
            }

            ans = gridExtra::grid.arrange(ct, pt)
        }
    }

    return(ans)
}

#' Plot estimates from a \code{BMSC} object.
#'
#' @usage plot(object, who = 'both', type = 'interval')
#'
#' @param object An object of class \code{BMSC}.
#' @param who parameter to choose the estimates to plot
#'     \describe{
#'         \item{both} {plot in the same graph both controls and the patient}
#'         \item{control} {only the controls}
#'         \item{patient} {only the patient \eqn{(\beta + \delta)}}
#'         \item{delta} {only the difference between the patient and controls}
#'     }
#' @param type a parameter to select the typology of graph
#'     \describe{
#'         \item{interval} {the estimates will be represented by means of pointrange, with median and the boundaries of the credible interval}
#'         \item{area} {a density plot}
#'         \item{hist} {a density histogram}
#'     }
#' @param CI the dimension of the Credible Interval
#'
#' @return a plot, a ggplot2 object, or a bayesplot object
#'
#' @export
plot.BMSC = function(mdl, who = "both", type = "interval", CI = 0.95) {

    if (class(mdl)[2] != "BMSC")
        stop("Not a valid BMSC object.")

    limits = c((1 - CI)/2, 1 - ((1 - CI)/2))

    if (type == "interval") {
        if (who == "both") {
            beta <- extract(mdl[[2]], pars = "b_Ctrl")
            controls <- apply(beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))
            delta <- extract(mdl[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]] + beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namC <- colnames(mdl[[5]]$XF_Ctrl)
            namP <- colnames(mdl[[5]]$XF_Pts)

            dat <- data.frame(low = c(patient[1, ], controls[1, ]), med = c(patient[2, ], controls[2, ]), high = c(patient[3, ], controls[3, ]), group = c(rep("Patient",
                ncol(patient)), rep("Controls", ncol(controls))), i = ordered(c(namP, namC), level = namC[length(namC):1]))

            ans <- ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med, colour = group)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "control") {
            beta <- extract(mdl[[2]], pars = "b_Ctrl")
            controls <- apply(beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namC <- colnames(mdl[[5]]$XF_Ctrl)

            dat <- data.frame(low = controls[1, ], med = controls[2, ], high = controls[3, ], i = ordered(namC, level = namC[length(namC):1]))

            ans <- ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "patient") {
            beta <- extract(mdl[[2]], pars = "b_Ctrl")
            delta <- extract(mdl[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]] + beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namP <- colnames(mdl[[5]]$XF_Pts)

            dat <- data.frame(low = patient[1, ], med = patient[2, ], high = patient[3, ], i = ordered(namP, level = namP[length(namP):1]))

            ans <- ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "delta") {
            delta <- extract(mdl[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namP <- colnames(mdl[[5]]$XF_Pts)

            dat <- data.frame(low = patient[1, ], med = patient[2, ], high = patient[3, ], i = ordered(namP, level = namP[length(namP):1]))

            ans <- ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        }
    } else if (type == "area") {
        if (requireNamespace("reshape2", quietly = TRUE)) {
            if (who == "both") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(mdl[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])
                controls$group = "Controls"

                delta <- extract(mdl[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])
                patient$group = "Patient"

                dat <- rbind(patient, controls)

                ans <- ggplot(dat, aes(x = value, colour = group)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "control") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(mdl[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])

                dat <- controls

                ans <- ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            } else if (who == "patient") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                delta <- extract(mdl[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "delta") {
                delta <- extract(mdl[[2]], pars = "b_Delta")
                patient <- reshape2::melt(delta[[1]])
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            }
        }
    } else if (type == "hist") {
        if (requireNamespace("reshape2", quietly = TRUE)) {
            if (who == "both") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(mdl[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])
                controls$group = "Controls"

                delta <- extract(mdl[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])
                patient$group = "Patient"

                dat <- rbind(patient, controls)

                ans <- ggplot(dat, aes(x = value, fill = group)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "control") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(mdl[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])

                dat <- controls

                ans <- ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            } else if (who == "patient") {
                beta <- extract(mdl[[2]], pars = "b_Ctrl")
                delta <- extract(mdl[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "delta") {
                delta <- extract(mdl[[2]], pars = "b_Delta")
                patient <- reshape2::melt(delta[[1]])
                namP <- colnames(mdl[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            }
        }
    }

    return(ans)
}
