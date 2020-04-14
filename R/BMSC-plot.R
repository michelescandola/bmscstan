#' Posterior predictive checks for BMSC objects
#'
#' \code{pp_check()} plots the posterior predictive check for BMSC objects.
#'
#' @param object a BMSC object
#' @param type
#' \describe{
#'         \item{dens}{density overlay plot}
#'         \item{hist}{histogram plot}
#'         \item{mean}{the distribution of the mean statistic, over the simulated datasets, compared to the mean of the real data}
#'     }
#' @param limited logical. TRUE if the output should be limited within 4 standard deviations, FALSE it should not. default FALSE
#' @param ... other arguments for \code{bayesplot} methods.
#'
#' @examples
#'  \dontrun{
#' # simulation of healthy controls data
#'
#' Sigma.ctrl <- matrix(cbind(1, .7,  .7, 1) ,nrow=2)
#'
#' U <- t(chol(Sigma.ctrl))
#'
#' numobs <- 100
#'
#' set.seed(123)
#'
#' random.normal <- matrix( rnorm( n = ncol(U) * numobs, mean = 3, sd = 1),
#'                          nrow = ncol(U), ncol = numobs)
#'
#' X = U %*% random.normal
#'
#' dat.ctrl <- as.data.frame(t(X))
#'
#' names(dat.ctrl) <- c("y","x")
#'
#' cor(dat.ctrl)
#'
#' # simulation of patient data
#'
#' Sigma.pt <- matrix(cbind(1, 0,  0, 1) ,nrow=2)
#'
#' U <- t(chol(Sigma.pt))
#'
#' numobs <- 20
#'
#' set.seed(0)
#'
#' random.normal <- matrix( rnorm( n = ncol(U) * numobs, mean = 3, sd = 1),
#'                  nrow = ncol(U), ncol = numobs)
#'
#' X = U %*% random.normal
#'
#' dat.pt <- as.data.frame(t(X))
#'
#' names(dat.pt) <- c("y","x")
#'
#' cor(dat.pt)
#'
#' # fit the single case model
#'
#' mdl.reg <- BMSC(y ~ x, data_ctrl = dat.ctrl, data_pt = dat.pt, seed = 10)
#'
#' # summarize the data
#'
#' summary(mdl.reg)
#'
#' # plot the posterior predictive checks
#'
#' pp_check(mdl.reg, limited = FALSE)
#'
#' pp_check(mdl.reg, limited = TRUE)
#'
#' pp_check(mdl.reg, type = "mean", limited = FALSE)
#'
#' pp_check(mdl.reg, type = "hist", limited = FALSE)
#' }
#'
#' @return a ggplot2 object
#' @export
pp_check.BMSC = function(object, type = "dens", limited = FALSE, ...) {

    if (class(object)[2] != "BMSC")
        stop("Not a valid BMSC object.")

    y_ct_rep = extract(object[[2]], pars = "y_ct_rep")
    y_pt_rep = extract(object[[2]], pars = "y_pt_rep")

    y_ct = object[[4]][, as.character(object[[1]][2])]
    y_pt = object[[3]][, as.character(object[[1]][2])]

    ans = list()

    if ((requireNamespace("bayesplot", quietly = TRUE)) && (requireNamespace("gridExtra", quietly = TRUE))) {
        if (type == "hist") {
            ct = bayesplot::ppc_hist(y_ct, y_ct_rep[[1]][1:8, ], ...) + ggtitle("Control Group")
            pt = bayesplot::ppc_hist(y_pt, y_pt_rep[[1]][1:8, ], ...) + ggtitle("Patient")

            if (limited) {
                ct = ct + xlim(c(-sd(y_ct) * 4, sd(y_ct) * 4))
                pt = pt + xlim(c(-sd(y_pt) * 4, sd(y_pt) * 4))
            }

            ans = gridExtra::grid.arrange(ct, pt)
        } else if (type == "mean") {
            ct = bayesplot::ppc_stat(y_ct, y_ct_rep[[1]], ...) + ggtitle("Control Group")
            pt = bayesplot::ppc_stat(y_pt, y_pt_rep[[1]], ...) + ggtitle("Patient")

            if (limited) {
                ct = ct + xlim(c(-sd(y_ct) * 4, sd(y_ct) * 4))
                pt = pt + xlim(c(-sd(y_pt) * 4, sd(y_pt) * 4))
            }

            ans = gridExtra::grid.arrange(ct, pt)
        } else {
            ct = bayesplot::ppc_dens_overlay(y_ct, y_ct_rep[[1]][1:200, ], ...) + ggtitle("Control Group")
            pt = bayesplot::ppc_dens_overlay(y_pt, y_pt_rep[[1]][1:200, ], ...) + ggtitle("Patient")

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
#' @param x An object of class \link{BMSC}.
#' @param who parameter to choose the estimates to plot
#' \describe{
#'         \item{both}{plot in the same graph both controls and the patient}
#'         \item{control}{only the controls}
#'         \item{patient}{only the patient \eqn{(\beta + \delta)}}
#'         \item{delta}{only the difference between the patient and controls}
#' }
#' @param type a parameter to select the typology of graph
#' \describe{
#'         \item{interval}{the estimates will be represented by means of pointrange, with median and the boundaries of the credible interval}
#'         \item{area}{a density plot}
#'         \item{hist}{a density histogram}
#' }
#' @param CI the dimension of the Credible Interval (or Equally Tailed Interval). Default 0.95.
#' @param ... other arguments are ignored.
#'
#' @examples
#'  \dontrun{
#'
#'  data(BSE)
#'
#'  # Normal robust regression of data coming from a body representation paradigm
#'  # with a control sample of 12 participants and one patient with
#'  # unilateral brachial plexus lesion
#'  mdl <- BMSC(formula = RT ~ Body.District * Congruency * Side +
#'              (Body.District + Congruency + Side | ID),
#'              data_ctrl = data.ctrl,
#'              data_pt = data.pt,
#'              cores = 4)
#'
#'  # generate a summary of the results
#'  summary(mdl)
#'
#'  # posterior predictive p-value checking
#'  pp_check(mdl, limited = FALSE)
#'
#'  # plot of the results
#'  plot(mdl)
#'
#' # simulation of healthy controls data
#'
#' Sigma.ctrl <- matrix(cbind(1, .7,  .7, 1) ,nrow=2)
#'
#' U <- t(chol(Sigma.ctrl))
#'
#' numobs <- 100
#'
#' set.seed(123)
#'
#' random.normal <- matrix( rnorm( n = ncol(U) * numobs, mean = 3, sd = 1),
#'                          nrow = ncol(U), ncol = numobs)
#'
#' X = U %*% random.normal
#'
#' dat.ctrl <- as.data.frame(t(X))
#'
#' names(dat.ctrl) <- c("y","x")
#'
#' cor(dat.ctrl)
#'
#' # simulation of patient data
#'
#' Sigma.pt <- matrix(cbind(1, 0,  0, 1) ,nrow=2)
#'
#' U <- t(chol(Sigma.pt))
#'
#' numobs <- 20
#'
#' set.seed(0)
#'
#' random.normal <- matrix( rnorm( n = ncol(U) * numobs, mean = 3, sd = 1),
#'                  nrow = ncol(U), ncol = numobs)
#'
#' X = U %*% random.normal
#'
#' dat.pt <- as.data.frame(t(X))
#'
#' names(dat.pt) <- c("y","x")
#'
#' cor(dat.pt)
#'
#' # fit the single case model
#'
#' mdl.reg <- BMSC(y ~ x, data_ctrl = dat.ctrl, data_pt = dat.pt, seed = 10)
#'
#' # summarize the data
#'
#' summary(mdl.reg)
#'
#' # plot the results of both patient and control group
#'
#' plot(mdl.reg)
#'
#' # plot the results of the patient
#'
#' plot(mdl.reg, who = "patient")
#'
#' # plot the results of the difference between the control group and the patient
#'
#' plot(mdl.reg, who = "delta")
#'
#' # density plots
#'
#' plot(mdl.reg, type = "area")
#'
#' # histograms
#'
#' plot(mdl.reg, type = "hist")
#' }
#'
#' @return a plot, a ggplot2 object, or a bayesplot object
#'
#' @method plot BMSC
#' @export
plot.BMSC = function(x, who = "both", type = "interval", CI = 0.95, ...) {

    if (class(x)[2] != "BMSC")
        stop("Not a valid BMSC object.")

    limits <- c((1 - CI)/2, 1 - ((1 - CI)/2))

    low <- i <- high <- med <- group <- value <- NULL

    if (type == "interval") {
        if (who == "both") {
            beta <- extract(x[[2]], pars = "b_Ctrl")
            controls <- apply(beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))
            delta <- extract(x[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]] + beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namC <- colnames(x[[5]]$XF_Ctrl)
            namP <- colnames(x[[5]]$XF_Pts)

            dat <- data.frame(low = c(patient[1, ], controls[1, ]), med = c(patient[2, ], controls[2, ]), high = c(patient[3, ], controls[3, ]), group = c(rep("Patient",
                ncol(patient)), rep("Controls", ncol(controls))), i = ordered(c(namP, namC), level = namC[length(namC):1]))

            ans <- ggplot2::ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med, colour = group)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "control") {
            beta <- extract(x[[2]], pars = "b_Ctrl")
            controls <- apply(beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namC <- colnames(x[[5]]$XF_Ctrl)

            dat <- data.frame(low = controls[1, ], med = controls[2, ], high = controls[3, ], i = ordered(namC, level = namC[length(namC):1]))

            ans <- ggplot2::ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "patient") {
            beta <- extract(x[[2]], pars = "b_Ctrl")
            delta <- extract(x[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]] + beta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namP <- colnames(x[[5]]$XF_Pts)

            dat <- data.frame(low = patient[1, ], med = patient[2, ], high = patient[3, ], i = ordered(namP, level = namP[length(namP):1]))

            ans <- ggplot2::ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        } else if (who == "delta") {
            delta <- extract(x[[2]], pars = "b_Delta")
            patient <- apply(delta[[1]], 2, quantile, probs = c(limits[1], 0.5, limits[2]))

            namP <- colnames(x[[5]]$XF_Pts)

            dat <- data.frame(low = patient[1, ], med = patient[2, ], high = patient[3, ], i = ordered(namP, level = namP[length(namP):1]))

            ans <- ggplot2::ggplot(dat, aes(x = i, ymin = low, ymax = high, y = med)) + geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
        }
    } else if (type == "area") {
        if (requireNamespace("reshape2", quietly = TRUE)) {
            if (who == "both") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(x[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])
                controls$group = "Controls"

                delta <- extract(x[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])
                patient$group = "Patient"

                dat <- rbind(patient, controls)

                ans <- ggplot2::ggplot(dat, aes(x = value, colour = group)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "control") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(x[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])

                dat <- controls

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            } else if (who == "patient") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                delta <- extract(x[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "delta") {
                delta <- extract(x[[2]], pars = "b_Delta")
                patient <- reshape2::melt(delta[[1]])
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_density() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            }
        }
    } else if (type == "hist") {
        if (requireNamespace("reshape2", quietly = TRUE)) {
            if (who == "both") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(x[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])
                controls$group = "Controls"

                delta <- extract(x[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])
                patient$group = "Patient"

                dat <- rbind(patient, controls)

                ans <- ggplot2::ggplot(dat, aes(x = value, fill = group)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "control") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                controls <- reshape2::melt(beta[[1]])
                namC <- colnames(x[[5]]$XF_Ctrl)
                controls$Var2 <- factor(controls$Var2)
                levels(controls$Var2) <- namC
                controls$Var2 <- ordered(controls$Var2, level = namC[length(namC):1])

                dat <- controls

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            } else if (who == "patient") {
                beta <- extract(x[[2]], pars = "b_Ctrl")
                delta <- extract(x[[2]], pars = "b_Delta")
                tmp <- delta[[1]] + beta[[1]]
                patient <- reshape2::melt(tmp)
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")

            } else if (who == "delta") {
                delta <- extract(x[[2]], pars = "b_Delta")
                patient <- reshape2::melt(delta[[1]])
                namP <- colnames(x[[5]]$XF_Pts)
                patient$Var2 <- factor(patient$Var2)
                levels(patient$Var2) <- namP
                patient$Var2 <- ordered(patient$Var2, level = namP[length(namP):1])

                dat <- patient

                ans <- ggplot2::ggplot(dat, aes(x = value)) + geom_histogram() + facet_grid(Var2 ~ .) + xlab("coefficients") + ylab("")
            }
        }
    }

    return(ans)
}

#' Plot estimates from a \code{pairwise.BMSC} object.
#'
#' @param x An object of class \link{pairwise.BMSC}.
#' @param type a parameter to select the typology of graph
#' \describe{
#'         \item{interval}{the estimates will be represented by means of pointrange, with median and the boundaries of the credible interval}
#'         \item{area}{a density plot}
#'         \item{hist}{a density histogram}
#' }
#' @param CI the dimension of the Credible Interval (or Equally Tailed Interval). Default 0.95.
#' @param ... other arguments are ignored.
#'
#' @examples
#'  \dontrun{
#'
#'  data(BSE)
#'
#'  # Normal regression of data coming from a body representation paradigm
#'  # with a control sample of 12 participants and one patient with
#'  # unilateral brachial plexus lesion
#'  mdl <- BMSC(formula = RT ~ Body.District * Congruency * Side +
#'              (Body.District + Congruency + Side | ID),
#'              data_ctrl = data.ctrl,
#'              data_pt = data.pt,
#'              cores = 4)
#'
#'  # generate a summary of the results
#'  summary(mdl)
#'
#'  # posterior predictive p-value checking
#'  pp_check(mdl, limited = FALSE)
#'
#'  # plot of the results
#'  plot(mdl)
#'
#'  # compute pairwise contrasts
#'  ph <- pairwise(mdl , contrast = "Body.District1:Side1")
#'
#'  ph
#'
#'  # plot pairwise comparisons
#'
#'  plot(ph)
#'
#'  plot(ph , type = "area")
#'
#'  # customization of pairiwse comparisons plot
#'
#'  plot(ph)[[1]]+theme_bw(base_size = 18)
#'
#'  plot(ph , type = "area")[[1]]+theme_bw(base_size = 18)+
#'    theme(strip.text.y = element_text( angle = 0))
#' }
#'
#' @return a list of two ggplot2 objects
#'
#' @method plot pairwise.BMSC
#' @export
plot.pairwise.BMSC = function(x, type = "interval", CI = 0.95, ...) {
  if (class(x)[2] != "pairwise.BMSC")
    stop("Not a valid pairwise.BMSC object.")

  limits <- c((1 - CI)/2, 1 - ((1 - CI)/2))

  low <- i <- high <- med <- group <- value <- y <- NULL

  # compute the posterior distributions of the contrasts

  deltas.dist     <- list()
  deltas.names    <- list()

  for(i in 1:(length(x[[5]])-1)){
    for(j in (i+1):length(x[[5]])){
      deltas.dist[[paste(i,j)]] <- x[[5]][[i]]$y - x[[5]][[j]]$y

      deltas.names[[paste(i,j)]] <- paste(x[[5]][[i]]$name[1],
                                       x[[5]][[j]]$name[1],
                                       sep = " - ")
    }
  }

  deltas.names <- do.call("c" , deltas.names)

  if (type == "interval") {
    betas  <- do.call("rbind" , lapply(x[[5]], function(x){ quantile(x$y , probs = c(limits[1], 0.5, limits[2]))}))
    deltas <- do.call("rbind" , lapply(deltas.dist, function(x){ quantile(x , probs = c(limits[1], 0.5, limits[2]))}))

    namC <- sapply(x[[5]], function(x){ as.character( x$name[1] ) } )

    dat1 <- data.frame(low = betas[ , 1 ], med = betas[ , 2 ], high = betas[ , 3 ],
                       i = ordered(namC, level = namC[length(namC):1]))
    dat2 <- data.frame(low = deltas[ , 1 ], med = deltas[ , 2 ], high = deltas[ , 3 ],
                       i = ordered(deltas.names, level = deltas.names[length(deltas.names):1]))

    ans <- list(
      ggplot2::ggplot(dat1, aes(x = i, ymin = low, ymax = high, y = med)) +
        geom_pointrange() + coord_flip() + xlab("coefficients") + ylab(""),
      ggplot2::ggplot(dat2, aes(x = i, ymin = low, ymax = high, y = med)) +
        geom_pointrange() + coord_flip() + xlab("coefficients") + ylab("")
    )
  } else if (type == "area") {
        betas  <- do.call("rbind", x[[5]])
        deltas <- data.frame(
          y = do.call("c", deltas.dist),
          contrast = rep(deltas.names , each = length(deltas.dist[[1]]))
          )

        betas$name <- ordered(betas$name, level = levels(betas$name)[length(levels(betas$name)):1])
        deltas$name <- ordered(deltas$contrast, level = levels(deltas$contrast)[length(levels(deltas$contrast)):1])

        dat1 <- betas
        dat2 <- deltas

        ans <- list(
          ggplot2::ggplot(dat1, aes(x = y)) + geom_density() + facet_grid(name ~ .) + xlab("coefficients") + ylab(""),
          ggplot2::ggplot(dat2, aes(x = y)) + geom_density() + facet_grid(contrast ~ .) + xlab("contrast") + ylab("")
        )
  } else if (type == "hist") {
    betas  <- do.call("rbind", x[[5]])
    deltas <- data.frame(
      y = do.call("c", deltas.dist),
      contrast = rep(deltas.names , each = length(deltas.dist[[1]]))
    )

    betas$name <- ordered(betas$name, level = levels(betas$name)[length(levels(betas$name)):1])
    deltas$name <- ordered(deltas$contrast, level = levels(deltas$contrast)[length(levels(deltas$contrast)):1])

    dat1 <- betas
    dat2 <- deltas

    ans <- list(
      ggplot2::ggplot(dat1, aes(x = y)) + geom_histogram() + facet_grid(name ~ .) + xlab("coefficients") + ylab(""),
      ggplot2::ggplot(dat2, aes(x = y)) + geom_histogram() + facet_grid(contrast ~ .) + xlab("contrast") + ylab("")
    )
  }

  return(ans)
}
