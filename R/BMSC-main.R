#' Bayesian Multilevel Single Case models using 'Stan'
#'
#' The \strong{bmscstan} package provides an interface to fit Bayesian
#' Multilevel Single Case models.
#' These models compare the performance of a Single Case against a control
#' group, combining
#' the flexibility of multilevel models and the potentiality of Bayesian
#' Statistics.
#'
#' The package is now limited to gaussian data only, but we will further expand
#' it to cover
#' binomial and ordinal (Likert scales) data.
#'
#' By means of \strong{bmscstan} the effects of the control group and the
#' effects of the deviance between the
#' Single Case and the group will be estimated.
#'
#' The model to estimate the controls parameters is:
#'
#'\ifelse{html}{\out{<center><i>y~N(&beta; X + b Z, &sigma;<sup>2</sup>)</i></center>}}{\eqn{y~N(\beta X + b Z, \sigma^2)}}
#'
#'
#' where \eqn{y} is the controls' dependent variable, \eqn{X} the contrast
#' matrix for Population-level (or Fixed)
#' Effects, and \eqn{\beta} are the unknown coefficients to be estimate. \eqn{Z}
#' is the contrast matrix for the
#' Varying (or Random, or Group-level) effects, and \eqn{b} are the unknown
#' estimates for the varying effects.
#' \eqn{\sigma^2} is the variance.
#'
#' In order to estimate the coefficients of the Single Case, the formula is the following:
#'
#' \ifelse{html}{\out{<center><i>y<sub>pt</sub>~N(&phi; X<sub>pt</sub>, &sigma;<sup>2</sup><sub>pt</<sub>)</i></center>}}{\eqn{y_{pt}~N(\phi X_{pt}, \sigma_{pt}^2)}}
#'
#' where \eqn{\phi = \beta + \delta}.
#'
#' The validation of the approach can be found here:
#' https://www.doi.org/10.31234/osf.io/sajdq
#'
#' @section Details:
#' The main function of \strong{bmscstan} is \code{\link{BMSC}}, which uses formula syntax to
#' specify your model.
#'
#' @import rstan logspline bayesplot LaplacesDemon stats ggplot2
#'
#' @docType package
#' @name bmscstan
NULL

#' Data from a Single Case with brachial plexious lesion
#'
#' A dataset containing the results from the Body Sidedness Task
#' from a single Single Case
#'
#' @format A data frame with 467 rows and 4 variables
#' \describe{
#'  \item{RT}{Reaction times, in milliseconds}
#'  \item{Body.District}{Body district, categorial factor of
#'        Body Sidedness Task: FOOT or HAND}
#'  \item{Congruency}{The trail was Congruent or Incongruent?}
#'  \item{Side}{The trial showed a left or right limb}
#' }
"data.pt"

#' Data from a control group of 16 participants
#'
#' A dataset containing the results from the Body Sidednedd Task
#' from a control group of 16 participants
#'
#' @format A data frame with 4049 rows and 5 variables
#' \describe{
#'  \item{RT}{Reaction times, in milliseconds}
#'  \item{Body.District}{Body district, categorial factor of
#'        Body Sidedness Task: FOOT or HAND}
#'  \item{Congruency}{The trail was Congruent or Incongruent?}
#'  \item{Side}{The trial showed a left or right limb}
#'  \item{ID}{The participant ID}
#' }
"data.ctrl"

#' Fit Bayesian Multilevel Single Case models
#'
#' \code{BMSC} fits the Bayesian Multilevel Single Case models.
#'
#'
#' @param formula An object of class \code{formula}: a symbolic description of the model to be fitted.
#' @param data_ctrl An object of class \code{data.frame} (or one that can be coerced to that class)
#' containing data of all variables used in the model for the control group.
#' @param data_sc An object of class \code{data.frame} (or one that can be coerced to that class)
#' containing data of all variables used in the model for the Single Case
#' @param cores The number of cores to use when executing the Markov chains in parallel. The default is 1.
#' @param chains Number of Markov chains (defaults to 4).
#' @param iter Number of total iterations per chain (including warmup; defaults to 4000).
#' @param warmup A positive integer specifying number of warmup (aka burnin) iterations.
#' This also specifies the number of iterations used for stepsize adaptation,
#' so warmup samples should not be used for inference. The number of warmup should
#' not be larger than iter and the default is 2000.
#' @param seed The seed for random number generation to make results reproducible.
#' If NA (the default), Stan will set the seed randomly.
#' @param typeprior Set the desired prior distribution for the fixed effects.
#' \describe{
#' \item{normal}{a normal distribution with \eqn{\mu = 0} and \eqn{\sigma = 10}}
#' \item{cauchy}{a cauchy distribution with \eqn{\mu = 0} and scale \eqn{\sqrt{2}/2}}
#' \item{student}{a Student's T distribution, with \eqn{\mu = 0}, \eqn{\nu = 3} and \eqn{\sigma = 10}}
#' }
#' The normal distribution is the default.
#' @param ... further arguments to be passed to \strong{stan} function.
#'
#' @examples
#'  \donttest{
#'
#'  data(BSE)
#'
#' # Linear regression of data coming from a body representation paradigm
#' # with a control sample of 12 participants and one patient with
#' # unilateral brachial plexus lesion
#' mdl <- BMSC(formula = RT ~ Body.District * Congruency * Side +
#'             (Body.District + Congruency + Side | ID),
#'             data_ctrl = data.ctrl,
#'             data_sc = data.pt,
#'             cores = 4)
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
#' mdl.reg <- BMSC(y ~ x, data_ctrl = dat.ctrl, data_sc = dat.pt, seed = 10)
#'
#' # summarize the data
#'
#' summary(mdl.reg)
#' }
#'
#' @return a \code{BMSC} object
#'
#' @export
BMSC <- function(formula, data_ctrl, data_sc,
                cores = 1, chains = 4, warmup = 2000,
                iter = 4000, seed = NA, typeprior = "normal",
                ...){

  mypaste <- function(x){
    out <- NULL
    if(length(x)>1){
      for(xx in 1:(length(x)-1)) out <- paste(out,x[xx],"+")
    }
    out <- paste(out,x[length(x)])
    return(out)
  }

  if(missing(formula)) stop("the argument \"formula\" is not specified")
  if(missing(data_ctrl)) stop("the dataframe \"data_ctrl\" is not specified")
  if(missing(data_sc)) stop("the dataframe \"data_sc\" is not specified")
  if(typeprior!="normal"&&typeprior!="cauchy"&&typeprior!="student")
      stop("Not a valid typeprior")

  # extract formula's terms
  form.terms      <- attributes(terms(formula))$term.labels

  # build contrasts matrices of fixed effects
  fix.terms       <- form.terms[!(grepl("\\|",form.terms))]

  fix.formula     <- paste0(" ~",mypaste(fix.terms))

  matrix.fix.ctrl <- model.matrix(as.formula(fix.formula),data_ctrl)

  matrix.fix.pt   <- model.matrix(as.formula(fix.formula),data_sc)

  # build contrasts matrices for random effects
  ran.terms       <- form.terms[(grepl("\\|",form.terms))]

  ran.matrices    <- list()
  grouping        <- list()
  for(ran in ran.terms){
    tmp <- unlist(strsplit(ran,"\\|"))
    grouping[[ran]] <- trimws(tmp[2])
    ran.matrices[[ran]] <- model.matrix(as.formula(paste0(" ~",mypaste(tmp[1]))),data_ctrl)
  }

  stancode <- .building.model(ran.matrices,typeprior)

  datalist <- .building.data.list(ran.matrices,grouping,matrix.fix.ctrl,
                                 matrix.fix.pt,data_ctrl,data_sc,formula)

  mdl <- suppressMessages(stan(model_code = stancode, data = datalist, iter = iter,
             chains = chains,cores = cores, warmup = warmup,
             seed = seed, ...))

  out <- list(formula,mdl,data_sc,data_ctrl,datalist,stancode,typeprior)

  class(out) <- append(class(out),"BMSC")

  return(out)
}

