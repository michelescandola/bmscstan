#' BMSC: Bayesian Multilevel Single Case models using 'Stan'
#'
#' The \strong{BMSC} package provides an interface to fit Bayesian Multilevel Single Case models.
#' These models compare the performance of a single patient against a control group, combining
#' the flexibility of multilevel models and the potentiality of Bayesian Statistics.
#'
#' The package is now limited to gaussian data only, but we will further expand it to cover
#' binomial and ordinal (Likert scales) data.
#'
#' By means of BMSC the effects of the control group and the effects of the deviance between the
#' patient and the group will be esimated.
#'
#' The model to estimate the controls parameters is:
#'
#' \deqn{y~N(\beta X + b Z, \sigma^2)}
#'
#' where \eqn{y} is the controls' dependent variable, \eqn{X} the contrast matrix for Population-level (or Fixed)
#' Effects, and \eqn{\beta} are the unknown coefficients to be estimate. \eqn{Z} is the contrast matrix for the
#' Varying (or Random, or Group-level) effects, and \eqn{b} are the unknown estimates for the varying effects.
#' \eqn{\sigma^2} is the variance.
#'
#' In order to estimate the coefficients of the patient, the formula is the following:
#'
#' \deqn{y_{pt}~\mathcal{N}(\phi X_{pt}, \sigma_{pt}^2)}
#'
#' where \eqn{\phi = \beta + \delta}.
#'
#' @section Details:
#' The main function of \strong{BMSC} is \code{\link{BMSC}}, which uses formula syntax to
#' specify your model.
#'
#' @import rstan logspline bayesplot LaplacesDemon stats ggplot2
#'
#' @docType package
#' @name BMSC-package
NULL

#' Data from a patient with brachial plexious lesion
#'
#' A dataset containing the results from the Body Sidednedd Task
#' from a single patient
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
#' @param data_pt An object of class \code{data.frame} (or one that can be coerced to that class)
#' containing data of all variables used in the model for the patient.
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
#'  \dontrun{
#'
#'  data(BSE)
#'
#'  # Normal robust regression of data coming from a body representation paradigm
#'  # with a control sample of 12 participants and one patient with
#'  # unilateral brachial plexus lesion
#'  mdl <- BMSC(formula = RT ~ Body.District * Congruency * Side + (Body.District + Congruency + Side | ID),
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
#'  }
#'
#' @return a \code{BMSC} object
#'
#' @export
BMSC <- function(formula, data_ctrl, data_pt,
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
  if(missing(data_pt)) stop("the dataframe \"data_pt\" is not specified")
  if(typeprior!="normal"&&typeprior!="cauchy"&&typeprior!="student")
    stop("Not a valid typeprior")

  # extract formula's terms
  form.terms      <- attributes(terms(formula))$term.labels

  # build contrasts matrices of fixed effects
  fix.terms       <- form.terms[!(grepl("\\|",form.terms))]

  fix.formula     <- paste0(" ~",mypaste(fix.terms))

  matrix.fix.ctrl <- model.matrix(as.formula(fix.formula),data_ctrl)

  matrix.fix.pt   <- model.matrix(as.formula(fix.formula),data_pt)

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
                                 matrix.fix.pt,data_ctrl,data_pt,formula)

  mdl <- stan(model_code = stancode, data = datalist, iter = iter,
             chains = chains,cores = cores, warmup = warmup,
             seed = seed, ...)

  out <- list(formula,mdl,data_pt,data_ctrl,datalist,stancode,typeprior)

  class(out) <- append(class(out),"BMSC")

  return(out)
}


.building.model <- function(ran.matrices=NULL,typeprior){
  code.data <- "
  data {
    int<lower=1> Obs_Controls; // the total number of observations in the control group
    int<lower=1> Obs_Patients; // the total numebr of observation in the patient
    int<lower=1> Nparameters;  // the total number of parameters for the independent variables
    real y_Ctrl[Obs_Controls];                  // the dependent variable for the control group
    real y_Pts[Obs_Patients];                   // the patient d.v.
    matrix[Obs_Controls,Nparameters] XF_Ctrl;   // the control matrix
    matrix[Obs_Patients,Nparameters] XF_Pts;    // the patient matrix"

  code.parameter <-   "parameters {
    vector[Nparameters] b_Ctrl;             //the regression parameters for controls
    vector[Nparameters] b_Delta;            //the regression parameters for the controls - patients difference
    real<lower=0> sigmaC;                    //the standard deviation for controls
    real<lower=0> sigmaP;                    //the standard deviation for patient"

  code.transformed.parameter <-   "transformed parameters{
    real mu_Pts[Obs_Patients];
    real mu_Ctrl[Obs_Controls];"

  last.string.code.transformed.parameter <- "

    for(i in 1:Obs_Patients){
      mu_Pts[i] = dot_product(b_Ctrl+b_Delta,XF_Pts[i,]);
    }
    for(i in 1:Obs_Controls){
      mu_Ctrl[i] = dot_product(b_Ctrl,XF_Ctrl[i,])"
  if(typeprior=="normal"){
    code.model <-   "

      target += cauchy_lpdf(sigmaC|0,1000);
      target += cauchy_lpdf(sigmaP|0,1000);

      target += normal_lpdf(b_Ctrl  | 0, 10);
      target += normal_lpdf(b_Delta | 0, 10);

      target += cauchy_lpdf(y_Pts|mu_Pts,sigmaP);
      target += cauchy_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
    }"
  }else if(typeprior=="cauchy"){
    code.model <-   "

      target += cauchy_lpdf(sigmaC|0,1000);
      target += cauchy_lpdf(sigmaP|0,1000);

      target += cauchy_lpdf(b_Ctrl  | 0, sqrt(2)/2);
      target += cauchy_lpdf(b_Delta | 0, sqrt(2)/2);

      target += cauchy_lpdf(y_Pts|mu_Pts,sigmaP);
      target += cauchy_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
    }"
  }else if(typeprior=="student"){
    code.model <-   "

      target += cauchy_lpdf(sigmaC|0,1000);
      target += cauchy_lpdf(sigmaP|0,1000);

      target += student_t_lpdf(b_Ctrl  | 3, 0, 10);
      target += student_t_lpdf(b_Delta | 3, 0, 10);

      target += cauchy_lpdf(y_Pts|mu_Pts,sigmaP);
      target += cauchy_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
    }"
  }

  code.generated.quantities <-   "generated quantities {
    real y_pt_rep[Obs_Patients];
    real y_ct_rep[Obs_Controls];

    for(i in 1:Obs_Patients){
      y_pt_rep[i] = cauchy_rng(mu_Pts[i], sigmaP);
    }
    for(i in 1:Obs_Controls){
      y_ct_rep[i] = cauchy_rng(mu_Ctrl[i], sigmaC);
    }
  }"


  if(!is.null(ran.matrices)){
    ir <- 1
    for(ran in ran.matrices){
      if(ncol(ran)>1){
        code.data <- paste(code.data,
                          paste0("    int<lower=1> Nrands",ir,";    //number of random coefficients for grouping factor ",ir),
                          paste0("    matrix[Obs_Controls,Nrands",ir,"] XR_Ctrl",ir,";    //the control random matrix for grouping factor ",ir),
                          paste0("    int grouping",ir,"[Obs_Controls];    //the index vector for the grouping factor ",ir),
                          paste0("    int<lower=1> Ngrouping",ir,";    // the total number of levels for grouping",ir),
                          sep ="\n"
                          )

        code.parameter <- paste(code.parameter,
                               paste0("    vector<lower=0>[Nrands",ir,"] sigma_u",ir,";      // random effects sd for grouping factor ",ir),
                               paste0("    cholesky_factor_corr[Nrands",ir,"] L_Omega",ir,";"),
                               paste0("    matrix[Nrands",ir,",Ngrouping",ir,"] z_u",ir,";"),
                               sep="\n"
                               )

        code.transformed.parameter <- paste(code.transformed.parameter,
                                           paste0("    matrix[Nrands",ir,",Ngrouping",ir,"] u",ir,";"),
                                           paste0("    u",ir," = (diag_pre_multiply(sigma_u",ir,", L_Omega",ir,") * z_u",ir,"); //random effects for grouping factor",ir),
                                           sep="\n")

        last.string.code.transformed.parameter <- paste(last.string.code.transformed.parameter,
                                                       paste0("+ dot_product(u",ir,"[,grouping",ir,"[i]],XR_Ctrl",ir,"[i,])"))

        code.model <- paste(paste0("target += lkj_corr_cholesky_lpdf(L_Omega",ir," | 1);"),
                           paste0("target += normal_lpdf(to_vector(z_u",ir,") | 0, 1);"),
                           code.model,sep="\n")

      }else if(dim(ran)[2]==1){
        code.data <- paste(code.data,
                          paste0("    int grouping",ir,"[Obs_Controls];    //the index vector for the grouping factor ",ir),
                          paste0("    int<lower=1> Ngrouping",ir,";    // the total number of levels for grouping",ir),
                          sep ="\n"
        )

        code.parameter <- paste(code.parameter,
                               paste0("    real<lower=0> sigma_u",ir,";      // random effects sd for grouping factor ",ir),
                               paste0("    real[Ngrouping",ir,"] u",ir,";"),
                               sep="\n"
        )

        last.string.code.transformed.parameter <- paste(last.string.code.transformed.parameter,
                                                       paste0("+ u",ir,"[grouping",ir,"[i]]"))

        code.model <- paste(paste0("target += normal_lpdf(u",ir," | 0, 10);"),
                           code.model,sep="\n")
      }

      ir <- ir +1
    }
  }

  code.transformed.parameter <- paste(code.transformed.parameter,
                                     paste0(last.string.code.transformed.parameter,";"),
                                     "    }",sep="\n")

  code.model <- paste("  model{
    //priors",code.model,sep="\n")

  out <- paste(code.data,"  }",
            code.parameter,"  }",
            code.transformed.parameter,"  }",
            code.model,
            code.generated.quantities,sep="\n")

  return(out)
}

.building.data.list <- function(ran.matrices = NULL, grouping,
                               matrix.fix.ctrl, matrix.fix.pt,
                               data_ctrl, data_pt, formula){
  data.list <- list(
    Nparameters = ncol(matrix.fix.ctrl),

    y_Ctrl=data_ctrl[,as.character(formula[2])],
    y_Pts =data_pt[,as.character(formula[2])],

    XF_Ctrl=matrix.fix.ctrl,
    XF_Pts =matrix.fix.pt,

    Obs_Controls = nrow(matrix.fix.ctrl),
    Obs_Patients = nrow(matrix.fix.pt)
  )

  if(!is.null(ran.matrices)){
    ir <- 1
    for(ran in ran.matrices){

      nn <- length(data.list)


      if(ncol(ran)>1){

        data.list[[paste0("Nrands",ir)]]    <- ncol(ran)
        data.list[[paste0("XR_Ctrl",ir)]]   <- ran
        data.list[[paste0("grouping",ir)]]  <- as.numeric(data_ctrl[,grouping[[ir]]])
        data.list[[paste0("Ngrouping",ir)]] <- length(unique(data_ctrl[,grouping[[ir]]]))

      }else if(dim(ran)[2]==1){
        data.list[[paste0("grouping",ir)]]  <- as.numeric(data_ctrl[,grouping[[ir]]])
        data.list[[paste0("Ngrouping",ir)]] <- length(unique(data_ctrl[,grouping[[ir]]]))
      }
      ir <- ir + 1
    }
  }

  return(data.list)
}


#' Computes log marginal likelihood via bridge sampling.
#'
#' @param object a BMSC object
#' @param ... further arguments passed to or from other methods.
#' @return an "psis_loo" "loo" object
#' @export

bridge_sampler.BMSC <- function(object, ...){

  if (requireNamespace("bridgesampling", quietly = TRUE)) {
    ans <- bridgesampling::bridge_sampler(object[[2]],...)
  }

  return(ans)

}

