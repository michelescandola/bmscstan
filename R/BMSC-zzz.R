.suppresslogspline <- function(x){
  suppressWarnings(logspline::logspline(x))
}

.suppressdlogspline <- function(d,x){
  suppressWarnings(logspline::dlogspline(d,x))
}

.building.model <- function(ran.matrices=NULL,typeprior,s,family){
  if(family == "gaussian"){
    code.data <- "
    data {
      int<lower=1> Obs_Controls; // the total number of observations in the control group
      int<lower=1> Obs_Patients; // the total numebr of observation in the patient
      int<lower=1> Nparameters;  // the total number of parameters for the independent variables
      real y_Ctrl[Obs_Controls];                  // the dependent variable for the control group
      real y_Pts[Obs_Patients];                   // the patient d.v.
      matrix[Obs_Controls,Nparameters] XF_Ctrl;   // the control matrix
      matrix[Obs_Patients,Nparameters] XF_Pts;    // the patient matrix
      real s;                    // the dispersion parameter for priors"

    code.parameter <-   "parameters {
      vector[Nparameters] b_Ctrl;             //the regression parameters for controls
      vector[Nparameters] b_Delta;            //the regression parameters for the controls - patients difference
      real<lower=0> sigmaC;                    //the standard deviation for controls
      real<lower=0> sigmaP;                    //the standard deviation for patient"

    code.transformed.parameter <-   "transformed parameters{
      real mu_Pts[Obs_Patients];
      real mu_Ctrl[Obs_Controls];"

    code.transformed.parameter2 <- ""

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

        target += normal_lpdf(b_Ctrl  | 0, s);
        target += normal_lpdf(b_Delta | 0, s);

        target += normal_lpdf(y_Pts|mu_Pts,sigmaP);
        target += normal_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
      }"
    }else if(typeprior=="cauchy"){
      code.model <-   "

        target += cauchy_lpdf(sigmaC|0,1000);
        target += cauchy_lpdf(sigmaP|0,1000);

        target += cauchy_lpdf(b_Ctrl  | 0, s);
        target += cauchy_lpdf(b_Delta | 0, s);

        target += normal_lpdf(y_Pts|mu_Pts,sigmaP);
        target += normal_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
      }"
    }else if(typeprior=="student"){
      code.model <-   "

        target += cauchy_lpdf(sigmaC|0,1000);
        target += cauchy_lpdf(sigmaP|0,1000);

        target += student_t_lpdf(b_Ctrl  | 3, 0, s);
        target += student_t_lpdf(b_Delta | 3, 0, s);

        target += normal_lpdf(y_Pts|mu_Pts,sigmaP);
        target += normal_lpdf(y_Ctrl|mu_Ctrl,sigmaC);
      }"
    }

    code.generated.quantities <-   "generated quantities {
      real y_pt_rep[Obs_Patients];
      real y_ct_rep[Obs_Controls];
      vector[Obs_Patients] log_lik_pt;
      vector[Obs_Controls] log_lik_ct;

      for(i in 1:Obs_Patients){
        y_pt_rep[i] = normal_rng(mu_Pts[i], sigmaP);
        log_lik_pt[i] = normal_lpdf(y_Pts[i] | mu_Pts[i], sigmaP);
      }
      for(i in 1:Obs_Controls){
        y_ct_rep[i] = normal_rng(mu_Ctrl[i], sigmaC);
        log_lik_ct[i] = normal_lpdf(y_Ctrl[i] | mu_Ctrl[i], sigmaC);
      }
    }"
  } else if (family == "binomial"){
    code.data <- "
    data {
      int<lower=1> Obs_Controls; // the total number of observations in the control group
      int<lower=1> Obs_Patients; // the total numebr of observation in the patient
      int<lower=1> Nparameters;  // the total number of parameters for the independent variables
      int y_Ctrl[Obs_Controls];                  // the dependent variable for the control group
      int y_Pts[Obs_Patients];                   // the patient d.v.
      int n_Ctrl[Obs_Controls];
      int n_Pts[Obs_Patients];
      matrix[Obs_Controls,Nparameters] XF_Ctrl;   // the control matrix
      matrix[Obs_Patients,Nparameters] XF_Pts;    // the patient matrix
      real s;                    // the dispersion parameter for priors"

    code.parameter <-   "parameters {
      vector[Nparameters] b_Ctrl;             //the regression parameters for controls
      vector[Nparameters] b_Delta;            //the regression parameters for the controls - patients difference"

    code.transformed.parameter <-   "transformed parameters{
      real mu_Pts[Obs_Patients];
      real mu_Ctrl[Obs_Controls];"

    code.transformed.parameter2 <- ""

    last.string.code.transformed.parameter <- "

      for(i in 1:Obs_Patients){
        mu_Pts[i] = dot_product(b_Ctrl+b_Delta,XF_Pts[i,]);
      }
      for(i in 1:Obs_Controls){
        mu_Ctrl[i] = dot_product(b_Ctrl,XF_Ctrl[i,])"
    if(typeprior=="normal"){
      code.model <-   "

        target += normal_lpdf(b_Ctrl  | 0, s);
        target += normal_lpdf(b_Delta | 0, s);

        target += binomial_logit_lpmf(y_Pts|n_Pts,mu_Pts);
        target += binomial_logit_lpmf(y_Ctrl|n_Ctrl,mu_Ctrl);
      }"
    }else if(typeprior=="cauchy"){
      code.model <-   "

        target += cauchy_lpdf(b_Ctrl  | 0, s);
        target += cauchy_lpdf(b_Delta | 0, s);

        target += binomial_logit_lpmf(y_Pts|n_Pts,mu_Pts);
        target += binomial_logit_lpmf(y_Ctrl|n_Ctrl,mu_Ctrl);
      }"
    }else if(typeprior=="student"){
      code.model <-   "

        target += student_t_lpdf(b_Ctrl  | 3, 0, s);
        target += student_t_lpdf(b_Delta | 3, 0, s);

        target += binomial_logit_lpmf(y_Pts|n_Pts,mu_Pts);
        target += binomial_logit_lpmf(y_Ctrl|n_Ctrl,mu_Ctrl);
      }"
    }

    code.generated.quantities <-   "generated quantities {
      int y_pt_rep[Obs_Patients];
      int y_ct_rep[Obs_Controls];
      vector[Obs_Patients] log_lik_pt;
      vector[Obs_Controls] log_lik_ct;
      real tmp;

      for(i in 1:Obs_Patients){
        tmp = inv_logit(mu_Pts[i]);
        y_pt_rep[i] = binomial_rng(n_Pts[i], tmp);
        log_lik_pt[i] = binomial_logit_lpmf(y_Pts[i] |n_Pts[i], mu_Pts[i]);
      }
      for(i in 1:Obs_Controls){
        tmp = inv_logit(mu_Ctrl[i]);
        y_ct_rep[i] = binomial_rng(n_Ctrl[i], tmp);
        log_lik_ct[i] = binomial_logit_lpmf(y_Ctrl[i] |n_Ctrl[i], mu_Ctrl[i]);
      }
    }"
  }

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
                                            sep="\n")

        code.transformed.parameter2 <- paste(code.transformed.parameter2,
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
                                paste0("    real u",ir,"[Ngrouping",ir,"];"),
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
                                      code.transformed.parameter2,
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
                                data_ctrl, data_sc, formula, s,
                                family){

  if(family == "gaussian"){
    data.list <- list(
      Nparameters = ncol(matrix.fix.ctrl),

      y_Ctrl=data_ctrl[,as.character(formula[2])],
      y_Pts =data_sc[,as.character(formula[2])],

      XF_Ctrl=matrix.fix.ctrl,
      XF_Pts =matrix.fix.pt,

      Obs_Controls = nrow(matrix.fix.ctrl),
      Obs_Patients = nrow(matrix.fix.pt),

      s = s
    )
  } else if(family == "binomial"){
    dv = as.character(formula[2])
    dv = gsub("cbind\\(","",dv)
    dv = gsub("\\)","",dv)
    dv1 = unlist(strsplit(dv,","))[1]
    dv2 = trimws(unlist(strsplit(dv,","))[2])

    data.list <- list(
      Nparameters = ncol(matrix.fix.ctrl),

      y_Ctrl=data_ctrl[,dv1],
      y_Pts =data_sc[,dv1],

      n_Ctrl=data_ctrl[,dv2],
      n_Pts =data_sc[,dv2],

      XF_Ctrl=matrix.fix.ctrl,
      XF_Pts =matrix.fix.pt,

      Obs_Controls = nrow(matrix.fix.ctrl),
      Obs_Patients = nrow(matrix.fix.pt),

      s = s
    )
  }


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
