#' Pairwise contrasts
#'
#' Calculate pairwise comparisons between marginal posterior distributions divided by group levels
#'
#'
#' @param mdl An object of class \code{BMSC}.
#' @param contrast Character value giving the name of the coefficient whose levels need to be compared.
#' @param covariate at the moment is silent
#' @param who parameter to choose the estimates to contrast
#' \describe{
#'         \item{control}{only the controls}
#'         \item{singlecase}{only the single case \eqn{(\beta + \delta)}}
#'         \item{delta}{only the difference between the single case and controls}
#' }
#'
#' @examples
#'  \donttest{
#'
#' ######################################
#' # simulation of controls' group data
#' ######################################
#'
#' # Number of levels for each condition and trials
#' NCond1  <- 2
#' NCond2  <- 2
#' Ntrials <- 8
#' NSubjs  <- 30
#'
#' betas <- c( 0 , 0 , 0 ,  0.2)
#'
#' data.sim <- expand.grid(
#'   trial      = 1:Ntrials,
#'   ID         = factor(1:NSubjs),
#'   Cond1      = factor(1:NCond1),
#'   Cond2      = factor(1:NCond2)
#' )
#'
#' contrasts(data.sim$Cond1) <- contr.sum(2)
#' contrasts(data.sim$Cond2) <- contr.sum(2)
#'
#' ### d.v. generation
#' y <- rep( times = nrow(data.sim) , NA )
#'
#' # cheap simulation of individual random intercepts
#' set.seed(1)
#' rsubj <- rnorm(NSubjs , sd = 0.1)
#'
#' for( i in 1:length( levels( data.sim$ID ) ) ){
#'
#'   sel <- which( data.sim$ID == as.character(i) )
#'
#'   mm  <- model.matrix(~ Cond1 * Cond2 , data = data.sim[ sel , ] )
#'
#'   set.seed(1 + i)
#'   y[sel] <- mm %*% as.matrix(betas + rsubj[i]) +
#'     rnorm( n = Ntrials * NCond1 * NCond2 )
#'
#' }
#'
#' data.sim$y <- y
#'
#' # just checking the simulated data...
#' boxplot(y~Cond1*Cond2, data = data.sim)
#'
#' ######################################
#' # simulation of patient data
#' ######################################
#'
#' betas.pt <- c( 0 , 0.8 , 0 ,  0)
#'
#' data.pt <- expand.grid(
#'   trial      = 1:Ntrials,
#'   Cond1      = factor(1:NCond1),
#'   Cond2      = factor(1:NCond2)
#' )
#'
#' contrasts(data.pt$Cond1) <- contr.sum(2)
#' contrasts(data.pt$Cond2) <- contr.sum(2)
#'
#' ### d.v. generation
#' mm  <- model.matrix(~ Cond1 * Cond2 , data = data.pt )
#'
#' set.seed(5)
#' data.pt$y <- (mm %*% as.matrix(betas.pt) +
#'                 rnorm( n = Ntrials * NCond1 * NCond2 ))[,1]
#'
#' # just checking the simulated data...
#' boxplot(y~Cond1*Cond2, data = data.pt)
#'
#' mdl <- BMSC(y ~ Cond1 * Cond2 + ( 1 | ID ),
#'             data_ctrl = data.sim, data_sc = data.pt, seed = 77,
#'             typeprior = "cauchy", s = 1)
#'
#' summary(mdl)
#'
#' pp_check(mdl)
#'
#'
#' pairwise.BMSC( mdl, contrast = "Cond11:Cond21")
#'
#' }
#'
#' @return a \code{pairwise.BMSC} object
#'
#' @export
pairwise.BMSC = function(mdl, contrast, covariate = NULL, who = "delta") {
  # function to find all the column of the posterior distribution involved in the contrast
  check_contrasts <- function( contr.names , contr.parts , contrast ){
    M <- strsplit(contr.names , ":")
    out <- NULL
    len <- NULL
    i <- 1

    ## in this loop I store in len the number of main effects present in every coefficients (main effects and interactions)
    ## and I check, for each coefficient, how many of the main effects involved in the contrast are present, storing it in out
    for(m in M){
      len <- c(len , length(m))
      for(cp in contr.parts){
        for(mm in m){
          if(grepl(cp,mm)) out <- c(out , i)
        }
      }
      i <- i +1
    }

    ## table for out
    tab.out <- table( out )

    ## return the indexes of all the coefficients involved in the contrast
    ## if the total number of main effects (in len)
    return(as.numeric(names(which(len[unique(out)] == tab.out))))
  }

  se <- function(object) {
    sd(object)/sqrt(length(object))
  }

  if (mdl[[7]] == "normal") {
    d0 <- dnorm(0, 0, 10)
  } else if (mdl[[7]] == "cauchy") {
    d0 <- dcauchy(0, 0, sqrt(2)/2)
  } else if (mdl[[7]] == "student") {
    d0 <- LaplacesDemon::dst(0, 10, 3)
  }


  if (class(mdl)[2] != "BMSC")
    stop("Not a valid BMSC object.")

  if (!exists("contrast"))
    stop("Not a valid contrast")

  # variable declaration
  tmp.post <- tmp.data <- contr.names <- contr.column <- contr.parts <- contr.table <-
    tmp.marginal <- findRow <- sum_logspl <- bf_sd <- tmp.y <- NULL

  # extracting data from BMSC model
  if(who == "delta"){
    tmp.post <- rstan::extract(mdl[[2]], pars = "b_Delta")
    tmp.data <- mdl[[3]]
    contr.names <- colnames(mdl[[5]]$XF_Pts)
    contr.table <- mdl[[5]]$XF_Pts
  } else if(who == "control"){
    tmp.post <- rstan::extract(mdl[[2]], pars = "b_Ctrl")
    tmp.data <- mdl[[4]]
    contr.names <- colnames(mdl[[5]]$XF_Ctrl)
    contr.table <- mdl[[5]]$XF_Ctrl
  } else if(who == "singlecase"){
    delta <- rstan::extract(mdl[[2]], pars = "b_Delta")
    ctrl  <- rstan::extract(mdl[[2]], pars = "b_Ctrl")
    tmp.post <- list(delta[[1]] + ctrl[[1]])
    tmp.data <- mdl[[3]]
    contr.names <- colnames(mdl[[5]]$XF_Pts)
    contr.table <- mdl[[5]]$XF_Pts
  } else stop("Not a valid \"who\" value")

  if(sum(grepl(contrast,contr.names))==0) stop("Not a valid contrast")

  # find the columns of the posterior distribution
  contr.parts <- unlist(strsplit(contrast,":"))

  contr.column <- which(grepl(contrast,contr.names))

  contr.column <- c(contr.column, check_contrasts(contr.names[1:(contr.column-1)] , contr.parts , contrast))

  contr.column <- c(contr.column, which(grepl("(Intercept)",contr.names)))

  contr.column <- contr.column[order(contr.column)]

  # find a name to each marginal distribution

  create.names <- matrix(nrow = nrow(tmp.data) , ncol = length(contr.parts))
  ricontrasts  <- matrix(nrow = nrow(tmp.data) , ncol = length(contr.parts))

  for( i in 1:length(contr.parts) ){
    cp <- contr.parts[i]
    create.names[ , i] <- as.character(tmp.data[,grepl(substr(cp , start = 1 , stop = (nchar(cp)-1) ), colnames(tmp.data) )])
    ricontrasts[ , i]  <- contr.table[,contr.column][,cp == colnames(contr.table[,contr.column])]
  }

  create.names <- as.data.frame( cbind( apply(create.names, 1, paste, collapse = " ") , ricontrasts ) )

  create.names <- unique(create.names)

  # compute the marginal distributions

  marginal_distribution <- list()

  tmp.table <- unique(contr.table[,contr.column])

  for(i in 1:nrow(create.names)){
    findRow <- NULL
    for(j in 1:nrow(tmp.table)){
      if(sum(tmp.table[j,colnames(tmp.table)%in%contr.parts] == create.names[i,2:(length(contr.parts)+1)])==length(contr.parts)) findRow <- j
    }
    marginal_distribution[[create.names[i,1]]] <- data.frame( y = colSums(t(tmp.post[[1]][,contr.column])*tmp.table[findRow,]),
                                            name = create.names[i,1])
  }


  # create a summary for each marginal distribution

  sum.unique <- list()

  for(marginal_name in marginal_distribution){
    sum_logspl <- .suppresslogspline(marginal_name$y)
    bf_sd <- d0/.suppressdlogspline(0, sum_logspl)

    sum.unique[[marginal_name$name[1]]] <- as.data.frame(cbind(mean(marginal_name$y),
                                                               se(marginal_name$y),
                                                               sd(marginal_name$y),
                                                               quantile(marginal_name$y, probs = 2.5/100),
                                                               quantile(marginal_name$y, probs = 25/100),
                                                               quantile(marginal_name$y, probs = 50/100),
                                                               quantile(marginal_name$y, probs = 75/100),
                                                               quantile(marginal_name$y, probs = 97.5/100),
                                                               bf_sd))

    colnames(sum.unique[[marginal_name$name[1]]]) <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%","BF10 (not zero)")
  }

  sum.unique <- do.call("rbind" , sum.unique)

  row.names(sum.unique) <- create.names[,1]

  # create a summary for each comparison

  sum.contrast <- list()
  sum.names    <- list()

  for(i in 1:(length(marginal_distribution)-1)){
    for(j in (i+1):length(marginal_distribution)){
      tmp.y <- marginal_distribution[[i]]$y - marginal_distribution[[j]]$y

      sum_logspl <- .suppresslogspline(tmp.y)
      bf_sd <- d0/.suppressdlogspline(0, sum_logspl)

      sum.contrast[[paste(i,j)]] <- as.data.frame(cbind(mean(tmp.y),
                                                                 se(tmp.y),
                                                                 sd(tmp.y),
                                                                 quantile(tmp.y, probs = 2.5/100),
                                                                 quantile(tmp.y, probs = 25/100),
                                                                 quantile(tmp.y, probs = 50/100),
                                                                 quantile(tmp.y, probs = 75/100),
                                                                 quantile(tmp.y, probs = 97.5/100),
                                                                 bf_sd))

      colnames(sum.contrast[[paste(i,j)]]) <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%","BF10")

      sum.names[[paste(i,j)]] <- paste(marginal_distribution[[i]]$name[1],
                                                 marginal_distribution[[j]]$name[1],
                                                 sep = " - ")
    }
  }

  sum.contrast <- do.call("rbind" , sum.contrast)

  row.names(sum.contrast) <- do.call("c" , sum.names)

  out <- list(sum.unique , sum.contrast , contrast , covariate , marginal_distribution , mdl[[7]])

  class(out) <- append(class(out),"pairwise.BMSC")

  return(out)

}

#' Print summaries of Pairwise Bayesian Multilevel Single Case objects
#'
#'
#' @param x An object of class \code{pairwise.BMSC}, resulting from the \link{pairwise.BMSC} function.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @method print pairwise.BMSC
#' @export
print.pairwise.BMSC = function(x, ...) {

  if(is.null(x[[4]])){
    cat("\nPairwise Bayesian Multilevel Single Case contrasts of coefficients divided by", x[[3]] , "\n\n")
  } else {
    cat("\nPairwise Bayesian Multilevel Single Case contrasts of", x[[4]] ,"covariate divided by", x[[3]] , "\n\n")
  }

  cat("\n\n  Marginal distributions\n\n")

  print(x[[1]], ...)
  cat("\n")
  cat("\n\n  Table of contrasts\n\n")

  print(x[[2]], ...)
}
