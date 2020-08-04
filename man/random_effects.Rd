\name{randomeffects}
\alias{randomeffect}
\title{Random Effects specification on
       Bayesian Multilevel Single Case models using 'Stan'}
\description{
  The \strong{BMSC} function allows the flexibility of multilevel (generalised)
  linear models on single case analysis.

  In particular, it is possible to specify the population-level (a.k.a. mixed
  effects) and the group-level (a.k.a. random effects) coefficients.

  The specification of the population- and group-level effects can be done
  using the well-known \strong{lme4} notation with specific limitations:

  \itemize{
  \item it is no possible to estimate uncorrelated group-level effects
  \item it is no possible to directly estimate nested effects. You need to use
    a trick that is specified in the \strong{Details} section.
  }
}
\details{

\tabular{lc}{
  lmer formulation                                \tab BMSC availability   \cr
  \code{(1 | grouping_factor)}                    \tab Yes                 \cr
  \code{(1 + slope | grouping_factor)}            \tab Yes                 \cr
  \code{(0 + slope | grouping_factor)}            \tab No                  \cr
  \code{(1 | grouping_factor1 : grouping_factor2)}\tab Yes[^1]             \cr
  \code{(1 | grouping_factor1 / grouping_factor2)}\tab Yes[^2]             \cr
 }

[^1]: The \strong{BMSC} function dose not allow the use of the interaction
  symbol ":", but this problem is easily solved by creating a new variable
  within your dataframe given by the interaction of the two factors.

[^2]: The \code{(1 | grouping_factor1 / grouping_factor2)  } syntax is the
  equivalent of the explicit version
  \code{(1 \| grouping_factor1:grouping_factor2) +  (1 | grouping_factor1)}.

  Therefore, you need to create a new grouping factor representing the
  interaction between \code{grouping_factor1} and \code{grouping_factor2},
  and use this in the explicit version
  \code{(1 | grouping_factor_interaction) +  (1 | grouping_factor1)}.
}
