\name{bmscstanNews}
\title{bmscstan News}
\encoding{UTF-8}

\section{Version 1.2.0.1}{
  \itemize{
    \item solved small, but annoying, bugs
    \item vignette updated
  }
}

\section{Version 1.2.0}{
  \itemize{
    \item added code for Pareto smoothed importance sampling (PSIS),
          Leaving-One-Out (LOO) and
          Widely applicable information criterion (WAIC)
    \item added function loo.BMSC
    \item added function loo_compare.loo_BMSC
    \item added function print.loo_BMSC
    \item added function print.loo_compare_BMSC
    \item added diagnostics for PSIS
    \item added functions for waic
    \item modfiied main functions for binomial data
    \item modified summary.BMSC function for binomial data
    \item added the possibility to fit only-intercept fixed effects models
    \item vignette updated
  }
}

\section{Version 1.1.1}{
  \itemize{
    \item solved small bugs in pairwise.BMSC
  }
}

\section{Version 1.1.0}{
  \itemize{
    \item added dispersion parameter "s" to \code{BMSC} function.
    \item Updated \code{summary.BMSC} function to take into account the "s" parameter
    \item Updated all examples.
    \item Updated documentation.
  }
}

\section{Version 1.0.0}{
  \itemize{
    \item New package bmscstan, for fitting Bayesian Multilevel Single Case models.
  }
}
