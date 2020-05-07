%% metaboGSE Version History

\name{NEWS}
\title{metaboGSE News}
\encoding{UTF-8}

\section{Changes in version 1.2.4 2020-05-07}{
  \itemize{
    \item Remove dependence on \code{GO.db}, \code{ctc}, \code{clpAPI}, \code{KernSmooth}.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.2.3 2020-01-03}{
  \itemize{
    \item Remove dependence on \code{topGO}.
    \item Move \code{GO.db} to Imports.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.2.2 2019-12-06}{
  \itemize{
    \item Change sys::eval_safe to unix::eval_safe
    \item Use \code{fixObjVal = F} in fluxVar
    \item Move \code{topGO} to Suggests.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.2.1 2018-12-12}{
  \itemize{
    \item Supply the list of metabolites in each submodel
    \item Supply the AUC difference in posthoc pairwise comparisons
    \item Rescue both subtrates and products of reactions in the
    objective function to reduce bugs of infinite loops in LP.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.2.0 2018-11-22}{
  \itemize{
    \item Introduce the option \code{timeout} to the function
    \code{rescue}, \code{fitness} and \code{submnet} to deal with the
    bug of infinite loop in \code{optimizeProb} and \code{fluxVar} from
    \code{sybil} with different solvers.
    \item Fix the function \code{bestRanking} to deal with the NA values
    obtained in infeasible LP.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.1.2 2018-10-16}{
  \itemize{
    \item Suggest in the vignette to improve the ranking for genes of
    identical expression in each sample with the overall expression in
    all the samples.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.1.1 2018-09-24}{
  \itemize{
    \item Rescue only substrates, but not products, of reactions in the
    objective function.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.1.0 2018-08-07}{
  \itemize{
    \item Modify the fitness formulation
    \item Introduce the function \code{fitness} to separate the
    evaluation of gene rankings and the construction of submodels with
    the best ranking.
    \item Restructure the function \code{metaboGSE}.
  }
}
% ---------------------------------------------------------------------------- %

\section{Changes in version 1.0.2 2018-04-13}{
  \itemize{
    \item Fix CRAN submission errors.
  }
}
% ---------------------------------------------------------------------------- %
