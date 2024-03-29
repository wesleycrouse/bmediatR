% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmediatR.R
\name{return_preset_odds_index}
\alias{return_preset_odds_index}
\title{Column indeces for commonly used posterior odds}
\usage{
return_preset_odds_index(
  odds_type = c("mediation", "partial", "complete", "colocal", "mediation_or_colocal",
    "y_depends_x", "reactive", "y_depends_m")
)
}
\arguments{
\item{odds_type}{The desired posterior odds.}
}
\value{
\code{return_preset_odds_index} returns a list with indices for individual models and combinations of models.
}
\description{
This helper function returns the columns of the log posterior case probabilities to be summed for
commonly desired log posterior odds summaries.
}
\examples{
return_preset_odds_index()
}
