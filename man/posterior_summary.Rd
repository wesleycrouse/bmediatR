% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmediatR.R
\name{posterior_summary}
\alias{posterior_summary}
\title{Summaries of posterior model probabilities function}
\usage{
posterior_summary(ln_prob_data, ln_prior_c, c_numerator)
}
\arguments{
\item{ln_prob_data}{Log posterior likelihoods under the various models, returned by bmediatR().}

\item{ln_prior_c}{Log prior case probabilities. If posterior_summary() is being used for a non-default posterior odds
summary, the log prior case probabilities used with bmediatR() are stored in its output.}

\item{c_numerator}{The index of cases to be summed in the numerator of the posterior odds. Cases, their order, and likelihoods
are provided in model_info().}
}
\value{
\code{posterior_summary} returns a list containing the following components:
 \item{ln_post_c}{a matrix with posterior probabilities of each causal model for each candidate mediator.}
 \item{ln_post_odds}{a matrix with posterior odds of individual models or combinations of models for each candidate mediator.}
 \item{ln_prior_odds}{a single row matrix with prior odds of individual models or combinations of models.}
 \item{ln_ml}{the natural log of the marginal likelihood.}
}
\description{
This function takes the log posterior probability of the data (posterior likelihood) for the various cases, the log prior case probabilities, and
returns log posterior odds.
}
\examples{
posterior_summary()
}
