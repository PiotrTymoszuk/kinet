% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling.R
\name{model_kinetic}
\alias{model_kinetic}
\title{Generate a kinetic object model.}
\usage{
model_kinetic(
  data,
  response,
  time = "time",
  ID = "ID",
  family = "gaussian",
  order = 2,
  ...
)
}
\arguments{
\item{data}{a data frame.}

\item{response}{name of the dependent variable.}

\item{time}{name of the time variable.}

\item{ID}{name of the pairing/grouping or block variable.}

\item{family}{modeling family.}

\item{...}{extra arguments passed to \code{\link[lme4]{lmer}} or
\code{\link[lme4]{glmer}}.}

\item{model}{order.}
}
\value{
a kinetic class object.
}
\description{
Generates a mixed-effect n-order model and embeds it
in a kinetic class object.
}
\details{
Technically, a wrapper around \code{\link[lme4]{lmer}}
(if family is 'gaussian') and \code{\link[lme4]{glmer}} functions.
}
