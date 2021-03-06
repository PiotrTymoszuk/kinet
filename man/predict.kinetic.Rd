% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OOP.R
\name{predict.kinetic}
\alias{predict.kinetic}
\title{Predict values for a kinetic object.}
\usage{
\method{predict}{kinetic}(object, ...)
}
\arguments{
\item{object}{a kinetic class object.}

\item{...}{extra arguments passed to
\code{\link[lme4]{predict.merMod}}.}
}
\value{
a vector with the predicted values.
}
\description{
Prediction of the fitted of new values for a finetic object.
}
