% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference_lrt.R
\name{normality}
\alias{normality}
\title{Test for normality of the model residuals.}
\usage{
normality(kinetic_object, ...)
}
\arguments{
\item{kinetic_object}{a kinetic class object.}

\item{...}{additional arguments passed to \code{\link{residuals.kinetic}}.}
}
\value{
a tibble with the test results: test statistic and p value.
}
\description{
Tests the normality of the model resioduals with Shapiro-Wilk
test.
}
