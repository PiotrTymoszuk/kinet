% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OOP.R
\name{plot.kinetic}
\alias{plot.kinetic}
\title{Plot a kinetic class object.}
\usage{
\method{plot}{kinetic}(
  x,
  type = c("fit", "frequency", "diagnostic"),
  cust_theme = ggplot2::theme_classic(),
  plot_title = NULL,
  plot_subtitle = NULL,
  plot_tag = NULL,
  ...
)
}
\arguments{
\item{x}{a kinetic class object.}

\item{type}{plot type, 'fit' by default.}

\item{cust_theme}{a ggplot theme.}

\item{plot_title}{plot title.}

\item{plot_subtitle}{plot subtitle.}

\item{plot_tag}{plot tag. If NULL, number or range of complete observations
is presented.}

\item{...}{extra arguments passed to \code{\link{plot_kin_numeric_}}
(fit plots), \code{\link{plot_kin_fct_}} (frequency plots) or
\code{\link{get_qc_plots}} (diagnostic plots).}
}
\value{
a ggplot object or a list of ggplot objects.
}
\description{
Generates a plot of the actual numeric values (type = 'fit'),
a plot with the frequency (type = 'frequency', for binomial-family objects)
or a series of diagnostic plots of model residuals (type = 'diagnostic').
}
