# Kinetic class and it's S3 OOP.

# Kinetic class constructor ------

#' Make a kinetic class object.
#'
#' @description Generates an instance of the 'kinetic' class given a merMod
#' model and it's base features.
#' @param model a merMod model.
#' @param response the name of the dependent variable.
#' @param ID the name of the mixed-effect variable, i.e. specifying the
#' data pairing/blocking.
#' @param time the name of the time variable.
#' @param family the model family name (GLM).
#' @param order model order.
#' @return an instance of the 'kinetic' class.
#' @export

  kinetic <- function(model,
                      response,
                      ID,
                      time,
                      family,
                      order) {

    ## entry control

    if(!class(model) %in% c('lmerModLmerTest', 'lmerMod', 'glmerMod')) {

      stop('Please provide a valid lmerMod or lmerModLmerTest class model.',
           call. = FALSE)

    }

    if(order > 0) {

      if(any(!c(response, ID, time) %in% names(model.frame(model)))) {

        stop('ID and time variable must be present in the modeling data.',
             call. = FALSE)

      }

    }

    structure(list(model = model,
                   response = response,
                   ID = ID,
                   time = time,
                   family = family,
                   order = order),
              class = 'kinetic')


  }
# Class testing -----

#' Check for the kinetic_object class.
#'
#' @description Checks if the object is an instance of the 'kinetic' class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_kinetic <- function(x) {

    any(class(x) == 'kinetic')

  }

# Extraction methods -------

#' Model frame of the kinetic object.
#'
#' @description Extracts the kinetic model data.
#' @param formula a kinetic class object.
#' @param ... extra arguments, currently none.
#' @return a data frame.
#' @export

  model.frame.kinetic <- function(formula, ...) {

    stopifnot(kinet::is_kinetic(formula))

    model.frame(formula$model)

  }

#' Formula of a kinetic object.
#'
#' @description Extracts the model formula of a kinetic object.
#' @param x a kinetic class object.
#' @param ... extra arguments, currently none.
#' @return a formula object.
#' @export

  formula.kinetic <- function(x, ...) {

    stopifnot(kinet::is_kinetic(x))

    formula(x$model)

  }

#' Get observation numbers for a kinetic object.
#'
#' @description Extracts the number of complete observations..
#' @param object a kinetic class object.
#' @param ... extra arguments, currently none.
#' @export

  nobs.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    nrow(model.frame(object$model))

  }

# Model residiuals and predictions ------

#' Residuals of a kinetic object.
#'
#' @description Extracts model residuals from a kinetic object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to
#' \code{\link[lme4]{residuals.merMod}}.
#' @return a tibble with the actual and predicted values and the residuals.
#' @details Potential missfits are identified with the +/- 2*SD rule.
#' @export

  residuals.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    res_tbl <- tibble::as_tibble(model.frame(object$model))

    res_tbl <- dplyr::mutate(res_tbl,
                  .observation = 1:nrow(res_tbl),
                  .outcome = res_tbl[[1]],
                  .fitted = predict(object$model, type = 'response'),
                  .resid = residuals(object$model, ...),
                  .std.resid = scale(.resid)[, 1],
                  .sq.std.resid = .std.resid^2,
                  .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 'yes', 'no'))

    ## expected normal

    res_tbl <- res_tbl[order(res_tbl$.std.resid), ]

    res_tbl$.expect.norm <- qnorm(ppoints(nrow(res_tbl)))

    res_tbl

  }

#' Predict values for a kinetic object.
#'
#' @description Prediction of the fitted of new values for a finetic object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to
#' \code{\link[lme4]{predict.merMod}}.
#' @return a vector with the predicted values.
#' @export

  predict.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    predict(object$model, ...)

  }

# Printing -----

#' Print a kinetic object.
#'
#' @description Prints a kinetic object.
#' @param x a kinetic class object.
#' @param ... extra arguments, currently none.
#' @return none, called for it's side effects.
#' @export

  print.kinetic <- function(x, ...) {

    stopifnot(kinet::is_kinetic(x))

    print(x$model)

  }

# Model features -----

#' AIC for a kinetic object.
#'
#' @description Retrieves the value of the Akaike Information Criterion for
#' a kinetic class object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to \code{\link[stats]{AIC}}.
#' @return a numeric value.
#' @export

  AIC.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    AIC(object$model, ...)

  }

#' BIC for a kinetic object.
#'
#' @description Retrieves the value of the Bayesian Information Criterion for
#' a kinetic class object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to \code{\link[stats]{BIC}}.
#' @return a numeric value.
#' @export

  BIC.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    BIC(object$model, ...)

  }

#' Deviance of a kinetic object.
#'
#' @description Calculates deviance for a kinetic class object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to \code{\link[lme4]{deviance.merMod}}.
#' @return a numeric value.
#' @export

  deviance.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    deviance(object$model, ...)

  }

# Inference -------

#' Model coefficients and inference stats for a kinetic class object.
#'
#' @description Calculates coefficient values, their errors, 95% confidence
#' intervals (CI) and p values, as appropriate for the modeling family.
#' @details In case of a CI calculation failure, the function is falling back
#' to computation based on the normal distribution of the estimate value.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to the downstream functions:
#' \code{\link[lmerTest]{contest.lmerModLmerTest}} (gaussian family) or
#' \code{\link[lme4]{summary.merMod}}.
#' @return a data frame with the values of model estimates, errors, CI, numbers
#' of complete observations and p values for the fixed-term estimates.
#' @export

  coef.kinetic <- function(object, ...) {

    ## entry control

    stopifnot(kinet::is_kinetic(object))

    ## inference stats

    if(object$family == 'gaussian') {

      inf_stats <- kinet:::get_lme_results_(object$model, ...)

    } else {

      inf_stats <- kinet:::get_glm_results_(object$model, ...)

    }

    ## model info

    coefs <- dplyr::mutate(inf_stats,
                           response = object$response,
                           time = object$time,
                           ID = object$ID,
                           family = object$family,
                           order = 0:object$order)

    ## output

    coefs[c('response',
            'time',
            'ID',
            'family',
            'order',
            names(inf_stats))]

  }

# ANOVA -------

#' ANOVA for a kinetic class object.
#'
#' @description Obtains an ANOVA table for a kinetic object.
#' @param object a kinetic class object.
#' @param ... extra arguments passed to \code{\link[lme4]{anova.merMod}}.
#' @return a data frame with the sum of squares, F statistic and the fraction
#' of explained variance.
#' @export

  anova.kinetic <- function(object, ...) {

    stopifnot(kinet::is_kinetic(object))

    aov_tbl <- as.data.frame(anova(object$model))

    aov_tbl <- dplyr::mutate(aov_tbl,
                             frac_explained = `Sum Sq`/sum(`Sum Sq`))

    tibble::as_tibble(tibble::rownames_to_column(aov_tbl, 'variable'))

  }

# Summary -----

#' Model summary for a kinetic class object.
#'
#' @description Retrieves inference or fit statistics or normality test results
#' for the model residuals from a kinetic class object.
#' @param object a kinetic class object.
#' @param type type of the statistic returned.
#' @param resid.type type of the residuals, passed to
#' \code{\link{residuals.kinetic}}.
#' @param ... extra arguments passed to \code{\link{coef.kinetic}}
#' or \code{\link{residuals.kinetic}}.
#' @export summary.kinetic
#' @export

  summary.kinetic <- function(object,
                              type = c('inference', 'fit', 'assumptions'),
                              resid.type = 'response',
                              ...) {

    stopifnot(kinet::is_kinetic(object))

    type <- match.arg(type, c('inference', 'fit', 'assumptions'))

    switch(type,
           inference = coef(object, ...),
           fit = kinet:::get_stats(object, type = resid.type, ...),
           assumptions = kinet:::normality(object, type = resid.type, ...))

  }

# Plotting -----

#' Plot a kinetic class object.
#'
#' @description Generates a plot of the actual numeric values (type = 'fit'),
#' a plot with the frequency (type = 'frequency', for binomial-family objects)
#' or a series of diagnostic plots of model residuals (type = 'diagnostic').
#' @return a ggplot object or a list of ggplot objects.
#' @param x a kinetic class object.
#' @param type plot type, 'fit' by default.
#' @param cust_theme a ggplot theme.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag. If NULL, number or range of complete observations
#' is presented.
#' @param ... extra arguments passed to \code{\link{plot_kin_numeric_}}
#' (fit plots), \code{\link{plot_kin_fct_}} (frequency plots) or
#' \code{\link{get_qc_plots}} (diagnostic plots).
#' @export plot.kinetic
#' @export

  plot.kinetic <- function(x,
                           type = c('fit', 'frequency', 'diagnostic'),
                           cust_theme = ggplot2::theme_classic(),
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           plot_tag = NULL,
                           ...) {

    ## entry control

    stopifnot(kinet::is_kinetic(x))

    type <- match.arg(type[1], c('fit', 'frequency', 'diagnostic'))

    ## plotting

    switch(type,
           fit = kinet:::plot_kin_numeric_(kinetic_object = x,
                                           cust_theme = cust_theme,
                                           plot_title = plot_title,
                                           plot_subtitle = plot_subtitle,
                                           plot_tag = plot_tag, ...),
           frequency = kinet:::plot_kin_fct_(kinetic_object = x,
                                             cust_theme = cust_theme,
                                             plot_title = plot_title,
                                             plot_subtitle = plot_subtitle,
                                             plot_tag = plot_tag, ...),
           diagnostic = kinet:::get_qc_plots(kinetic_object = x,
                                             cust_theme = cust_theme, ...))

  }

# END -----
