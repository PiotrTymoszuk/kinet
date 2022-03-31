# Non-exported utils.

# Modeling utils ------

#' Generate an identity matrix.
#'
#' @description Creates an identity matrix of the given dimension.
#' @param k matrix dimension.
#' @return a matrix.

  identity_matrix_ <- function(k) {

    out_matrix <- matrix(0, k, k)

    for (i in 1:k) {

      out_matrix[i, i] <- 1

    }

    out_matrix

  }

#' Retrieve inference stats or ANOVA from a merMod object.
#'
#' @description Obtains inference statistics or ANOVA table for a mixed-effect
#' model (class merMod).
#' @details a wrapper around
#' \code{\link[lmerTest]{contest.lmerModLmerTest}} function.
#' @param model a merMod mixed-effect model object.
#' @param L a contrast vector or a matrix. If NULL, inference stats for all
#' model estimates are returned.
#' @param joint logical, make an F-test of potentially several contrast vectors?
#' If FALSE single DF t-tests are applied to each vector or each row of
#' contrasts matrices.
#' @param ... additional arguments passed to
#' \code{\link[lmerTest]{contest.lmerModLmerTest}}.
#' @return a data frame with model estimates, confidence intervals and their
#' significance or, if joint = TRUE, an ANOVA table.

  get_lme_results_ <- function(model, L = NULL, joint = FALSE, ...){

    stopifnot(class(model) %in% c('lmerModLmerTest', 'lmerMod'))

    if(is.null(L)) {

      L <- kinet:::identity_matrix_(length(model@beta))

    }

    if(joint) {

      return(lmerTest:::contest(model = model,
                                L = L,
                                joint = joint))

    } else {

      coef_names <- rownames(summary(model)$coefficients)

      output <- lmerTest:::contest.lmerModLmerTest(model = model,
                                                   L = L,
                                                   joint = joint, ...)

      output <- dplyr::mutate(output,
                              parameter = coef_names)

      output <- output[c('parameter',
                         names(output)[names(output)!= 'parameter'])]

      output <- tibble::as_tibble(rlang::set_names(output,
                                                   c('parameter',
                                                     'estimate',
                                                     'se',
                                                     'df',
                                                     't',
                                                     'lower_ci',
                                                     'upper_ci',
                                                     'p_value')))

      output <- dplyr::mutate(output,
                              n = nrow(model.frame(model)))

      return(tibble::tibble(output[c('parameter',
                                     'estimate',
                                     'se',
                                     'df',
                                     't',
                                     'lower_ci',
                                     'upper_ci',
                                     'n',
                                     'p_value')]))

    }

  }

#' Inference statistic for a glmerMod model.
#'
#' @description Calculates model estimate inference stats for a GLM-like model,
#' more specifically for the glmerMod class.
#' @details In case, the default confidence interval calculation method fails,
#' falls back to the computation based on the normal distribution of the
#' estimate values.
#' @param model a glmerMod mixed-effect model object.
#' @param exponentiate logical, should the model extimates be returned in
#' an exponent form?
#' @param ... extra arguments passed to \code{\link[lme4]{summary.merMod}}.
#' @return a data frame with model estimates, confidence intervals and their
#' significance.

  get_glm_results_ <- function(model,
                               exponentiate = FALSE, ...) {

    ## entry control

    stopifnot(any(class(model) == 'glmerMod'))
    stopifnot(is.logical(exponentiate))

    ## model coefficients

    n_number <- nrow(model.frame(model))

    coefs <- data.frame(summary(model, ...)$coefficients)

    coefs <- rlang::set_names(coefs, c('estimate', 'se', 'z', 'p_value'))

    coefs <- tibble::rownames_to_column(coefs, 'parameter')

    coefs <- dplyr::mutate(coefs,
                           n = nrow(model.frame(model)))

    ## confidence intervals

    ci <- try(tibble::aas_tibble(confint(model)),
              silent = TRUE)

    if(any(class(ci) == 'try-error')) {

      ci <- tibble::tibble(lower_ci = coefs[['estimate']] + coefs[['se']] * qnorm(0.025),
                           upper_ci = coefs[['estimate']] + coefs[['se']] * qnorm(0.975))

      warning('Falling back to CI calculation from normal distribution',
              call. = FALSE)

    } else {

      ci <- ci[-1, ] ## dropping out the CI estimates for sigma

      ci <- rlang::set_names(c('lower_ci', 'upper_ci'))

    }

    ## common inference table

    if(exponentiate) {

      coefs <- dplyr::mutate(coefs,
                             estimate = exp(estimate),
                             se = exp(se))

      ci <- dplyr::mutate(ci,
                          lower_ci = exp(lower_ci),
                          upper_ci = exp(upper_ci))

    }

    coefs <- cbind(coefs, ci)

    tibble::tibble(coefs[c('parameter',
                           'estimate',
                           'se',
                           'z',
                           'lower_ci',
                           'upper_ci',
                           'n',
                           'p_value')])


  }


# Plotting utils: numeric feature plot -----

#' Generate a trajectory plot for a kinetic object.
#'
#' @description Plots a numeric kinetic object. Points represent single
#' observations, medians and IQR (interquartile ranges) at particular time
#' points are displayed as lines and colored regions.
#' @return a ggplot object.
#' @param kinetic_object a kinetic class object.
#' @param point_color color of the points (observations).
#' @param point_alpha point alpha.
#' @param outcome_color color of the median and IQR representations.
#' @param fitted_color color of the model fit line.
#' @param fitted_line style of the model fit line.
#' @param cust_theme a ggplot theme.
#' @param jitter_w_perc horizontal jittering of the point, in percents.
#' @param jitter_h_perc vertical jittering of the points, in percents.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param IQR logical, display the IQR of the outcome?
#' @param fitted logical, display the model fit?
#' @param IQR_fitted logical, display the model fit IQR?

  plot_kin_numeric_ <- function(kinetic_object,
                                point_color = 'gray60',
                                point_alpha = 0.5,
                                outcome_color = 'steelblue',
                                fitted_color = 'coral3',
                                fitted_line = 'dashed',
                                cust_theme = ggplot2::theme_classic(),
                                jitter_w_perc = 2,
                                jitter_h_perc = 1,
                                plot_title = NULL,
                                plot_subtitle = NULL,
                                plot_tag = NULL,
                                x_lab = kinetic_object$time,
                                y_lab = kinetic_object$response,
                                IQR = TRUE,
                                fitted = TRUE,
                                IQR_fitted = FALSE) {

    ## entry control

    stopifnot(kinet::is_kinetic(kinetic_object))

    stopifnot(any(class(cust_theme) == 'theme'))
    stopifnot(jitter_w_perc <= 100)
    stopifnot(jitter_w_perc >= 0)
    stopifnot(jitter_h_perc <= 100)
    stopifnot(jitter_h_perc >= 0)

    stopifnot(is.logical(IQR))
    stopifnot(is.logical(fitted))
    stopifnot(is.logical(IQR_fitted))

    ## plotting tables

    plot_tbl <- residuals(kinetic_object)

    ## handling the categorized responses/binomial modeling

    if(is.factor(plot_tbl$.outcome)) {

      plot_tbl <- dplyr::mutate(plot_tbl, .outcome = as.numeric(.outcome) - 1)

    }

    ## setting up the point jitters and the tag

    jitter_w = jitter_w_perc * diff(range(plot_tbl[[kinetic_object$time]],
                                          na.rm = TRUE))/100

    jitter_h = jitter_h_perc * diff(range(plot_tbl$.outcome,
                                          na.rm = TRUE))/100

    if(is.null(plot_tag)) {

      n_numbers <- dplyr::count(plot_tbl,
                                .data[[kinetic_object$time]])$n

      n_numbers <- range(n_numbers)

      if(diff(n_numbers) != 0) {

        plot_tag <- paste0('\nn = ',
                           n_numbers[1],
                           ' - ',
                           n_numbers[2])

      } else {

        plot_tag <- paste('\nn =',
                          n_numbers[1])

      }

    }

    ## summary tables with the outcome and fitted medians/IQR

    summary_real <- dplyr::group_by(plot_tbl, .data[[kinetic_object$time]])

    summary_real <- dplyr::summarise(summary_real,
                                     median = median(.outcome),
                                     perc25 = quantile(.outcome, 0.25, na.rm = TRUE),
                                     perc75 = quantile(.outcome, 0.75, na.rm = TRUE),
                                     type = 'outcome')

    if(!fitted) {

      summary_tbl <- summary_real

    } else {

      summary_fitted <- dplyr::group_by(plot_tbl, .data[[kinetic_object$time]])

      summary_fitted <- dplyr::summarise(summary_fitted,
                                         median = median(.fitted),
                                         perc25 = quantile(.fitted, 0.25, na.rm = TRUE),
                                         perc75 = quantile(.fitted, 0.75, na.rm = TRUE),
                                         type = 'fitted')

      summary_tbl <- rbind(summary_real,
                           summary_fitted)

    }

    ## kinetic plot

    kinet_plot <- ggplot2::ggplot(plot_tbl,
                                  ggplot2:::aes(x = .data[[kinetic_object$time]],
                                                y = .outcome)) +
      ggplot2::geom_line(data = summary_tbl,
                         ggplot2::aes(y = median,
                                      color = type,
                                      linetype = type))

    if(IQR) {

      kinet_plot <- kinet_plot +
        ggplot2::geom_ribbon(data = dplyr::filter(summary_tbl,
                                                  type == 'outcome'),
                             ggplot2::aes(y = median,
                                          ymin = perc25,
                                          ymax = perc75,
                                          fill = type),
                             alpha = 0.1,
                             color = 'gray80')

    }

    if(all(IQR_fitted, fitted)) {

      kinet_plot <- kinet_plot +
        ggplot2::geom_ribbon(data = dplyr::filter(summary_tbl,
                                                  type == 'fitted'),
                             ggplot2::aes(y = median,
                                          ymin = perc25,
                                          ymax = perc75,
                                          fill = type),
                             alpha = 0.1,
                             color = 'gray80')

    }

    kinet_plot +
      ggplot2::geom_point(shape = 21,
                          size = 2,
                          alpha = point_alpha,
                          fill = point_color,
                          position = ggplot2::position_jitter(width = jitter_w,
                                                              height = jitter_h)) +
      ggplot2::scale_color_manual(values = c(outcome = outcome_color,
                                             fitted = fitted_color),
                                  labels = c(outcome = 'Outcome',
                                             fitted = 'Fitted'),
                                  name = '') +
      ggplot2::scale_fill_manual(values = c(outcome = outcome_color,
                                            fitted = fitted_color),
                                 labels = c(outcome = 'Outcome',
                                            fitted = 'Fitted'),
                                 name = '') +
      ggplot2::scale_linetype_manual(values = c(outcome = 'solid',
                                                fitted = fitted_line),
                                     labels = c(outcome = 'Outcome',
                                                fitted = 'Fitted'),
                                     name = '') +
      cust_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = 'gray90'),
                     plot.tag.position = 'bottom') +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = x_lab,
                    y = y_lab)

  }

# Plotting utils: prevalence plot -------

#' Plot frequency of a binary feature.
#'
#' @description Plots the frequency of a binary feature together with the
#' model predictions.
#' @return a ggplot object.
#' @param kinetic_object a kinetic class object.
#' @param point_color color of the points and text labels.
#' @param outcome_color color of the median and IQR representations.
#' @param fitted_color color of the model fit line.
#' @param fitted_line style of the model fit line.
#' @param cust_theme a ggplot theme.
#' @param plot_title plot title.
#' @param plot_subtitle plot subtitle.
#' @param plot_tag plot tag.
#' @param x_lab x axis title.
#' @param y_lab y axis title.
#' @param labels logical, display percent labels for the time points?
#' @param repel_labels logical, place the labels to reduce their overlap?
#' @param signif_digits significant digits for the percent value rounding.
#' @param label_size size of the percent labels.
#' @param fitted logical, display the model fit?

  plot_kin_fct_ <- function(kinetic_object,
                            point_color = 'steelblue',
                            outcome_color = 'steelblue',
                            fitted_color = 'coral3',
                            fitted_line = 'dashed',
                            cust_theme = ggplot2::theme_classic(),
                            plot_title = NULL,
                            plot_subtitle = NULL,
                            plot_tag = NULL,
                            x_lab = kinetic_object$time,
                            y_lab = '% cohort',
                            labels = TRUE,
                            repel_labels = FALSE,
                            signif_digits = 2,
                            label_size = 2.3,
                            fitted = FALSE) {

    ## entry control

    stopifnot(kinet::is_kinetic(kinetic_object))
    stopifnot(any(class(cust_theme) == 'theme'))
    stopifnot(is.logical(labels))
    stopifnot(is.logical(repel_labels))
    stopifnot(is.logical(fitted))

    signif_digits <- as.integer(signif_digits)

    ## plotting tables

    plot_tbl <- residuals(kinetic_object)

    ## handling the categorized responses/binomial modeling

    if(is.factor(plot_tbl$.outcome)) {

      plot_tbl <- dplyr::mutate(plot_tbl,
                                .outcome = as.numeric(.outcome) - 1)

    }

    ## setting up the plot tag

    if(is.null(plot_tag)) {

      n_numbers <- dplyr::count(plot_tbl,
                                .data[[kinetic_object$time]])$n

      n_numbers <- range(n_numbers)

      if(diff(n_numbers) != 0) {

        plot_tag <- paste0('\nn = ',
                           n_numbers[1],
                           ' - ',
                           n_numbers[2])

      } else {

        plot_tag <- paste('\nn =',
                          n_numbers[1])

      }

    }

    ## summary tables with the outcome prevalence in the subsequent time points

    summary_real <- dplyr::group_by(plot_tbl, .data[[kinetic_object$time]])

    summary_real <- dplyr::summarise(summary_real,
                                     prevalence = mean(.outcome, na.rm = TRUE),
                                     type = 'outcome')

    if(!fitted) {

      summary_tbl <- summary_real

    } else {

      summary_fitted <- dplyr::group_by(plot_tbl, .data[[kinetic_object$time]])

      summary_fitted <- dplyr::summarise(summary_real,
                                       prevalence = mean(.fitted, na.rm = TRUE),
                                       type = 'fitted')

      summary_tbl <- rbind(summary_real,
                           summary_fitted)

    }

    summary_tbl <- dplyr::mutate(summary_tbl,
                                 point_lab = ifelse(type == 'outcome',
                                                    signif(prevalence * 100, signif_digits),
                                                    NA))

    ## kinetic plot

    kinet_plot <- ggplot2::ggplot(summary_tbl,
                                  ggplot2::aes(x = .data[[kinetic_object$time]],
                                               y = prevalence * 100)) +
      ggplot2::geom_line(ggplot2::aes(color = type,
                                      linetype = type)) +
      ggplot2::geom_point(data = dplyr::filter(summary_tbl,
                                               type == 'outcome'),
                          shape = 21,
                          size = 2,
                          fill = point_color) +
      ggplot2::scale_color_manual(values = c(outcome = outcome_color,
                                             fitted = fitted_color),
                                  labels = c(outcome_color = 'Actual',
                                             fitted = 'Fitted'),
                                  name = '') +
      ggplot2::scale_fill_manual(values = c(outcome_color = outcome_color,
                                            fitted = fitted_color),
                                 labels = c(outcome = 'Actual',
                                            fitted = 'Fitted'),
                                 name = '') +
      ggplot2::scale_linetype_manual(values = c(outcome = 'solid',
                                                fitted = fitted_line),
                                     labels = c(outcome = 'Actual',
                                                fitted = 'Fitted'),
                                     name = '') +
      cust_theme +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = 'gray90'),
                     plot.tag.position = 'bottom') +
      ggplot2::labs(title = plot_title,
                    subtitle = plot_subtitle,
                    tag = plot_tag,
                    x = x_lab,
                    y = y_lab)

    if(labels) {

      if(!repel_labels) {

        kinet_plot <- kinet_plot +
          ggplot2::geom_text(ggplot2::aes(label = point_lab),
                             size = label_size,
                             hjust = 0.25,
                             vjust = -1,
                             color = point_color)

      } else {

        kinet_plot <- kinet_plot +
          ggplot2::geom_text_repel(ggplot2::aes(label = point_lab),
                                   size = label_size,
                                   hjust = 0.25,
                                   vjust = -1,
                                   color = point_color,
                                   box.padding = 0.1,
                                   direction = 'x')

      }

    }

    return(kinet_plot)

  }







# Residual plots -----

#' Plots of model residuals.
#'
#' @description Generates a series of residual plots: residuals vs fitted,
#' standardized residuals vs fitted, squared residuals vs fitted and residuals
#' quantile-quantile plot.
#' @param kinet_object a kinetic object.
#' @param cust_theme a ggplot theme.
#' @param resid.type type of the residuals, passed to
#' \code{\link{residuals.kinetic}}.
#' @param point_wjitter horizontal point jittering.
#' @param point_hjitter vertical point jittering.
#' @param point_alpha plot point alpha.
#' @param ... extra arguments passed to \code{\link{residuals.kinetic}}.

  get_qc_plots <- function(kinetic_object,
                           cust_theme = ggplot2::theme_classic(),
                           resid.type = 'response',
                           point_wjitter = 0.01,
                           point_hjitter = 0.01,
                           point_alpha = 0.75, ...) {


    ## entry control

    stopifnot(kinet::is_kinetic(kinetic_object))
    stopifnot(any(class(cust_theme) == 'theme'))

    ## plotting table

    plot_tbl <- residuals(kinetic_object, type = resid.type, ...)

    ## residual plots

    resid_plots <- list(x = c('.fitted',
                              '.fitted',
                              '.fitted',
                              '.std.resid'),
                        y = c('.resid',
                              '.std.resid',
                              '.sq.std.resid',
                              '.expect.norm'),
                        plot_title = c('Residuals vs. fitted',
                                       'Standardized residuals vs. fitted',
                                       'Sqared residuals vs. fitted',
                                       'QQ standardized residuals vs expected normal'))

    resid_plots <- purrr::pmap(resid_plots,
                               function(x, y, plot_title) ggplot2::ggplot(plot_tbl,
                                                                          ggplot2::aes(x = .data[[x]],
                                                                                       y = .data[[y]],
                                                                                       fill = .candidate_missfit)) +
                                 ggplot2::labs(title = plot_title))

    resid_plots <- purrr::map(resid_plots,
                              ~.x +
                                ggplot2::geom_point(size = 2,
                                                    shape = 21,
                                                    position = ggplot2::position_jitter(width = point_wjitter,
                                                                                        height = point_hjitter)))

    resid_plots[1:3] <- purrr::map(resid_plots[1:3],
                                   ~.x +
                                     ggplot2::geom_smooth(method = 'loess',
                                                          color = 'black',
                                                          fill = 'dodgerblue2'))

    resid_plots[[4]] <- resid_plots[[4]] +
      ggplot2::geom_smooth(method = 'lm',
                           color = 'black',
                           fill = 'dodgerblue2')

    resid_plots <- purrr::map(resid_plots,
                              ~.x +
                                cust_theme +
                                ggplot2::scale_fill_manual(values = c(no = 'cornflowerblue',
                                                                      yes = 'firebrick4'),
                                                           name = 'Candidate outlier'))

    rlang::set_names(resid_plots,
                     c('resid_fitted',
                       'std.resid_fitted',
                       'sq.resid_fitted',
                       'qq.std.resid'))

  }
